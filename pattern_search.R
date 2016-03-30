# Three-dimensional parameter fitting of a siminf model using Pattern Search optimization
# Using observed data to compute the goal function.
# P.Bauer, 2015

#loads some kind of fitting data
data(region)

# max. number of retries if error occur
maxerrs=10;
# number of siminf trajectories
ntr=80;
#grow ratio
growratio=1.5;

#parallel settings
#Sys.setenv(GOMP_CPU_AFFINITY = "0-32");
#Sys.setenv(OMP_NUM_THREADS = "12");

#generate goal function
compres=TRUE;
if(compres) {
  #generate one trajectory of the model at goal parameter values (this could be given by data, too)
  ntr_res=1;
  
  vtec  = vtec_model(init      = init,
                      tspan     = tspan,
                      events    = events,
                      initial_infectious_pressure = initial_infectious_pressure,
                      response_calves = response_calves,
                      response_young_stock = response_young_stock,
                      response_adults = response_adults,
                      recover_calves = recover_calves,
                      recover_young_stock = recover_young_stock,
                      recover_adults = recover_adults,
                      alpha = ntr_res, #this field is modified to tell the solver how many trajectories we want to compute.
                      beta_q1 = beta_q1,
                      beta_q2 = beta_q2,
                      beta_q3 = beta_q3,
                      beta_q4 = beta_q4,
                      epsilon = epsilon)

presult <- siminf(vtec, solver="psiminf", nthreads=12, report_level=0, use_rmatio = FALSE,
link_matio=FALSE,inputfile='input.mat',outfile='output.mat')
res=presult@U/ntr_res;
}

#perturbation of parameters
perturbation=runif(1, min = 0, max = 5); 
k1=response_calves*perturbation;
k2=response_young_stock*perturbation;
k3=response_adults*perturbation;

#stoping criteria:
# 1) tolerance
tol=1000;
# 2) max. no. of iterations
max_iter=200;

#iteration count
iter=0;

#initial step-length
delta1=k1*0.1;
delta2=k2*0.1;
delta3=k3*0.1;

#lists used for progress logging
total=prev_total=Inf;
totals=c(NULL,NULL,NULL,NULL,NULL,NULL);
totals_hist=NULL;
delta_hist=NULL;
k_hist=NULL;
errors_hist=NULL;
coretime_hist=NULL;

#store initial values
k_hist=rbind(k_hist,c(k1,k2,k3));
delta_hist=rbind(delta_hist,c(delta1,delta2,delta3));

#growing flag
grow=0;
#0: didnt grow
#1: just grew
#2: finished growing

#start optimization loop
while(total>tol || iter<max_iter) {
  
  #parameters for this iteration: k_1,2,3 +/- delta
  found=0;
  for(i in 1:6) {
    if(i==1) {
      k1_new=k1+delta1;
      k2_new=k2;
      k3_new=k3;
    }
    if(i==2) {
      k1_new=k1-delta1;
      k2_new=k2;
      k3_new=k3;
    }
    if(i==3) {
      k1_new=k1;
      k2_new=k2+delta2;
      k3_new=k3;
    }
    if(i==4) {
      k1_new=k1;
      k2_new=k2-delta2;
      k3_new=k3;
    }
    if(i==5) {
      k1_new=k1;
      k2_new=k2;
      k3_new=k3+delta3;
    }
    if(i==6) {
      k1_new=k1;
      k2_new=k2;
      k3_new=k3-delta3;
    }
    
    if(k1_new>0.0 && k2_new>0.0 & k3_new>0.0) {
    #create vtec moel
    vtec  = vtec_model(init      = init,
                     tspan     = tspan,
                     events    = events,
                     initial_infectious_pressure = initial_infectious_pressure,
                     response_calves = k1_new,
                     response_young_stock = k2_new,
                     response_adults = k3_new,
                     recover_calves = recover_calves,
                     recover_young_stock = recover_young_stock,
                     recover_adults = recover_adults,
                     alpha = ntr, #seed
                     beta_q1 = beta_q1,
                     beta_q2 = beta_q2,
                     beta_q3 = beta_q3,
                     beta_q4 = beta_q4,
                     epsilon = epsilon)
    
    #execute using exception handling (for some weird parameters, solver can give negative values)
    errcnt=0;
    repeat {
     feedb=try(newres <- siminf(vtec, solver="psiminf", nthreads=12, report_level=0, use_rmatio = FALSE,
	   link_matio=FALSE,inputfile='/home/pavpa354/tests/input.mat',outfile='/home/pavpa354/tests/output.mat')); 
  
      if ( inherits( feedb,  "try-error")){
        print("Error");
        errcnt=errcnt+1;
        
        #if error happened, perturbe parameters slightly
        k1_new=k1_new*1.001;
        k2_new=k2_new*1.001;
        k3_new=k3_new*1.001;
        vtec  = vtec_model(init      = init,
                     tspan     = tspan,
                     events    = events,
                     initial_infectious_pressure = initial_infectious_pressure,
                     response_calves = k1_new,
                     response_young_stock = k2_new,
                     response_adults = k3_new,
                     recover_calves = recover_calves,
                     recover_young_stock = recover_young_stock,
                     recover_adults = recover_adults,
                     alpha = ntr, #seed
                     beta_q1 = beta_q1,
                     beta_q2 = beta_q2,
                     beta_q3 = beta_q3,
                     beta_q4 = beta_q4,
                     epsilon = epsilon)
        if(errcnt>=maxerrs) {
          stop("Maximum amount of error retries reached.");
        }
      } else {
        break;
      }
    }
    
    #compute average
    avg=newres@U/ntr; #E[X]
    
    #compute std.dev using E[X^2] stored in Uvar-vector (availble from modified solver file)
    x2=Uvar/ntr; #E[X^2]
    std=sqrt(x2-avg^2)/sqrt(ntr); #E[X^2] - E[X]^2
    
    #log error
    errors_hist=rbind(errors_hist,sum(std)/sqrt(ntr)); # std.err.in.mean
    
    #log coretime
    coretime_hist=rbind(coretime_hist,coretime);
    
    #compute the goal function, respecting counties
    rss=c(0,0,0,0,0,0);
    resc1=res[seq(1, 223326, 6),];
    resc2=res[seq(2, 223326, 6),];
    resc3=res[seq(3, 223326, 6),];
    resc4=res[seq(4, 223326, 6),];
    resc5=res[seq(5, 223326, 6),];
    resc6=res[seq(6, 223326, 6),];
    
    nresc1=avg[seq(1, 223326, 6),];
    nresc2=avg[seq(2, 223326, 6),];
    nresc3=avg[seq(3, 223326, 6),];
    nresc4=avg[seq(4, 223326, 6),];
    nresc5=avg[seq(5, 223326, 6),];
    nresc6=avg[seq(6, 223326, 6),];  
    
    fun=sum;
    for(j in 1:21) {
      ids=which(region$region==j,arr.ind=TRUE);
      
      rss[1]=rss[1]+fun((nresc1[ids,]-resc1[ids,])^2);
      rss[2]=rss[2]+fun((nresc2[ids,]-resc2[ids,])^2);
      rss[3]=rss[3]+fun((nresc3[ids,]-resc3[ids,])^2);
      rss[4]=rss[4]+fun((nresc4[ids,]-resc4[ids,])^2);
      rss[5]=rss[5]+fun((nresc5[ids,]-resc5[ids,])^2);
      rss[6]=rss[6]+fun((nresc6[ids,]-resc6[ids,])^2);
    }
    totals[i]=sum(rss);
    
    } else {
      #if k_1,2,3<0
      totals[i]=Inf;
    }
  }
  totals_hist=rbind(totals_hist,totals);
  
  i=which.min(totals);
  total=min(totals);
  
    #found better value -> move the new parameters there
    if(total<prev_total) {
      if(i==1) {
        k1=k1+delta1;
        k2=k2;
        k3=k3;
      }
      if(i==2) {
        k1=k1-delta1;
        k2=k2;
        k3=k3;
      }
      if(i==3) {
        k1=k1;
        k2=k2+delta2;
        k3=k3;
      }
      if(i==4) {
        k1=k1;
        k2=k2-delta2;
        k3=k3;
      }
      if(i==5) {
        k1=k1;
        k2=k2;
        k3=k3+delta3;
      }
      if(i==6) {
        k1=k1;
        k2=k2;
        k3=k3-delta3;
      }
      #report and log
      cat(sprintf("Iteration %d, RSS: %f<%f MOVING (i=%d).\n", iter, total,prev_total,i));
      cat(sprintf("k1: %f k2: %f k3: %f.\n",k1,k2,k3));
      k_hist=rbind(k_hist,c(k1,k2,k3));
      delta_hist=rbind(delta_hist,c(delta1,delta2,delta3));
      prev_total=total;
      
      #moving due to a growing region
      if(grow==1) {
        #successfull grow, back to former delta
        delta1=delta1/growratio;
        delta2=delta2/growratio;
        delta3=delta3/growratio;
        grow=0;
      } else if (grow==2){
        #move after shrink, reset the posibility to grow
        grow=0;
      }
    } else {
      #before shrinking, move 2*delta to detect if we're in a local minimum
      if(grow==0) {
        #grow!
        delta1=delta1*growratio;
        delta2=delta2*growratio;
        delta3=delta3*growratio;   
        cat(sprintf("Iteration %d, RSS: %f>%f GROWING.\n", iter,total,prev_total));
        grow=1;
      } else if(grow==2) {
        #regular srhink
        delta1=delta1/2;
        delta2=delta2/2;
        delta3=delta3/2;
        cat(sprintf("Iteration %d, RSS: %f>%f SHRINKING 2x.\n", iter,total,prev_total));
        grow=2;
      } else if(grow==1) {
        #unsuccesfull grow -> shrink
        delta1=delta1/(2*growratio);
        delta2=delta2/(2*growratio);
        delta3=delta3/(2*growratio);
        cat(sprintf("Iteration %d, RSS: %f>%f SHRINKING 4x.\n", iter,total,prev_total));
        grow=2;
      }
          delta_hist=rbind(delta_hist,c(delta1,delta2,delta3));
    }
  iter=iter+1;

  #save all log-files during the execution
  save(totals_hist,delta_hist,k_hist,errors_hist,coretime_hist,file="debug_out.RData")
}
#end of optimization loop

#plot convergence
convergence=apply(totals_hist,1,min);
points=1:length(convergence);
convergence=log(convergence);
df=data.frame(points,convergence)
ggplot(data=df,aes(x=points,y=convergence))+geom_line()+geom_point()

#plot k's
kdim=dim(k_hist);
k1optim=array(response_calves,c(kdim[1],1));
k2optim=array(response_young_stock,c(kdim[1],1));
k3optim=array(response_adults,c(kdim[1],1));

qplot(x=1:kdim[1],y=k_hist[,1],colour="k1",geom="line")+
  geom_line(aes(y = k_hist[,2], colour = "k2"))+
  geom_line(aes(y = k_hist[,3], colour = "k3"))+
  geom_line(aes(y = k1optim, colour = "k1"),linetype="dashed")+
  geom_line(aes(y = k2optim, colour = "k2"),linetype="dashed")+
  geom_line(aes(y = k3optim, colour = "k3"),linetype="dashed")+
  scale_size_area()+xlab("iteration")+ylab("parameter value")+
  ggtitle("k1,k2,k3 convergence against optimum")
