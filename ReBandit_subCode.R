######################################
######################################
##### updating variance estimates ######
update_variance<-function(mean_r,mean_r2,size,M_r){
  ## element in mean_r: reward mean of each arm
  ## element in mean_r2: mean of reward^2 of each arm
  ## element in size: size of each arm
    N<-max(sum(size),1);
    N1<-max(sum(size[size>=1]-1),1);
    sigmar2<-sum(size*mean_r2-size*mean_r^2)/N1;
    if(N1==0){sigmar2<-0}
    u<-mean_r-M_r;
    Nstar<-max(N-sum(size^2)/N,1);
    m<-sum(size*u^2)/Nstar;
    sigmamu2<-max(m,0);
    if(N==0){sigmamu2<-0;sigmar2<-0}
  ###
  output=list("sigmamu2"=sigmamu2,"sigmar2"=sigmar2)
  return(output)
}


####
blup<-function(sigmamu2,sigmar2,mean_r,size,M_r,Known=F,mu0=NULL){
  ####### blup estimator
  ## Known: priors are known or unknown
  ## mu0: set mu_0 (need to estimate it by default)
  K<-length(size);
    if((sigmamu2==0)&&(sigmar2==0)){
      w<-rep(0,K);
    } else {
      w<-size*sigmamu2/(size*sigmamu2+sigmar2);
    }
    w[size==0]<-0;
    AA<-(1-w)*size;
    sum_1ws<-max(sum((1-w)*size),1);
    if(Known==F){
      if(sum_1ws==0){bar<-M_r;
      } else {
      bar<-sum(mean_r*AA)/sum_1ws;
      }
    } else {bar<-mu0;}
    muhat<-w*mean_r+(1-w)*bar;
    if(Known==F){
      if(sum_1ws==0){
        tau2<-max(sigmar2,0)*w/size;
      } else {
      tau2<-max(sigmar2,0)*(w/size+(1-w)^2/sum_1ws);
      }
    } else {tau2<-w/size*sigmar2;}
  ###
  output=list("muhat"=muhat,"tau2"=tau2)
  return(output)
}

#########
arms<-function(theta,B=50,method="BLUP",N=10^5,mu0=NULL,Known=F,DEG=1,MODEL="Gaussian",GPrior=NULL,BerINIT=c(1,1),sigmaR=0.5,KnownMu0=F){
  K<-length(theta);
  regretb<-matrix(0,N,B);
  ##
  for(b in 1:B){
    alpha<-rep(BerINIT[1],K);beta<-rep(BerINIT[2],K);
    timek<-rep(0,K);mean_r<-rep(0,K);mean_r2<-rep(0,K);
    sigmamu2<-0;sigmar2<-0;
    M_r<-0;
    for(t in 1:N){
      thetahat<-rep(10000,K);
       #####
      if(method=="BLUP"){ ### BLUP estimates
         size<-timek;
         update.variance<-update_variance(mean_r,mean_r2,size,M_r);
         sigmamu2<-update.variance$sigmamu2;
         sigmar2<-update.variance$sigmar2;
         estimation<-blup(sigmamu2,sigmar2,mean_r,size,M_r);
         muhat<-estimation$muhat;
         tau2<-estimation$tau2;
         bound<-(tau2)^0.5;
      } else if(method=="BLUP-Known"){ ### BLUP with everything is known
        size<-timek;
        mu0<-GPrior[1];
        sigmamu2<-GPrior[2];
        sigmar2<-GPrior[3];
        if(KnownMu0==F){
          estimation<-blup(sigmamu2,sigmar2,mean_r,size,M_r);
        } else {
          estimation<-blup(sigmamu2,sigmar2,mean_r,size,M_r,mu0,Known=T);
        }
        muhat<-estimation$muhat;
        tau2<-estimation$tau2;
        bound<-(tau2)^0.5;
       }
    ####   
    for(k in 1:K){
       #################### begin 
        if(method=="BLUP"){ #ReUCB
            if(timek[k]<=0){
              thetahat[k]<-10000;
            } else {thetahat[k]<-muhat[k]+bound[k]*(DEG*log(max(t-1,1)))^0.5;}
        } else if(method=="BLUP-Known"){ #ReUCB*
            if(timek[k]<=0){
              thetahat[k]<-10000;
            } else {thetahat[k]<-muhat[k]+bound[k]*(DEG*log(max(t-1,1)))^0.5;}
        } else if(method=="UCB1"){
           thek<-mean_r[k];
           thetahat[k]<-thek+DEG*(8*sigmaR^2*log(max(t-1,1))/max(timek[k],1))^0.5;
           if(timek[k]==0){thetahat[k]<-10000;}
        } else if(method=="G-TS"){ #TS for Gaussian
            mu0<-GPrior[1];
            sigma02<-GPrior[2];
            sigma2<-GPrior[3];
            nk<-timek[k];
            mupost<-1/(1/sigma02+nk/sigma2)*(mu0/sigma02+nk*mean_r[k]/sigma2);
            sdpost<-1/(1/sigma02+nk/sigma2)^0.5;
            thetahat[k]<-rnorm(1,mean=mupost,sd=DEG^0.5*sdpost);
        } else if(method=="B-TS"){ #TS for Bernoulli
            thetahat[k]<-rbeta(1,shape1=alpha[k],shape2=beta[k]);
        } else if(method=="G-BUCB"){ #BayesUCB for Gaussian
            mu0<-GPrior[1];
            sigma02<-GPrior[2];
            sigma2<-GPrior[3];
            nk<-timek[k];
            mupost<-1/(1/sigma02+nk/sigma2)*(mu0/sigma02+nk*mean_r[k]/sigma2);
            sdpost<-1/(1/sigma02+nk/sigma2)^0.5;
            AAA<-rnorm(100,mean=mupost,sd=sdpost);
            if(t<100){
            thetahat[k]<-quantile(AAA,1-1/max(t,1));
            } else {thetahat[k]<-max(AAA[-which.max(AAA)])}
        } else if(method=="B-BUCB"){ #BayesUCB for Bernoulli
            AAA<-rbeta(100,shape1=alpha[k],shape2=beta[k])
            if(t<100){
            thetahat[k]<-quantile(AAA,1-1/max(t,1));
            } else {thetahat[k]<-max(AAA[-which.max(AAA)])}
        }
       }
       ################## end
       K2<-c(1:K)[thetahat==max(thetahat)]
       if(length(K2)>1){
        xidx<-sample(K2,1);
       } else {xidx<-which.max(thetahat)}
      #####
      ###########
      if(MODEL=="Bernoulli"){
        r<-rbinom(n=1,1,prob=theta[xidx]);
      } else if (MODEL=="Unif") {
        rx<-runif(1,0,1);
        thetax<-theta[xidx];
        r<-sum(rx<thetax)
      } else if (MODEL=="Gaussian"){
        r<-rnorm(n=1,theta[xidx],sd=sigmaR);
      } else if (MODEL=="TGaussian"){
        r<-rnorm(n=1,theta[xidx],sd=sigmaR);
        if(r<0){r<-0;}
        if(r>1){r<-1;}
      } else if (MODEL=="Beta4"){
        thetax<-rbeta(n=1,shape1=4*theta[xidx],shape2=4*(1-theta[xidx]));
        r<-thetax;
      }
      M_r<-M_r+1/t*(r-M_r);
      alpha[xidx]<-alpha[xidx]+r;
      beta[xidx]<-beta[xidx]+1-r;
      timek[xidx]<-timek[xidx]+1;
      mean_r[xidx]<-mean_r[xidx]+1/timek[xidx]*(r-mean_r[xidx]);
      mean_r2[xidx]<-mean_r2[xidx]+1/timek[xidx]*(r^2-mean_r2[xidx]);
      ## regret
      regretb[t,b]<-max(theta)-theta[xidx];
    }
    print(b);
  }
  ###
  M2<-matrix(0,N,B);
  for(b in 1:B){
    M2[,b]<-cumsum(regretb[,b])
  }
  avg<-c(M2%*%rep(1/B,B));
  regretT<-M2[N,];
  ###
  output=list("avg"=avg,"regretT"=regretT)
  return(output)
}



############################################
##### Binary bandits (One Example) #########
############################################
B<-1000;
N<-10^4;
theta<-runif(50,0.2,0.5);##
GPrior<-c(0.35,0.3^2/12,0.5^2)
## ReUCB
arms_regret<-arms(theta=theta,method="BLUP",N=N,B=B,DEG=1,MODEL="Bernoulli");
regretT.RE<-arms_regret$regretT;avg.RE<-arms_regret$avg;
## ReUCB*
arms_regret<-arms(theta=theta,method="BLUP-Known",N=N,B=B,DEG=1,MODEL="Bernoulli",GPrior=GPrior,KnownMu0=F);
regretT.Known<-arms_regret$regretT;avg.Known<-arms_regret$avg;
## UCB1
arms_regret<-arms(theta=theta,method="UCB1",N=N,B=B,DEG=1,MODEL="Bernoulli");
regretT.UCB1<-arms_regret$regretT;avg.UCB1<-arms_regret$avg;
## UCB1 with a=1
arms_regret<-arms(theta=theta,method="UCB1",N=N,B=B,DEG=(1/8)^0.5,MODEL="Bernoulli");
regretT.UCB0<-arms_regret$regretT;avg.UCB0<-arms_regret$avg;
## BayesUCB
arms_regret<-arms(theta=theta,method="B-BUCB",N=N,B=B,DEG=1,MODEL="Bernoulli",BerINIT=c(1,1));
regretT.BUCB<-arms_regret$regretT;avg.BUCB<-arms_regret$avg;
## TS
arms_regret<-arms(theta=theta,method="B-TS",N=N,B=B,DEG=1,MODEL="Bernoulli",BerINIT=c(1,1));
regretT.TS<-arms_regret$regretT;avg.TS<-arms_regret$avg;



