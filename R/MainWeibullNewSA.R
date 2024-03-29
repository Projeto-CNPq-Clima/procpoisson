
#'Algoritmo MCMC Poisson Process
#'
#' Esta é a função para realizar o algoritmo Monte Carlo via Cadeias de Markov para os parâmetros de um processo de Poisson.
#'
#' @param data Matriz de dados simulados de um processo de Poisson
#' @param sites Malha de pontos
#' @param iter número inteiro positivo, interações do algoritmo MCMC
#' @param bar número inteiro positivo, período de burn-in do algoritmo MCMC
#'
#' @export
#'

MainWeibullNewSA<-function(data,sites,iter=1000,bar=900){
  #Valores iniciais
  bw=1
  vw=1
  bm=1
  vm=1
  theta=pi/2
  delta=0.001
  f=1/365
  M=as.matrix(rep(log(0.89),ncol(data)))
  W=as.matrix(rep(0,ncol(data)))

  #Hiperparametros
  aa1=0.001
  bb1=0.001

  aa2=2.01
  bb2=1.005

  c3=(-2*log(0.05)/max(dist(sites)))*0.1
  d3=0.1


  SU1=0.000001
  SU2=0.000001
  SU3=100
  SU4=100

  X=cbind(as.matrix(rep(1,ncol(data))),sites)
  Z=cbind(as.matrix(rep(1,ncol(data))),sites)

  Psi=as.matrix(rep(0,ncol(X)))
  A1=as.matrix(rep(0,ncol(X)))
  B1=diag(100,ncol(X))

  Beta=as.matrix(rep(0,ncol(Z)))
  A=as.matrix(rep(0,ncol(Z)))
  B=diag(100,ncol(Z))



  n=ncol(data)
  m=nrow(data)
  tempdados=is.na(data)
  nj=m-apply(tempdados,2,sum)

  Tt=array(NA,dim=c(1,n))
  for(y in 1:n){
    Tt[1,y]=data[nj[y],y]
  }


  MMj=NULL
  MMT=NULL

  MW=NULL
  MWT=NULL

  MPsi=NULL
  MBeta=NULL

  Mvw=NULL
  Mbw=NULL
  MbwT=NULL

  Mvm=NULL
  Mbm=NULL
  MbmT=NULL

  Mdelta=NULL
  MdeltaT=NULL

  Mtheta=NULL
  MthetaT=NULL

  Mf=NULL
  MfT=NULL

  for(j in 1:iter){


    if(j<=bar){
      temp=amostrarWsa(delta,theta,W,M,sites,X,Psi,bw,vw,nj,Tt,data,SU1,f)
      W=as.matrix(temp[[1]])
      MWT=c(MWT,temp[[2]])

      temp=amostrarMsa(delta,theta,W,M,sites,X,Beta,bm,vm,nj,Tt,data,SU2,f)
      M=as.matrix(temp[[1]])
      MMT=c(MMT,temp[[2]])

      temp=amostrardelta(theta,delta,W,M,data,nj,Tt,0.01,f,100)
      delta=temp[[1]]
      MdeltaT=c(MdeltaT,temp[[2]])

      temp=amostrartheta(theta,delta,W,M,data,nj,Tt,0.05,f)
      theta=temp[[1]]
      MthetaT=c(MthetaT,temp[[2]])

      temp=amostrarf(theta,delta,W,M,data,nj,Tt,0.00001/2,1/(365+10),1/(365-10),f)
      f=temp[[1]]
      MfT=c(MfT,temp[[2]])

      RR=gCorr(bw,sites)
      aa=(n/2)+aa2
      bb=0.5*t(W-X%*%Psi)%*%solve(RR)%*%(W-X%*%Psi)+bb2

      vw=1/rgamma(1,shape=aa, rate = bb)

      RR1=gCorr(bm,sites)
      aa=(n/2)+aa2
      bb=0.5*t(M-Z%*%Beta)%*%solve(RR1)%*%(M-Z%*%Beta)+bb2

      vm=1/rgamma(1,shape=aa, rate = bb)

      varPsi=solve(solve(B)+t(X)%*%solve(gSigma(bw,vw,sites))%*%X)
      medPsi=(t(A)%*%solve(B)+t(W)%*%solve(gSigma(bw,vw,sites))%*%X)%*%varPsi
      Psi=as.matrix(MASS:mvrnorm(1,medPsi,varPsi))


      varBeta=solve(solve(B)+t(Z)%*%solve(gSigma(bm,vm,sites))%*%Z)
      medBeta=(t(A)%*%solve(B)+t(M)%*%solve(gSigma(bm,vm,sites))%*%Z)%*%varBeta
      Beta=as.matrix(MASS:mvrnorm(1,medBeta,varBeta))

      temp=amostrarb(W,vw,bw,sites,c3,d3,X,Psi,SU3)
      bw=temp[[1]]
      MbwT=c(MbwT,temp[[2]])

      temp=amostrarb(M,vm,bm,sites,c3,d3,Z,Beta,SU4)
      bm=temp[[1]]
      MbmT=c(MbmT,temp[[2]])


      if((j%%50)==0){
        SU1=sintonizarN(bar,0.15,SU1,MWT,j)
        SU2=sintonizarN(bar,0.15,SU2,MMT,j)
        SU3=sintonizar(bar,0.44,SU3,MbwT,j)
        SU4=sintonizar(bar,0.44,SU4,MbmT,j)

      }else{

      }





      print(j)

    }else{

      temp=amostrarWsa(delta,theta,W,M,sites,X,Psi,bw,vw,nj,Tt,data,SU1,f)
      W=as.matrix(temp[[1]])
      MW=rbind(MW,t(W))
      MWT=c(MWT,temp[[2]])

      temp=amostrarMsa(delta,theta,W,M,sites,X,Beta,bm,vm,nj,Tt,data,SU2,f)
      M=as.matrix(temp[[1]])
      MMj=rbind(MMj,t(M))
      MMT=c(MMT,temp[[2]])

      temp=amostrardelta(theta,delta,W,M,data,nj,Tt,0.01,f,100)
      delta=temp[[1]]
      Mdelta=c(Mdelta,delta)
      MdeltaT=c(MdeltaT,temp[[2]])

      temp=amostrartheta(theta,delta,W,M,data,nj,Tt,0.05,f)
      theta=temp[[1]]
      Mtheta=c(Mtheta,theta)
      MthetaT=c(MthetaT,temp[[2]])

      temp=amostrarf(theta,delta,W,M,data,nj,Tt,0.00001/2,1/(365+10),1/(365-10),f)
      f=temp[[1]]
      Mf=c(Mf,f)
      MfT=c(MfT,temp[[2]])

      RR=gCorr(bw,sites)
      aa=(n/2)+aa2
      bb=0.5*t(W-X%*%Psi)%*%solve(RR)%*%(W-X%*%Psi)+bb2

      vw=1/rgamma(1,shape=aa, rate = bb)
      Mvw=c(Mvw,vw)

      RR1=gCorr(bm,sites)
      aa=(n/2)+aa2
      bb=0.5*t(M-Z%*%Beta)%*%solve(RR1)%*%(M-Z%*%Beta)+bb2

      vm=1/rgamma(1,shape=aa, rate = bb)
      Mvm=c(Mvm,vm)

      varPsi=solve(solve(B)+t(X)%*%solve(gSigma(bw,vw,sites))%*%X)
      medPsi=(t(A)%*%solve(B)+t(W)%*%solve(gSigma(bw,vw,sites))%*%X)%*%varPsi
      Psi=as.matrix(MASS:mvrnorm(1,medPsi,varPsi))
      MPsi=rbind(MPsi,t(Psi))

      varBeta=solve(solve(B)+t(Z)%*%solve(gSigma(bm,vm,sites))%*%Z)
      medBeta=(t(A)%*%solve(B)+t(M)%*%solve(gSigma(bm,vm,sites))%*%Z)%*%varBeta
      Beta=as.matrix(MASS:mvrnorm(1,medBeta,varBeta))
      MBeta=rbind(MBeta,t(Beta))

      temp=amostrarb(W,vw,bw,sites,c3,d3,X,Psi,SU3)
      bw=temp[[1]]
      Mbw=c(Mbw,bw)
      MbwT=c(MbwT,temp[[2]])

      temp=amostrarb(M,vm,bm,sites,c3,d3,Z,Beta,SU4)
      bm=temp[[1]]
      Mbm=c(Mbm,bm)
      MbmT=c(MbmT,temp[[2]])



      print(j)

    }

  }
  return(list(W,MW,MWT,M,MMj,MMT,Mvw,Mvm,Beta,MbwT,MbmT))
}
