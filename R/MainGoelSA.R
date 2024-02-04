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

MainGoelSA<-function(data,sites,iter=100000,bar=90000){
  #Valores iniciais
  theta=3.525654
  delta=0.01
  f=1/ 366.8474
  W=as.matrix(rep(9,ncol(data)))
  b=1
  v=1

  X=cbind(as.matrix(rep(1,ncol(data))),sites)
  Z=cbind(as.matrix(rep(1,ncol(data))),sites)
  Fj=Z
  Psi=as.matrix(rep(0,ncol(X)))

  BB1=diag(100,ncol(Z))
  AA1=as.matrix(rep(0,ncol(Z)))
  A=as.matrix(rep(0,ncol(Z)))
  B=diag(100,ncol(Z))

  Psi=as.matrix(rep(0,ncol(X)))
  lgama=as.matrix(c(-3.17022118,0.06602419,1.33314203))#as.matrix(rep(0.01,ncol(Z)))
  leta=as.matrix(c(1.384896960,0.032667601,0.004153729))#as.matrix(rep(1,ncol(Z)))

  n=ncol(data)
  m=nrow(data)
  tempdados=is.na(data)
  nj=m-apply(tempdados,2,sum)

  Tt=array(NA,dim=c(1,n))
  for(y in 1:n){
    Tt[1,y]=data[nj[y],y]
  }

  aa1=2.01
  bb1=1.005

  c3=(-2*log(0.05)/max(dist(sites)))*0.1
  d3=0.1

  SU1=0.00001
  SU2=0.0001
  SU3=0.0001
  SU4=100
  MW=NULL
  MWT=NULL
  Meta=NULL
  MetaT=NULL
  Mgamma=NULL
  MgammaT=NULL
  MPsi=NULL
  Mv=NULL
  Mb=NULL
  MbT=NULL
  Mf=NULL
  MfT=NULL
  MdeltaT=NULL
  Mdelta=NULL
  MthetaT=NULL
  Mtheta=NULL

  for(j in 1:iter){

    if(j<=bar){
      temp=amostrarWsaGOEL(lgama,leta,delta,f,theta,W,Z,Fj,data,t(Tt),b,v,sites,SU1,X,Psi)
      W=temp[[1]]
      MWT=c(MWT,temp[[2]])

      temp=amostrardeltaGOEL(lgama,leta,delta,f,theta,W,Z,Fj,data,t(Tt),0.01,100)
      delta=temp[[1]]
      MdeltaT=c(MdeltaT,temp[[2]])

      temp=amostraretaGOELSA(lgama,leta,delta,f,theta,W,Z,Fj,data,t(Tt),AA1,BB1,sites,SU2)
      leta=temp[[1]]
      MetaT=c(MetaT,temp[[2]])

      temp=amostrargamaGOELSA(lgama,leta,delta,f,theta,W,Z,Fj,data,t(Tt),AA1,BB1,sites,SU3)
      lgama=temp[[1]]
      MgammaT=c(MgammaT,temp[[2]])

      varPsi=solve(solve(B)+t(X)%*%solve(gSigma(b,v,sites))%*%X)
      medPsi=(t(A)%*%solve(B)+t(W)%*%solve(gSigma(b,v,sites))%*%X)%*%varPsi
      Psi=as.matrix(MASS::mvrnorm(1,medPsi,varPsi))

      RR=gCorr(b,sites)
      aa=(n/2)+aa1
      bb=0.5*t(W-X%*%Psi)%*%solve(RR)%*%(W-X%*%Psi)+bb1

      v=1/rgamma(1,shape=aa, rate = bb)

      temp=amostrarb(W,v,b,sites,c3,d3,X,Psi,SU4)
      b=temp[[1]]
      MbT=c(MbT,temp[[2]])

      temp=amostrarfGOEL(lgama,leta,delta,f,theta,W,Z,Fj,data,t(Tt),0.00001/2,1/(365+10),1/(365-10))
      f=temp[[1]]
      MfT=c(MfT,temp[[2]])

      temp=amostrarthetaGOEL(lgama,leta,delta,f,theta,W,Z,Fj,data,t(Tt),0.05)
      theta=temp[[1]]
      MthetaT=c(MthetaT,temp[[2]])

      if((j%%50)==0){
        SU1=sintonizarN(bar,0.20,SU1,MWT,j)
        SU2=sintonizarN(bar,0.20,SU2,MetaT,j)
        SU3=sintonizarN(bar,0.20,SU3,MgammaT,j)
        SU4=sintonizar(bar,0.44,SU4,MbT,j)

      }else{

      }
      print(j)
    }else{

      temp=amostrarWsaGOEL(lgama,leta,delta,f,theta,W,Z,Fj,data,t(Tt),b,v,sites,SU1,X,Psi)
      W=temp[[1]]
      MWT=c(MWT,temp[[2]])
      MW=rbind(MW,t(W))

      temp=amostrardeltaGOEL(lgama,leta,delta,f,theta,W,Z,Fj,data,t(Tt),0.01,100)
      delta=temp[[1]]
      MdeltaT=c(MdeltaT,temp[[2]])
      Mdelta=c(Mdelta,delta)

      temp=amostraretaGOELSA(lgama,leta,delta,f,theta,W,Z,Fj,data,t(Tt),AA1,BB1,sites,SU2)
      leta=temp[[1]]
      MetaT=c(MetaT,temp[[2]])
      Meta=rbind(Meta,t(leta))

      temp=amostrargamaGOELSA(lgama,leta,delta,f,theta,W,Z,Fj,data,t(Tt),AA1,BB1,sites,SU2)
      lgama=temp[[1]]
      MgammaT=c(MgammaT,temp[[2]])
      Mgamma=rbind(Mgamma,t(lgama))

      varPsi=solve(solve(B)+t(X)%*%solve(gSigma(b,v,sites))%*%X)
      medPsi=(t(A)%*%solve(B)+t(W)%*%solve(gSigma(b,v,sites))%*%X)%*%varPsi
      Psi=as.matrix(MASS::mvrnorm(1,medPsi,varPsi))
      MPsi=rbind(MPsi,t(Psi))

      RR=gCorr(b,sites)
      aa=(n/2)+aa1
      bb=0.5*t(W-X%*%Psi)%*%solve(RR)%*%(W-X%*%Psi)+bb1

      v=1/rgamma(1,shape=aa, rate = bb)
      Mv=c(Mv,v)

      temp=amostrarb(W,v,b,sites,c3,d3,X,Psi,SU4)
      b=temp[[1]]
      MbT=c(MbT,temp[[2]])
      Mb=c(Mb,b)

      temp=amostrarfGOEL(lgama,leta,delta,f,theta,W,Z,Fj,data,t(Tt),0.00001/2,1/(365+10),1/(365-10))
      f=temp[[1]]
      MfT=c(MfT,temp[[2]])
      Mf=c(Mf,f)

      temp=amostrarthetaGOEL(lgama,leta,delta,f,theta,W,Z,Fj,data,t(Tt),0.05)
      theta=temp[[1]]
      MthetaT=c(MthetaT,temp[[2]])
      Mtheta=c(Mtheta,theta)


      print(j)
    }

  }
  return(list(MW,MWT,Mdelta,MdeltaT,Mgamma,MgammaT,MPsi,Mv,Mb,MbT,Mf,MfT,Mtheta,MthetaT))
}
