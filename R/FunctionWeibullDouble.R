######################################################
gCorr<-function(b,def){
  n=nrow(def)
  R=exp(-b*(as.matrix(dist(def))))
  mat=R
  mat


}
######################################################
gSigma<-function(b,v,def){
  n<-nrow(def)
  R<-exp(-b*(as.matrix(dist(def))))
  mat<-v*R
  mat


}
######################################################
logvero=function(aalpha,WW,TT,nn,yT){

  res=sum(nn)*log(aalpha)+sum(WW*nn)-sum(exp(WW)*TT^(aalpha))+aalpha*sum(log(yT),na.rm =T)
  res
}
######################################################
logverosa=function(aalpha,WW,ddelta,ttheta,yT,TT,f){

  suma=0
  for(i in 1:ncol(yT)){
    res=sum(log( aalpha*exp(WW[i,])*as.matrix(yT[,i])^(aalpha-1)-ddelta*2*pi*f*sin(2*pi*f*data[,i]+ttheta) ),na.rm =T)
    suma=suma+res
  }

  res1=suma-sum( exp(WW)*t(TT)^(aalpha)+ddelta*cos(2*pi*f*t(TT)+ttheta) )
  res1

}
######################################################
sintonizar=function(bar,taxa,tau,mat,i){

  mat=as.matrix(mat)



  mater=(1/50)*sum(mat[(i-49):i,1])

  if(mater>=taxa){
    delta=min(0.01,(i/50+1)^(-0.5))
    temp4=log(tau)-delta
    temp5=exp(temp4)
    return(temp5)
  }else{
    delta=min(0.01,(i/50+1)^(-0.5))
    temp4=log(tau)+delta
    temp5=exp(temp4)
    return(temp5)
  }



}
######################################################
sintonizarN=function(bar,taxa,tau,mat,i){

  mat=as.matrix(mat)



  mater=(1/50)*sum(mat[(i-49):i,1])

  if(mater>=taxa){
    delta=min(0.01,(i/50+1)^(-0.5))
    temp4=log(tau)+delta
    temp5=exp(temp4)
    return(temp5)
  }else{
    delta=min(0.01,(i/50+1)^(-0.5))
    temp4=log(tau)-delta
    temp5=exp(temp4)
    return(temp5)
  }



}
######################################################
amostrarW=function(WW,MM,S,XX,PPs,bb,vv,nn,TT,ff){

  n=nrow(WW)
  WWprop=as.matrix(mvrnorm(1,WW,ff*diag(1,n)))


  SSig=gSigma(bb,vv,S)

  postWW=sum(WW*nn)-sum(exp(WW)*t(TT)^exp(MM))-0.5*t(WW-XX%*%PPs)%*%solve(SSig)%*%(WW-X%*%PPs)
  postWWprop=sum(WWprop*nn)-sum(exp(WWprop)*t(TT)^exp(MM))-0.5*t(WWprop-XX%*%PPs)%*%solve(SSig)%*%(WWprop-XX%*%PPs)

  prob=min(exp((postWWprop)-(postWW)),1)


  u=runif(1,0,1)

  if(u<prob){

    Wprox=WWprop

    rejei=1


  }else{

    Wprox=WW
    rejei=0
  }





  res=as.matrix(Wprox)
  res=list(Wprox,rejei)
  res

}
######################################################
#amostrarW(W,M,sites,X,Psi,b,v,nj,Tt,SU2)
amostrarM=function(WW,MM,S,ZZ,BBta,bb,vv,data,nn,TT,ff){

  n=nrow(MM)
  MMprop=as.matrix(mvrnorm(1,MM,ff*diag(1,n)))
  logdata=log(data)
  SSig=gSigma(bb,vv,S)

  sum1=0
  sum2=0

  for(j in 1:ncol(data)){
    res1=sum(exp(MM[j,])*logdata[,j],na.rm=T)
    res2=sum(exp(MMprop[j,])*logdata[,j],na.rm=T)
    sum1=sum1+res1
    sum2=sum2+res2
  }

  postMM=sum(MM*nn)-sum(exp(WW)*t(TT)^exp(MM))+sum1-0.5*t(MM-ZZ%*%BBta)%*%solve(SSig)%*%(MM-ZZ%*%BBta)
  postMMprop=sum(MMprop*nn)-sum(exp(WW)*t(TT)^exp(MMprop))+sum2-0.5*t(MMprop-ZZ%*%BBta)%*%solve(SSig)%*%(MMprop-ZZ%*%BBta)

  prob=min(exp((postMMprop)-(postMM)),1)


  u=runif(1,0,1)

  if(u<prob){

    MMprox=MMprop

    rejei=1


  }else{

    MMprox=MM
    rejei=0
  }





  res=as.matrix(MMprox)
  res=list(MMprox,rejei)
  res



}
#amostrarM(W,M,sites,Z,Beta,b,v,data,nj,Tt,SU2)
######################################################
amostrarb=function(W,v,b,loca,ab,bb,X,Psi,u1){

  bprop=rgamma(1,shape=b*u1, rate = u1)

  SSigprop=gSigma(bprop,v,loca)

  if((det(SSigprop)==0)|(bprop< 0.005)){
    return(list(b,0))
  }

  SSig=gSigma(b,v,loca)
  SSigprop=gSigma(bprop,v,loca)


  logp=-0.5*t(W-X%*%Psi)%*%solve(SSig)%*%(W-X%*%Psi)-0.5*log(det(SSig))+(ab-1)*log(b)-bb*b

  logpprop=-0.5*t(W-X%*%Psi)%*%solve(SSigprop)%*%(W-X%*%Psi)-0.5*log(det(SSigprop))+(ab-1)*log(bprop)-bb*bprop

  logprob=logpprop+log(dgamma(b,shape=bprop*u1,rate=u1))-(logp+log(dgamma(bprop,shape=b*u1,rate=u1)))
  prob<-min(c(1,exp(logprob)))

  u=runif(1,0,1)

  if(u<prob){

    bprox=bprop

    rejei=1


  }else{

    bprox=b

    rejei=0;

  }



  res=list(bprox,rejei)
  res
}

######################################################
logverosa=function(MM,WW,ddelta,ttheta,yT,TT,f){

  suma=0
  for(i in 1:ncol(yT)){
    res=sum(log( exp(WW[i,]+MM[i,])*as.matrix(yT[,i])^(exp(MM[i,])-1)-ddelta*2*pi*f*sin(2*pi*f*data[,i]+ttheta) ),na.rm =T)
    suma=suma+res
  }

  res1=suma-sum( exp(WW)*t(TT)^(exp(MM))+ddelta*cos(2*pi*f*t(TT)+ttheta) )
  res1

}
######################################################
amostrarWsa=function(ddelta,ttheta,WW,MM,loca,XX,PPs,bb,vv,nn,TT,yT,ff,f){

  n=nrow(WW)
  WWprop=as.matrix(mvrnorm(1,WW,ff*diag(1,n)))

  tema1=ifelse( ( exp(WWprop+MM)*t(TT)^(exp(MM)-1) )<=(2*pi*f*ddelta),1,0)
  tema2=ifelse( ( exp(WWprop+MM)*(yT[1,])^(exp(MM)-1) )<=(2*pi*f*ddelta),1,0)



  if((sum(tema1)+sum(tema2))>=1){
    return(list(WW,0))
  }else{

  }


  SSig=gSigma(bb,vv,loca)

  postWW=logverosa(MM,WW,ddelta,ttheta,yT,TT,f)-0.5*t(WW-XX%*%PPs)%*%solve(SSig)%*%(WW-XX%*%PPs)
  postWWprop=logverosa(MM,WWprop,ddelta,ttheta,yT,TT,f)-0.5*t(WWprop-XX%*%PPs)%*%solve(SSig)%*%(WWprop-XX%*%PPs)

  prob=min(exp((postWWprop)-(postWW)),1)


  u=runif(1,0,1)

  if(u<prob){

    Wprox=WWprop

    rejei=1


  }else{

    Wprox=WW
    rejei=0
  }





  res=as.matrix(Wprox)
  res=list(Wprox,rejei)
  res

}
#amostrarWsa(delta,theta,W,M,sites,X,Psi,bw,vw,nj,Tt,data,SU1,f)
######################################################
amostrarMsa=function(ddelta,ttheta,WW,MM,loca,XX,PPs,bb,vv,nn,TT,yT,ff,f){

  n=nrow(MM)
  MMprop=as.matrix(mvrnorm(1,MM,ff*diag(1,n)))

  tema1=ifelse( ( exp(WW+MMprop)*t(TT)^(exp(MMprop)-1) )<=(2*pi*f*ddelta),1,0)
  tema2=ifelse( (exp(WW+MMprop)*(yT[1,])^(exp(MMprop)-1) )<=(2*pi*f*ddelta),1,0)


  if((sum(tema1)+sum(tema2))>=1){
    return(list(MM,0))
  }else{

  }


  SSig=gSigma(bb,vv,loca)

  postMM=logverosa(MM,WW,ddelta,ttheta,yT,TT,f)-0.5*t(MM-XX%*%PPs)%*%solve(SSig)%*%(MM-X%*%PPs)
  postMMprop=logverosa(MMprop,WW,ddelta,ttheta,yT,TT,f)-0.5*t(MMprop-XX%*%PPs)%*%solve(SSig)%*%(MMprop-XX%*%PPs)

  prob=min(exp((postMMprop)-(postMM)),1)


  u=runif(1,0,1)

  if(u<prob){

    Mprox=MMprop

    rejei=1


  }else{

    Mprox=MM
    rejei=0
  }





  res=as.matrix(Mprox)
  res=list(Mprox,rejei)
  res

}
#amostrarMsa(delta,theta,W,M,sites,X,Beta,bm,vm,nj,Tt,data,SU2,f)
amostrardelta=function(ttheta,ddelta,WW,MM,yT,nn,TT,u1,f,d){


  deltaprop=runif(1,max(0,ddelta-u1),min(ddelta+u1,d))


  tema1=ifelse( ( exp(WW+MM)*t(TT)^(exp(MM)-1) )<=(2*pi*f*deltaprop),1,0)
  tema2=ifelse( (exp(WW+MM)*(yT[1,])^(exp(MM)-1) )<=(2*pi*f*deltaprop),1,0)


  if((sum(tema1)+sum(tema2))>=1){
    return(list(ddelta,0))
  }else{

  }



  logp=logverosa(MM,WW,ddelta,ttheta,yT,TT,f)-0.5*(log(ddelta)+log(d-ddelta))

  logpprop=logverosa(MM,WW,deltaprop,ttheta,yT,TT,f)-0.5*(log(deltaprop)+log(d-deltaprop))

  logprob=logpprop+log( dunif( ddelta, max(0,deltaprop-u1), min(d,deltaprop+u1) ) )-( logp+log( dunif( deltaprop, max(0,ddelta-u1), min(d,ddelta+u1) ) ) )

  prob<-min(c(1,exp(logprob)))

  u=runif(1,0,1)

  if(u<prob){

    bprox=deltaprop

    rejei=1


  }else{

    bprox=ddelta

    rejei=0;

  }



  res=list(bprox,rejei)
  res

}
#amostrardelta(theta,delta,W,M,data,nj,Tt,0.01,f,100)
########################################################################
amostrartheta=function(ttheta,ddelta,WW,MM,yT,nn,TT,u1,f){


  thetaprop=runif(1,max(0,ttheta-u1),min(theta+u1,2*pi))

  logp=logverosa(MM,WW,ddelta,ttheta,yT,TT,f)-0.5*(log(ttheta)+log(2*pi-ttheta))

  logpprop=logverosa(MM,WW,ddelta,thetaprop,yT,TT,f)-0.5*(log(thetaprop)+log(2*pi-thetaprop))

  logprob=logpprop+log( dunif( ttheta, max(0,thetaprop-u1), min(2*pi,thetaprop+u1) ) )-( logp+log( dunif( thetaprop, max(0,ttheta-u1), min(2*pi,ttheta+u1) ) ) )

  prob<-min(c(1,exp(logprob)))

  u=runif(1,0,1)

  if(u<prob){

    bprox=thetaprop

    rejei=1


  }else{

    bprox=ttheta

    rejei=0;

  }



  res=list(bprox,rejei)
  res

}

#amostrartheta(theta,delta,W,M,data,nj,Tt,0.01,f)
########################################################################
amostrarf=function(ttheta,ddelta,WW,MM,yT,nn,TT,u1,a,b,f){


  fprop=runif(1,max(a,f-u1),min(f+u1,b))

  logp=logverosa(MM,WW,ddelta,ttheta,yT,TT,f)-0.5*(log(f-a)+log(b-f))

  logpprop=logverosa(MM,WW,ddelta,ttheta,yT,TT,fprop)-0.5*(log(fprop-a)+log(b-fprop))

  logprob=logpprop+log( dunif( f, max(a,fprop-u1), min(b,fprop+u1) ) )-( logp+log( dunif( fprop, max(a,f-u1), min(b,f+u1) ) ) )

  prob<-min(c(1,exp(logprob)))

  u=runif(1,0,1)

  if(u<prob){

    bprox=fprop

    rejei=1


  }else{

    bprox=f

    rejei=0;

  }



  res=list(bprox,rejei)
  res

}
#amostrarf(theta,delta,W,M,data,nj,Tt,0.001,1/(365+10),1/(365-10),f)


