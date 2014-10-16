sgmconvex<-function(A,B,m,start,iteration){
  # m seeds, assumed 1:m in each graph 
  # solves the convex relaxation of the
  # GM problem ||AD-DB||_F^2 using FW
  require('clue')
  totv<-ncol(A)
  n<-totv-m
  if (m != 0){
    A12<-rbind(A[1:m,(m+1):(m+n)])
    A21<-cbind(A[(m+1):(m+n),1:m])
    B12<-rbind(B[1:m,(m+1):(m+n)])
    B21<-cbind(B[(m+1):(m+n),1:m])
  }
  if (m==0){
    A12<-matrix(0,n,n)
    A21<-matrix(0,n,n)
    B12<-matrix(0,n,n)
    B21<-matrix(0,n,n)}
  A22<-A[(m+1):(m+n),(m+1):(m+n)]
  B22<-B[(m+1):(m+n),(m+1):(m+n)]
  patience<-iteration
  tol<-1
  P<-start
  toggle<-1
  iter<-0
  a1 <- t(A12)%*%A12
  a2 <- t(A22)%*%A22
  b1 <- B21%*%t(B21)
  b2 <- B22%*%t(B22)
  c1 <- A21%*%t(B21)
  c2 <- t(A12)%*%B12
  while (toggle==1 && iter<patience)
  {
    iter<-iter+1
    c3 <- A22%*%P%*%t(B22)
    c4 <- t(A22)%*%P%*%B22
    Grad<- (a1+a2)%*%P+P%*%(b1+b2)-c1-c2-c3-c4;
    mm=max(abs(Grad))
    ind<-matrix(solve_LSAP(Grad+matrix(mm,totv-m,totv-m), maximum =FALSE))
    T<-diag(n)
    T<-T[ind,]
    cc4<-t(A22)%*%T%*%B22
    x<-P%*%t(P)-P%*%t(T)-T%*%t(P)+T%*%t(T)
    y<-t(P)%*%P-t(P)%*%T-t(T)%*%P+t(T)%*%T 
    alpha2<-(a1+a2)%*%x+(b1+b2)%*%y-2*c4%*%t(P)-2*cc4%*%t(T)+2*c4%*%t(T)+2*cc4%*%t(P)
    alpha2<-sum(diag(alpha2))
    xx<-P%*%t(T)+T%*%t(P)-2*T%*%t(T)
    yy<-t(P)%*%T+t(T)%*%P-2*t(T)%*%T
    alpha<-(a1+a2)%*%xx+(b1+b2)%*%yy-2*t(P)%*%(c1+c2)+2*t(T)%*%(c1+c2)+4*cc4%*%t(T)-2*c4%*%t(T)-2*cc4%*%t(P)
    alpha<-sum(diag(alpha))
    f0<-0
    f1<-alpha2+alpha
    if( alpha2==0){
      if (alpha<0){toggle<-0}else{P<-T}
    }else{
    f0<-0
    crit<- -alpha/(2*alpha2)
    f1<- alpha2+alpha
    fcrit<-alpha2*crit^2+alpha*crit
    if(crit < 1 && crit > 0 && fcrit < f0 && fcrit < f1){
      P<- crit*P+(1-crit)*T
    }else if(f0 < f1){
      P<-T
    }else{
      toggle<-0}
  }}
  D<-P
  corr<-matrix(solve_LSAP(P, maximum = TRUE))
  P=diag(n)
  P=rbind(cbind(diag(m),matrix(0,m,n)),cbind(matrix(0,n,m),P[corr,]))
  return(list(P=P, D=D))}