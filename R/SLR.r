library(MASS)
library(ibd)

bslr = function(v)
{
	if(v%%2 == 0) {
		if(v > 2) bslr.even(v) else stop("v should be greater than 2")
	} else {
		 if (v>=3) bslr.odd(v) else stop("v should be greater than 2")
		}
}


pbslr = function(v)
{
	if(v%%2 == 0) {
		if(v > 2) pbslr.even(v) else stop("v should be greater than 2")
	} else {
		 if (v>=3) pbslr.odd(v) else stop("v should be greater than 2")
		}
}


bslr.odd = function(v)
{
  n = (v-3)/2
  h<-2*n+3
  p<-(2*n+3)*(n+1)
  k<-2
  nr<-2*n+2
  nc<-2
  r<-h*nr
  mat<-cbind(matrix(1,(v-1)/2,1),matrix(2:((v+1)/2),(v-1)/2,1))
  final=NULL
  for(i in 1:nrow(mat)){
    vec <- mat[i,]
    mat1 <- matrix(nrow = 1, ncol = 2*v)
    mat1[, 1:2] <- vec
    for (i in 3:(2*v)) {
      mat1[, i] <- (mat1[, i-2])%%v + 1
    }
    SLR_1<-matrix(data=NA,nrow=v,ncol=2*v)
    SLR_1[1,]<-mat1
    for(t in 2:(v)){
      for(l in 1:(v-1)){
        SLR_1[t,(2*l-1):(2*l)]<-(SLR_1[t-1,((2*l+1):(2*l+2))])
        SLR_1[t,(2*v-1):(2*v)]<-(SLR_1[t-1,1:2])
      }
    }
    final=cbind(final,SLR_1)
  }
  SLR <-final
  SLR_2 <-matrix(NA,v,v*(v-1)/2)
  for (n in 1:(v*(v-1)/2)) {
    SLR_2[, n] <- paste(SLR[,2*n-1], SLR[,2*n], sep = ",")
  }  
  vect<-c(t(SLR))
  N1<-design_to_N(matrix(vect,nrow = h,byrow = TRUE))
  N2<-nc*matrix(1,v,p)
  N3<-design_to_N(matrix(vect,nrow = h*p,byrow = TRUE))
  X1_prim_X1<-p*nc*diag(1,v)
  X1_prim_X2<-cbind(h*nr*c(rep.int(1,v)), N1, N2, N3)
  X2_prim_X2<-rbind(cbind(h*p*k,p*k*t(c(rep.int(1,h))),h*k*t(c(rep.int(1,p))),k*t(c(rep.int(1,h*p)))),cbind(p*k*c(rep.int(1,h)),p*k*diag(1,h),k*matrix(1,h,p),k*matrix(1,h,h*p)),cbind(h*k*c(rep.int(1,p)),k*matrix(1,p,h),h*k*diag(1,p),k*matrix(1,p,h*p)),cbind(k*c(rep.int(1,h*p)),k*matrix(1,h*p,h),k*matrix(1,h*p,p),k*diag(1,h*p)))
  c_matrix<-X1_prim_X1-X1_prim_X2%*%ginv(X2_prim_X2)%*%t(X1_prim_X2)
  eg=eigen(c_matrix)$values
  cef<-eg[eg>0.0001]/r    
  aef<-1/(mean(1/cef))  
  output = list(design = noquote(SLR_2), Avg.Effi = aef)
  return(output)
}  

bslr.even = function(v)
{
  n = (v-2)/2
  v <- 2*n+2
  h<-n+1
  p<-(2*n+1)*(n+1)
  k<-2
  nr<-2*n+1
  nc<-1
  r<-h*nr
  pairs <- cbind(1:(v-1), rep(v, v-1))
  final=NULL
  for(m in 1:nrow(pairs)){
      new_pairs<-matrix(NA, 1, v)
      new_pairs[,1:2] <- pairs[m, ]
      for (i in 2:(v/2)) {
        for (j in 2:(v/2)) {
          new_pairs[,(2*i-1)] <- (new_pairs[1, 1] + i-2) %% (v-1)+1  
          new_pairs[,(v-(2*j-4))] <- (new_pairs[, v-1] +j-2) %% (v-1) + 1  
        }
      }
      SLR_1<-matrix(NA, v/2, v)
      SLR_1[1,]<-new_pairs
      for(t in 2:(v/2)){
         for(l in 1:(v/2-1)){
           SLR_1[t,(2*l-1):(2*l)]<-(SLR_1[t-1,((2*l+1):(2*l+2))])
           SLR_1[t,(v-1):(v)]<-(SLR_1[t-1,1:2])
          }
      }
      final=cbind(final,SLR_1)
  }
  SLR<- final
  SLR_2 <-matrix(NA,v/2,v*(v-1)/2)
  for (n in 1:(v*(v-1)/2)) {
    SLR_2[, n] <- paste(SLR[,2*n-1], SLR[,2*n], sep = ",")
  }
  vect<-c(t(SLR))
  N1<-design_to_N(matrix(vect,nrow = h,byrow = TRUE))
  N2<-nc*matrix(1,v,p)
  N3<-design_to_N(matrix(vect,nrow = h*p,byrow = TRUE))
  X1_prim_X1<-p*nc*diag(1,v)
  X1_prim_X2<-cbind(h*nr*c(rep.int(1,v)), N1, N2, N3)
  X2_prim_X2<-rbind(cbind(h*p*k,p*k*t(c(rep.int(1,h))),h*k*t(c(rep.int(1,p))),k*t(c(rep.int(1,h*p)))),cbind(p*k*c(rep.int(1,h)),p*k*diag(1,h),k*matrix(1,h,p),k*matrix(1,h,h*p)),cbind(h*k*c(rep.int(1,p)),k*matrix(1,p,h),h*k*diag(1,p),k*matrix(1,p,h*p)),cbind(k*c(rep.int(1,h*p)),k*matrix(1,h*p,h),k*matrix(1,h*p,p),k*diag(1,h*p)))
  c_matrix<-X1_prim_X1-X1_prim_X2%*%ginv(X2_prim_X2)%*%t(X1_prim_X2)
  eg=eigen(c_matrix)$values
  cef<-eg[eg>0.0001]/r    
  aef<-1/(mean(1/cef))  
  output = list(design = noquote(SLR_2), Avg.Effi = aef)
  return(output)
}

pbslr.even = function(v)
{
  n = (v-2)/2
  h<-n+1
  p<-2*n+2
  k<-2
  nr<-2
  nc<-1
  r<-nr*h
  mat<-cbind(matrix(1:(v/2),v/2,1),matrix(v:(v/2+1),v/2,1))
  SLR<- matrix(nrow = nrow(mat), ncol = 2*v)
  SLR[, 1:2]<- mat
  for (i in 3:(2*v)) {
    SLR[, i] <- (SLR[, i-2]) %% v+1
  }
  SLR_2 <- matrix(nrow = nrow(SLR), ncol = ncol(SLR)/2)
  for (i in 1:v) {
    SLR_2[, i] <-paste(SLR[, 2*i-1], SLR[, 2*i], sep=",")
  }
  vect<-c(t(SLR))
  N1<-design_to_N(matrix(vect,nrow = h,byrow = TRUE))
  N2<-nc*matrix(1,v,p)
  N3<-design_to_N(matrix(vect,nrow = h*p,byrow = TRUE))
  X1_prim_X1<-p*nc*diag(1,v)
  X1_prim_X2<-cbind(h*nr*c(rep.int(1,v)), N1, N2, N3)
  X2_prim_X2<-rbind(cbind(h*p*k,p*k*t(c(rep.int(1,h))),h*k*t(c(rep.int(1,p))),k*t(c(rep.int(1,h*p)))),cbind(p*k*c(rep.int(1,h)),p*k*diag(1,h),k*matrix(1,h,p),k*matrix(1,h,h*p)),cbind(h*k*c(rep.int(1,p)),k*matrix(1,p,h),h*k*diag(1,p),k*matrix(1,p,h*p)),cbind(k*c(rep.int(1,h*p)),k*matrix(1,h*p,h),k*matrix(1,h*p,p),k*diag(1,h*p)))
  c_matrix<-X1_prim_X1-X1_prim_X2%*%ginv(X2_prim_X2)%*%t(X1_prim_X2)
  eg=eigen(c_matrix)$values
  cef<-eg[eg>0.0001]/r    
  aef<-1/(mean(1/cef))  
  output = list(design = noquote(SLR_2), Avg.Effi = aef)
  return(output)
}

pbslr.odd = function(v)
{
  n = (v-3)/2
  h<-2*n+3
  p<-2*n+3
  k<-2
  nr<-2
  nc<-2
  r<-nr*h
  pov<-c(1:v)
  random_pair <- sample(pov,2)
  mat <- matrix(nrow = 1, ncol = 2*v)
  mat[, 1:2] <- random_pair
  for (i in 3:(2*v)) {
    mat[, i] <- (mat[, i-2])%%v + 1
  }
  SLR<-matrix(data=NA,nrow=v,ncol=2*v)
  SLR[1,]<-mat[1,]
  for (l in 2:v) {
    for (j in 1:(2*v)) {
      SLR[l,j]=((SLR[1,j] +(l-2))%%v)+1
    }
  }    
  SLR_4 <- matrix(nrow = nrow(SLR), ncol = ncol(SLR)/2)
  for (t in 1:v) {
    SLR_4[, t] <-paste(SLR[, 2*t-1], SLR[, 2*t], sep=",")
  }
  vect<-c(t(SLR))
  N1<-design_to_N(matrix(vect,nrow = h,byrow = TRUE))
  N2<-nc*matrix(1,v,p)
  N3<-design_to_N(matrix(vect,nrow = h*p,byrow = TRUE))
  X1_prim_X1<-p*nc*diag(1,v)
  X1_prim_X2<-cbind(h*nr*c(rep.int(1,v)), N1, N2, N3)
  X2_prim_X2<-rbind(cbind(h*p*k,p*k*t(c(rep.int(1,h))),h*k*t(c(rep.int(1,p))),k*t(c(rep.int(1,h*p)))),cbind(p*k*c(rep.int(1,h)),p*k*diag(1,h),k*matrix(1,h,p),k*matrix(1,h,h*p)),cbind(h*k*c(rep.int(1,p)),k*matrix(1,p,h),h*k*diag(1,p),k*matrix(1,p,h*p)),cbind(k*c(rep.int(1,h*p)),k*matrix(1,h*p,h),k*matrix(1,h*p,p),k*diag(1,h*p)))
  c_matrix<-X1_prim_X1-X1_prim_X2%*%ginv(X2_prim_X2)%*%t(X1_prim_X2)
  eg=eigen(c_matrix)$values
  cef<-eg[eg>0.0001]/r    
  aef<-1/(mean(1/cef))  
  output = list(design = noquote(SLR_4), Avg.Effi = aef)
  return(output)
}

