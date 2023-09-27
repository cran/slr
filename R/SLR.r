library(MASS)
library(ibd)

bslr = function(v, k)
{
	if(v%%2 == 0) {
  		if(v > 2) bslr.even(v, k) else stop("v should be greater than 2")
  	} else {
  		 if(v >= 3) bslr.odd(v, k) else stop("v should be greater than 2")
  	}
}

pbslr = function(v, k)
{
  if(v%%2 == 0) {
      if(v > 2) pbslr.even(v, k) else stop("v should be greater than 2")
    } else {
      if (v >= 3) pbslr.odd(v, k) else stop("v should be greater than 2")
    }
}

bslr.odd = function(v, k)
{
  if(!(k ==2 | k == 3)) output = "Facilty is available for k = 2 or 3 only"
  if(k == 2)
  {
      n = (v-3)/2
      h<-2*n+3
      p<-(2*n+3)*(n+1)
      # k<-2
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
  }  
  if(k == 3)
  {
    if(v < 5) stop("v should be an odd integer greater than 4")
    n = (v - 3)/2
    h<-2*n+3
    p<-(n+1)*(2*n+3)
    # k<-3
    nr<-3*(n+1)
    nc<-3
    r<-h*nr
    mat<-cbind(matrix(1,(v-1)/2,1),matrix(2:((v+1)/2),(v-1)/2,1),matrix(v:((v+3)/2),(v-1)/2,1))
    final=NULL
    for(i in 1:nrow(mat)){
      vec <- mat[i,]
      mat1 <- matrix(nrow = 1, ncol = 3*v)
      mat1[, 1:3] <- vec
      for (j in 4:(3*v)) {
        mat1[, j] <- (mat1[, j-3])%%v + 1
      }
      mat1
      SLR_1<-matrix(data=NA,nrow=v,ncol=3*v)
      SLR_1[1,]<-mat1
      for(t in 2:(v)){
        for(l in 1:(v-1)){
          SLR_1[t,(3*l-2):(3*l)]<-(SLR_1[t-1,((3*l+1):(3*l+3))])
          SLR_1[t,(3*v-2):(3*v)]<-(SLR_1[t-1,1:3])
        }
      }
      final=cbind(final,SLR_1)
    }
    SLR <-final
    SLR_2 <-matrix(NA,v,v*(v-1)/2)
    for (x in 1:(v*(v-1)/2)) {
      SLR_2[, x] <- paste(SLR[,3*x-2], SLR[,3*x-1],SLR[,3*x], sep = ",")
    }
    # noquote(SLR_2)
    ## C_matrix
    vect<-c(t(SLR))
    N1<-design_to_N(matrix(vect,nrow = h,byrow = TRUE))
    N2<-nc*matrix(1,v,p)
    N3<-design_to_N(matrix(vect,nrow = h*p,byrow = TRUE))
    X1_prim_X1<-p*nc*diag(1,v)
    X1_prim_X2<-cbind(h*nr*c(rep.int(1,v)), N1, N2, N3)
    X2_prim_X2<-rbind(cbind(h*p*k,p*k*t(c(rep.int(1,h))),h*k*t(c(rep.int(1,p))),k*t(c(rep.int(1,h*p)))),cbind(p*k*c(rep.int(1,h)),p*k*diag(1,h),k*matrix(1,h,p),k*matrix(1,h,h*p)),cbind(h*k*c(rep.int(1,p)),k*matrix(1,p,h),h*k*diag(1,p),k*matrix(1,p,h*p)),cbind(k*c(rep.int(1,h*p)),k*matrix(1,h*p,h),k*matrix(1,h*p,p),k*diag(1,h*p)))
    c_matrix<-X1_prim_X1-X1_prim_X2%*%ginv(X2_prim_X2)%*%t(X1_prim_X2)
    eg=eigen(c_matrix)$values
    cef<-eg[eg>0.0001]/r    ##Canonical efficiency factors
    aef<-1/(mean(1/cef))  ##Average efficiency
    # aef
    output = list(design = noquote(SLR_2), Avg.Effi = aef)
  }  
  return(output)
}  

bslr.even = function(v, k)
{
  if(!(k ==2 | k == 3)) output = "Facilty is available for k = 2 or 3 only"
  if(k == 2)
  {
    n = (v-2)/2
    v <- 2*n+2
    h<-n+1
    p<-(2*n+1)*(n+1)
    # k<-2
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
  }
  if(k == 3)
  {
    # n<-2
    # v<- 2*n+2
    if(v == 2) stop("v should be an even integer greater than 2")
    n = (v-2)/2 
    h<-2*n+2
    p<-(v*(v-1)*(v-2)/6)
    # k<-3
    nr<-((v-1)*(v-2)/2)
    nc<-3
    r<-h*nr
    mat<-combn(v,3)
    Init<-matrix(NA, 1, v*(v-1)*(v-2)/2)
    Init[,1:3]<-mat[,1]
    for (i in 2:(v*(v-1)*(v-2)/6)) {
      Init[,(3*i-2):(3*i)]<-mat[,i]
    }
    SLR<-matrix(NA, v, v*(v-1)*(v-2)/2)
    SLR[1,]<-Init[1,]
    for (j in 2:v) {
      SLR[j,]<-SLR[j-1,]%%v +1
    }
    SLR_1 <-matrix(NA,v,(v*(v-1)*(v-2)/6))
    for (x in 1:(v*(v-1)*(v-2)/6)) {
      SLR_1[, x] <- paste(SLR[,3*x-2], SLR[,3*x-1],SLR[,3*x], sep = ",")
    }
    # noquote(SLR_1)
    ##C_matrix
    vect<-c(t(SLR))
    N1<-design_to_N(matrix(vect,nrow = h,byrow = TRUE))
    N2<-nc*matrix(1,v,p)
    N3<-design_to_N(matrix(vect,nrow = h*p,byrow = TRUE))
    X1_prim_X1<-p*nc*diag(1,v)
    X1_prim_X2<-cbind(h*nr*c(rep.int(1,v)), N1, N2, N3)
    X2_prim_X2<-rbind(cbind(h*p*k,p*k*t(c(rep.int(1,h))),h*k*t(c(rep.int(1,p))),k*t(c(rep.int(1,h*p)))),cbind(p*k*c(rep.int(1,h)),p*k*diag(1,h),k*matrix(1,h,p),k*matrix(1,h,h*p)),cbind(h*k*c(rep.int(1,p)),k*matrix(1,p,h),h*k*diag(1,p),k*matrix(1,p,h*p)),cbind(k*c(rep.int(1,h*p)),k*matrix(1,h*p,h),k*matrix(1,h*p,p),k*diag(1,h*p)))
    c_matrix<-X1_prim_X1-X1_prim_X2%*%ginv(X2_prim_X2)%*%t(X1_prim_X2)
    eg=eigen(c_matrix)$values
    cef<-eg[eg>0.0001]/r    ## Canonical efficiency factors
    aef<-1/(mean(1/cef))  ## Average efficiency
    # aef
    output = list(design = noquote(SLR_1), Avg.Effi = aef)
  }
  return(output)
}

pbslr.even = function(v, k)
{
  if(!(k ==2 | k == 3)) output = "Facilty is available for k = 2 or 3 only"
  if(k == 2)
  {
    n = (v-2)/2
    h<-n+1
    p<-2*n+2
    # k<-2
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
  }
  if(k == 3)
  {
    if(v < 6) stop("v should be an even integer greater than 5")
    n = (v-4)/2
    h<-2*n+4
    p<-(n+1)*(2*n+4)
    # k<-3
    nr<-3*(n+1)
    nc<-3
    r<-h*nr
    mat<-cbind(matrix(1,n+1,1),matrix(2,n+1,1),matrix(3:(n+3),n+1,1))
    final=NULL
    for(i in 1:nrow(mat)){
      vec <- mat[i,]
      mat1 <- matrix(nrow = 1, ncol = 3*v)
      mat1[, 1:3] <- vec
      for (j in 4:(3*v)) {
        mat1[, j] <- (mat1[, j-3])%%v + 1
      }
      mat1
      SLR_1<-matrix(data=NA,nrow=v,ncol=3*v)
      SLR_1[1,]<-mat1
      for(t in 2:(v)){
        for(l in 1:(v-1)){
          SLR_1[t,(3*l-2):(3*l)]<-(SLR_1[t-1,((3*l+1):(3*l+3))])
          SLR_1[t,(3*v-2):(3*v)]<-(SLR_1[t-1,1:3])
        }
      }
      final=cbind(final,SLR_1)
    }
    SLR <-final
    SLR_3 <-matrix(NA,v,(n+1)*v)
    for (x in 1:(v*(n+1))) {
      SLR_3[, x] <- paste(SLR[,3*x-2], SLR[,3*x-1],SLR[,3*x], sep = ",")
    }
    # noquote(SLR_3)
    ## C_matrix
    vect<-c(t(SLR))
    N1<-design_to_N(matrix(vect,nrow = h,byrow = TRUE))
    N2<-nc*matrix(1,v,p)
    N3<-design_to_N(matrix(vect,nrow = h*p,byrow = TRUE))
    X1_prim_X1<-p*nc*diag(1,v)
    X1_prim_X2<-cbind(h*nr*c(rep.int(1,v)), N1, N2, N3)
    X2_prim_X2<-rbind(cbind(h*p*k,p*k*t(c(rep.int(1,h))),h*k*t(c(rep.int(1,p))),k*t(c(rep.int(1,h*p)))),cbind(p*k*c(rep.int(1,h)),p*k*diag(1,h),k*matrix(1,h,p),k*matrix(1,h,h*p)),cbind(h*k*c(rep.int(1,p)),k*matrix(1,p,h),h*k*diag(1,p),k*matrix(1,p,h*p)),cbind(k*c(rep.int(1,h*p)),k*matrix(1,h*p,h),k*matrix(1,h*p,p),k*diag(1,h*p)))
    c_matrix<-X1_prim_X1-X1_prim_X2%*%ginv(X2_prim_X2)%*%t(X1_prim_X2)
    eg=eigen(c_matrix)$values
    cef<-eg[eg>0.0001]/r    ##Canonical efficiency factors
    aef<-1/(mean(1/cef))  ##Average efficiency
    # aef
    output = list(design = noquote(SLR_3), Avg.Effi = aef)
  }
  return(output)
}

pbslr.odd = function(v, k)
{
  if(!(k ==2 | k == 3)) output = "Facilty is available for k = 2 or 3 only"
  if(k == 2)
  {
    n = (v-3)/2
    h<-2*n+3
    p<-2*n+3
    # k<-2
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
  }
  if(k == 3)
  {
    if(v <= 3) stop("v should be an odd integer greater than 4")
    n = (v + 1)/6
    if(n - floor(n) == 0)
    {
      h<-6*n-1
      p<-n*(6*n-1)
      # k<-3
      nr<-3*n
      nc<-3
      r<-h*nr
      mat<-cbind(matrix(1,n,1),matrix((n+1),n,1),matrix((n+2):(2*n+1),n,1))
      final=NULL
      for(i in 1:nrow(mat)){
        vec <- mat[i,]
        mat1 <- matrix(nrow = 1, ncol = 3*v)
        mat1[, 1:3] <- vec
        for (j in 4:(3*v)) {
          mat1[, j] <- (mat1[, j-3])%%v + 1
        }
        mat1
        SLR_1<-matrix(data=NA,nrow=v,ncol=3*v)
        SLR_1[1,]<-mat1
        for(t in 2:(v)){
          for(l in 1:(v-1)){
            SLR_1[t,(3*l-2):(3*l)]<-(SLR_1[t-1,((3*l+1):(3*l+3))])
            SLR_1[t,(3*v-2):(3*v)]<-(SLR_1[t-1,1:3])
          }
        }
        final=cbind(final,SLR_1)
      }
      SLR <-final
      SLR_4 <-matrix(NA,v,n*v)
      for (x in 1:(v*n)) {
        SLR_4[, x] <- paste(SLR[,3*x-2], SLR[,3*x-1],SLR[,3*x], sep = ",")
      }
      # noquote(SLR_41)
    }  
    ##method4_Series2
    
    # n<-4 ##n should be any integer
    # v<-6*n+1
    n = (v - 1)/6
    if(n - floor(n) == 0)
    {
      h<-6*n+1
      p<-n*(6*n+1)
      # k<-3
      nr<-3*n
      nc<-3
      r<-h*nr
      mat<-cbind(matrix(1,n,1),matrix((n+1),n,1),matrix((n+2):(2*n+1),n,1))
      final=NULL
      for(i in 1:nrow(mat)){
        vec <- mat[i,]
        mat1 <- matrix(nrow = 1, ncol = 3*v)
        mat1[, 1:3] <- vec
        for (j in 4:(3*v)) {
          mat1[, j] <- (mat1[, j-3])%%v + 1
        }
        mat1
        SLR_1<-matrix(data=NA,nrow=v,ncol=3*v)
        SLR_1[1,]<-mat1
        for(t in 2:(v)){
          for(l in 1:(v-1)){
            SLR_1[t,(3*l-2):(3*l)]<-(SLR_1[t-1,((3*l+1):(3*l+3))])
            SLR_1[t,(3*v-2):(3*v)]<-(SLR_1[t-1,1:3])
          }
        }
        final=cbind(final,SLR_1)
      }
      SLR <-final
      SLR_4 <-matrix(NA,v,n*v)
      for (x in 1:(v*n)) {
        SLR_4[, x] <- paste(SLR[,3*x-2], SLR[,3*x-1],SLR[,3*x], sep = ",")
      }
      # noquote(SLR_4)
    }
    ## method4_Series3
    # n<-1 ##n should be any integer
    # v<-6*n+3
    n = (v-3)/6
    if(n - floor(n) == 0)
    {
      h<-6*n+3
      p<-(n+1)*(6*n+3)
      # k<-3
      nr<-3*(n+1)
      nc<-3
      r<-h*nr
      mat<-cbind(matrix(1,n+1,1),matrix(n+2,n+1,1),matrix((n+3):(2*n+3),n+1,1))
      final=NULL
      for(i in 1:nrow(mat)){
        vec <- mat[i,]
        mat1 <- matrix(nrow = 1, ncol = 3*v)
        mat1[, 1:3] <- vec
        for (j in 4:(3*v)) {
          mat1[, j] <- (mat1[, j-3])%%v + 1
        }
        mat1
        SLR_1<-matrix(data=NA,nrow=v,ncol=3*v)
        SLR_1[1,]<-mat1
        for(t in 2:(v)){
          for(l in 1:(v-1)){
            SLR_1[t,(3*l-2):(3*l)]<-(SLR_1[t-1,((3*l+1):(3*l+3))])
            SLR_1[t,(3*v-2):(3*v)]<-(SLR_1[t-1,1:3])
          }
        }
        final=cbind(final,SLR_1)
      }
      SLR <-final
      SLR_4 <-matrix(NA,v,(n+1)*v)
      for (x in 1:(v*(n+1))) {
        SLR_4[, x] <- paste(SLR[,3*x-2], SLR[,3*x-1],SLR[,3*x], sep = ",")
      }
      # noquote(SLR_43)
    }
    ##C_matrix
    vect<-c(t(SLR))
    N1<-design_to_N(matrix(vect,nrow = h,byrow = TRUE))
    N2<-nc*matrix(1,v,p)
    N3<-design_to_N(matrix(vect,nrow = h*p,byrow = TRUE))
    X1_prim_X1<-p*nc*diag(1,v)
    X1_prim_X2<-cbind(h*nr*c(rep.int(1,v)), N1, N2, N3)
    X2_prim_X2<-rbind(cbind(h*p*k,p*k*t(c(rep.int(1,h))),h*k*t(c(rep.int(1,p))),k*t(c(rep.int(1,h*p)))),cbind(p*k*c(rep.int(1,h)),p*k*diag(1,h),k*matrix(1,h,p),k*matrix(1,h,h*p)),cbind(h*k*c(rep.int(1,p)),k*matrix(1,p,h),h*k*diag(1,p),k*matrix(1,p,h*p)),cbind(k*c(rep.int(1,h*p)),k*matrix(1,h*p,h),k*matrix(1,h*p,p),k*diag(1,h*p)))
    c_matrix<-X1_prim_X1-X1_prim_X2%*%ginv(X2_prim_X2)%*%t(X1_prim_X2)
    eg=eigen(c_matrix)$values
    cef<-eg[eg>0.0001]/r    ##Canonical efficiency factors
    aef<-1/(mean(1/cef))  ##Average efficiency
    # aef
    output = list(design = noquote(SLR_4), Avg.Effi = aef)
  }
  return(output)
}