
sim <- function(n,f,Xaux,Xspread,D,pik,y){
  
  
  # only used for  HIP
  findIndex <- function(x1,X){
    final <- c()
    for(i in 1:nrow(x1)){
      id <- seq(1,nrow(X),1)
      for(j in 1:ncol(x1)){
        out = which(x1[i,j] == X[id,j])
        if(length(out) > 1){
          id <- out
        }else if(length(out) == 1 & j == 1){
          final <- c(final,out)
          break;
        }else if(length(out) == 1){
          final <- c(final,id[out])
          break;
        }
      }
    }
    return(final)
  }
  
  # initialization 
  eps <- 1e-8
  np <- length(pik) # equal unequal
  ny <- length(y[[1]]) # number of variable of interest ( first level is for eq or unequal)
  
  X_Sp  <- sp::SpatialPoints(Xspread)
  
  W <- rep(list(0),np)
  for(k in 1:np){
    W[[k]] <- wpik(Xspread,pik[[k]])
    W[[k]] <- W[[k]] - diag(diag(W[[k]]))
  }
  print(W[[k]])
  # initialize containers 
  out_v <- rep(list(rep(list(0),length(f))),np)
  names(out_v) <- names(pik)
  
  out_spread <- rep(list(rep(list(0),length(f))),np)
  names(out_spread) <- names(pik)
  
  out_dev <- rep(list(rep(list(0),length(f))),np)
  names(out_dev) <- names(pik)
  
  for(k in 1:np){# loop for equal unequl pik
   
     N <- length(pik[[k]])
    names(out_spread[[k]]) <- names(out_v[[k]]) <- names(out_dev[[k]]) <- names(f)
    
    for(i in 1:length(f)){  # loop on the design
      
      # sample selection
      if(names(f)[i] == "lpm1"){
        
        s <- f[[i]](prob = pik[[k]],x = Xspread)
        s_01 <- rep(0,N)
        s_01[s] <- 1 
      }else if(names(f)[i] == "srswor"){
        
        # print("lpm1")
        s_01 <- sampling::srswor(sum(pik[[k]]),length(pik[[k]]))
        s <- which(s_01 > (1-eps))
        
        # s_01 <- f[[i]](pik[[k]])
        # s <- which(s_01 > (1- eps))
      }else if(names(f)[i] == "pwd") {
        con <- rep(0,nrow(D))
        stand_D <- stprod(mat = D, con = con)$mat
        s_01 <- rep(0,N)
        s <- f[[i]](stand_D,sum(pik[[k]]))$s
        s_01[s] <- 1
      } else if (names(f)[i] == "stream"){
        s <- f[[i]](pik[[k]],Xaux[[k]],Xspread)
        s_01 <- rep(0,N)
        s_01[s] <- 1
      }else if(names(f)[i] == "lcube"){
        s <- f[[i]](pik[[k]],Xspread,Xaux[[k]])
        s_01 <- rep(0,N)
        s_01[s] <- 1 
      }else if(names(f)[i] == "wave"){
       s_01 <- f[[i]](Xspread,pik[[k]])
       s <- which(s_01 == 1)
      }else if(names(f)[i] == "hip"){
        tmp1 = SDraw::hip.point(x = X_Sp,n = sum(pik[[k]]),plot.lattice = FALSE)
        tmp1 <- sp::coordinates(tmp1)
        s_01 <- rep(0,N)
        s_01[findIndex(tmp1,Xspread)] <- 1
        s <- which(s_01 == 1)
      }
      
      tot <- colSums(as.matrix(Xaux[[k]]))
      est <- colSums(as.matrix(Xaux[[k]][s,])/pik[[k]][s])
      out_dev[[k]][[i]] <- 100*(est - tot)/tot
      
      
      # spread
      out_spread[[k]][[i]] <- rep(list(0),2)
      names(out_spread[[k]][[i]]) <- c("sb","IB")
      
      out_spread[[k]][[i]][[1]] <- sb(pik[[k]],Xspread,s)
      out_spread[[k]][[i]][[2]] <- IB(W[[k]],s_01)
      
      # variance and variable of interest
      out_v[[k]][[i]] <- rep(list(0),ny)
      names(out_v[[k]][[i]]) <- names(y[[k]])
      
      for(j in 1:length(y[[k]])){ # loop on the variable of interests
        out_v[[k]][[i]][[j]] <- (sum(y[[k]][[j]][s]/pik[[k]][s]) - sum(y[[k]][[j]]))^2
      }
    }
  }
  return(list(spread = out_spread,
              var = out_v,
              dev = out_dev))
}




#' Title
#'
#' @param n 
#' @param f list of functions 
#' @param Xaux matrx
#' @param Xspread matrix
#' @param D matrix of distance
#' @param pik vector of probabilities
#' @param y list of vartiable of interest
#'
sim_design <- function(n,f,Xaux,Xspread,D,pik,y){
  
  
  # only used for  HIP
  findIndex <- function(x1,X){
    final <- c()
    for(i in 1:nrow(x1)){
      id <- seq(1,nrow(X),1)
      for(j in 1:ncol(x1)){
        out = which(x1[i,j] == X[id,j])
        if(length(out) > 1){
          id <- out
        }else if(length(out) == 1 & j == 1){
          final <- c(final,out)
          break;
        }else if(length(out) == 1){
          final <- c(final,id[out])
          break;
        }
      }
    }
    return(final)
  }
  
  # initialization 
  eps <- 1e-8
  ny <- length(y) # number of variable of interest ( first level is for eq or unequal)
  
  X_Sp  <- sp::SpatialPoints(Xspread)
  
  W <- wpik(Xspread,pik)
  W <- W - diag(diag(W))
  
  # initialize containers 
  out_v <- rep(list(0),length(f))
  out_spread <- rep(list(0),length(f))
  out_dev <- rep(list(0),length(f))
  out_dev <- matrix(rep(0,length(f)*ncol(Xaux)),ncol = ncol(Xaux))
  rownames(out_dev) <- names(f)

  out_vEst <- rep(list(0),length(f))
  out_HT <- rep(list(0),length(f))

  N <- length(pik)
  names(out_HT) <- names(out_vEst) <- names(out_spread) <- names(out_v) <- names(out_dev) <- names(f)
  
  
  # Approximated Variance does not depends on the sampling design 
  out_vApp <- rep(list(0),ny)
  names(out_vApp) <- names(y)
  for(j in 1:ny){
    out_vApp[[j]] <- StratifiedSampling::vApp(Xaux,pik,y[[j]])
  }
  
  # loop on the different design
  for(i in 1:length(f)){  # loop on the design
    
    # sample selection
    if(names(f)[i] == "lpm1"){
      
      s <- BalancedSampling::lpm1(prob = pik,x = Xspread)
      s_01 <- rep(0,N)
      s_01[s] <- 1 
    }else if(names(f)[i] == "srswor"){
    
      s_01 <- sampling::srswor(sum(pik),length(pik))
      s <- which(s_01 > (1-eps))
      
      # s_01 <- f[[i]](pik[[k]])
      # s <- which(s_01 > (1- eps))
    }else if(names(f)[i] == "cps"){
      
      s_01 <- cps(pik)
      s <- which(s_01 > (1-eps))
      
    }else if(names(f)[i] == "pwd") {
      con <- rep(0,nrow(D))
      stand_D <- stprod(mat = D, con = con)$mat
      s_01 <- rep(0,N)
      s <- f[[i]](stand_D,sum(pik))$s
      s_01[s] <- 1
    } else if (names(f)[i] == "stream"){
      s <- balseq(pik,Xaux,Xspread)
      s_01 <- rep(0,N)
      s_01[s] <- 1
    }else if(names(f)[i] == "lcube"){
      s <- f[[i]](pik,Xspread,Xaux)
      s_01 <- rep(0,N)
      s_01[s] <- 1 
    }else if(names(f)[i] == "wave"){
      s_01 <- f[[i]](Xspread,pik)
      s <- which(s_01 == 1)
    }else if(names(f)[i] == "hip"){
      tmp1 = try(SDraw::hip.point(x = X_Sp,n = sum(pik),plot.lattice = FALSE),
                         silent = TRUE)
      while(class(tmp1)=="try-error"){
        print("err in hip")
        tmp1 = try(SDraw::hip.point(x = X_Sp,n = sum(pik),plot.lattice = FALSE),
                        silent = TRUE)
      }
      tmp1 <- sp::coordinates(tmp1)
      s_01 <- rep(0,N)
      s_01[findIndex(tmp1,Xspread)] <- 1
      s <- which(s_01 == 1)
    }else if(names(f)[i] == "cube"){
      s <- BalancedSampling::cube(pik,Xaux)
      s_01 <- rep(0,N)
      s_01[s] <- 1
    }
    
    # dev
    tot <- colSums(as.matrix(Xaux))
    est <- colSums(as.matrix(Xaux[s,])/pik[s])
    out_dev[i,] <- 100*(est - tot)/tot
    
    
    # spread
    out_spread[[i]] <- rep(list(0),2)
    names(out_spread[[i]]) <- c("sb","IB")
    out_spread[[i]][[1]] <- sb(pik,Xspread,s)
    out_spread[[i]][[2]] <- IB(W,s_01)
    
    # # variance and variable of interest, vEst, vApp
    # out_HT[[i]] <- out_vEst[[i]] <- out_vDbs[[i]] <- out_v[[i]] <- rep(list(0),ny)
    # names(out_HT[[i]]) <- names(out_vEst[[i]]) <- names(out_vDbs[[i]]) <- names(out_v[[i]]) <- names(y)
    # 
    # for(j in 1:length(y)){ # loop on the variable of interests
    #   out_HT[[i]] <- sum(y[[j]][s]/pik[s])
    #   out_vEst[[i]][[j]] <- StratifiedSampling::vEst(Xaux,pik,y[[j]],s_01)
    #   out_vDbs[[i]][[j]] <- StratifiedSampling::varDBS(Xaux[s,],Xspread[s,],pik[s],y[[j]][s])
    #   out_v[[i]][[j]] <- (sum(y[[j]][s]/pik[s]) - sum(y[[j]]))^2
    # }
    
    # variance simulated
    out_HT[[i]] <- out_v[[i]] <- rep(list(0),ny)
    names(out_HT[[i]]) <- names(out_v[[i]]) <- names(y)
    
    for(j in 1:length(y)){ # loop on the variable of interests
      out_HT[[i]] <- sum(y[[j]][s]/pik[s])
      out_v[[i]][[j]] <- (sum(y[[j]][s]/pik[s]) - sum(y[[j]]))^2
    }
    
    # variance estimator
    out_vEst[[i]] <- rep(list(0),ny)
    names(out_vEst[[i]]) <-  names(y)
    
    if(names(f)[i] == "lpm1" | names(f)[i] == "wave" | names(f)[i] == "hip" | names(f)[i] == "pwd" ){
      for(j in 1:length(y)){ 
        out_vEst[[i]][[j]] <- BalancedSampling::vsb(pik[s],y[[j]][s],Xspread[s,])
      }
    }else if (names(f)[i] == "cps" | names(f)[i] == "srswor"){
      for(j in 1:length(y)){ 
        out_vEst[[i]][[j]] <- WaveSampling::varHAJ(y[[j]][s],pik[s],s)
      }
    }else if(names(f)[i] == "stream" | names(f)[i] == "lcube"){
      for(j in 1:length(y)){ 
        out_vEst[[i]][[j]] <- StratifiedSampling::vDBS(Xaux[s,],Xspread[s,],pik[s],y[[j]][s])
      }
    }else if(names(f)[i] == "cube"){
      for(j in 1:length(y)){ 
      # out_vEst[[i]][[j]] <- StratifiedSampling::vEst(Xaux,pik,y[[j]],s_01)
      out_vEst[[i]][[j]] <- StratifiedSampling::vEst(Xaux[s,],pik[s],y[[j]][s])
      }
    }
    
    
    
  }

  return(list(spread = out_spread,
              HT = out_HT,
              v = out_v,
              vEst = out_vEst,
              vApp = out_vApp,
              dev = as.data.frame(out_dev)))
}
