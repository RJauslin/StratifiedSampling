
sim <- function(n,f,Xaux,Xspread,D,pik,y){
  
  eps <- 1e-8
  np <- length(pik) # equal unequal
  ny <- length(y[[1]]) # number of variable of interest ( first level is for eq or unequal)
  
  
  W <- rep(list(0),np)
  for(k in 1:np){
    W[[k]] <- wpik(Xspread[[k]],pik[[k]])
    W[[k]] <- W[[k]] - diag(diag(W[[k]]))
  }
  
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
        
        s <- f[[i]](prob = pik[[k]],x = Xspread[[k]])
        s_01 <- rep(0,N)
        s_01[s] <- 1 
      }else if(names(f)[i] == "srswor"){
        
        # print("lpm1")
        s_01 <- sampling::srswor(sum(pik[[k]]),length(pik[[k]]))
        s <- which(s_01 > (1-eps))
        
        # s_01 <- f[[i]](pik[[k]])
        # s <- which(s_01 > (1- eps))
      }else if(names(f)[i] == "pwd") {
        s_01 <- rep(0,N)
        s <- f[[i]](D,sum(pik[[k]]))$s
        s_01[s] <- 1
      } else if (names(f)[i] == "stream"){
        s <- f[[i]](pik[[k]],Xaux[[k]],Xspread[[k]])
        s_01 <- rep(0,N)
        s_01[s] <- 1
      }else if(names(f)[i] == "lcube"){
        s <- f[[i]](pik[[k]],Xspread[[k]],Xaux[[k]])
        s_01 <- rep(0,N)
        s_01[s] <- 1 
      }
      
      
      tot <- colSums(as.matrix(Xaux[[k]]))
      est <- colSums(as.matrix(Xaux[[k]][s,])/pik[[k]][s])
      out_dev[[k]][[i]] <- 100*(est - tot)/tot
      
      
      # spread
      out_spread[[k]][[i]] <- rep(list(0),2)
      names(out_spread[[k]][[i]]) <- c("sb","IB")
      
      out_spread[[k]][[i]][[1]] <- sb(pik[[k]],Xspread[[k]],s)
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