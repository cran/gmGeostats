##### simulation algorithms ------------

### turning bands -----------

## move from geostats.R here


### LU Decomposition -----------

## TODO

### direct sampling ------




gsi.DS4CoDa <- function(n, f, t, n_realiz, nx_TI, ny_TI, nx_SimGrid, ny_SimGrid, TI_input, SimGrid_input,
                        V = "ilr", ivars_TI = 3:ncol(TI_input), 
                        W = t(gsi.produceV(V=V, D=length(ivars_TI), giveInv = TRUE)@W), 
                        SimGrid_mask = ncol(SimGrid_input), invertMask = TRUE
                        ){
  .Deprecated(new = "gsi.DS", 
              msg="gsi.DS4CoDa is deprecated; use gsi.DS or go via make.gm*-functions followed by DSpars and predict for gmSpatialModel")
  ### extract elements
  # mask
  if(length(SimGrid_mask)==1){
    mask = as.logical(SimGrid_input[, SimGrid_mask])
    if(is.numeric(SimGrid_mask)) SimGrid_input = SimGrid_input[,-SimGrid_mask]
    if(is.character(SimGrid_mask)) SimGrid_input = SimGrid_input[,setdiff(SimGrid_mask,colnames(SimGrid_input))]
  }else if(length(SimGrid_mask)==nrow(SimGrid_input)){
    mask = as.logical(SimGrid_mask)
  }else stop("gsi.DS4CoDa: SimGrid_mask not interpretable")
  if(invertMask) mask = !mask
  
  ### Prepare data matrices and eventually, transform
  if(V=="I"){
    ## case "no transformation"
    D=length(ivars_TI)
    # TI
    TI_ilr <- matrix(data = NA, nrow = nrow(TI_input), ncol = D)
    TI_ilr[which(complete.cases(TI_input)),] <- TI_input[which(complete.cases(TI_input)), ivars_TI]
    # conditioning data  + simgrid
    SimGrid_ilr <- matrix(data = NA, nrow = nrow(SimGrid_input), ncol = D)
    if (length(which(complete.cases(SimGrid_input)))>=1){
      SimGrid_ilr[which(complete.cases(SimGrid_input)),] <- 
        SimGrid_input[which(complete.cases(SimGrid_input)), ivars_TI]
      }
  }else{
    ## case "transform"
    # interpret V
    V = gsi.produceV(V=V, D=length(ivars_TI))
    # make space
    TI_ilr <- matrix(data = NA,nrow = nrow(TI_input),ncol = ncol(V))
    SimGrid_ilr <- matrix(data = NA, nrow = nrow(SimGrid_input), ncol = length(ivars_TI)-1)
    # transform TI
    TI_ilr[which(complete.cases(TI_input)),] <- ilr(TI_input[which(complete.cases(TI_input)),ivars_TI], V=V)
    # transform conditioning data
    if (length(which(complete.cases(SimGrid_input)))>=1){
      SimGrid_ilr[which(complete.cases(SimGrid_input)),] <- 
        ilr(SimGrid_input[which(complete.cases(SimGrid_input)), ivars_TI])
      }
  }
  
  # Array to store realizations
  SimGrid_ilr <- replicate(n_realiz, SimGrid_ilr)
  
  # If conditioning data is not provided, proceed with nonconditional simulation   
  tk0 = complete.cases(SimGrid_input)
  if (sum(tk0)<n){
    for (i in 1:n_realiz){
      nmissing = n-sum(tk0)
      stk = sample(x=which(mask &!tk0), size=nmissing)
      SimGrid_ilr[stk,,i] <- TI_ilr[sample(x=which(complete.cases(TI_ilr)),size = nmissing),]
    }
  }
  
  # Compositional range
  tkTI = complete.cases(TI_ilr)
  CRange <- max(dist(TI_ilr[tkTI,]))
  
  # Change TI to an array
  TI_ilr_array <- array(as.vector(TI_ilr),dim = c(nx_TI,ny_TI, ncol(TI_ilr)))
  
  # List to store realization
  SimGrid_ilr_list <- list()
  for (i in 1:n_realiz){
    SimGrid_ilr_list[[i]] <- array(as.vector(SimGrid_ilr[,,i]), dim = c(nx_SimGrid,ny_SimGrid,ncol(SimGrid_ilr)))
  }
  
  # matrix of the informed nodes in the training image
  mInformedTI  <- which(!is.na(TI_ilr_array[,,1]), arr.ind = TRUE)
  
  # array of nodes to be simulated
  maskArray <- array(mask, dim = c(nx_SimGrid,ny_SimGrid,1))
  
  #pb = list()
  #myfun = function(ii){
  for(ii in 1:n_realiz){
    cat(paste("\n Realization number #",ii, "\n"))
    
    # Defining a fully random path for simulation
    list_sim <- which(maskArray[,,1] & is.na(SimGrid_ilr_list[[ii]][,,1]), arr.ind = TRUE)
    path_sim <- list_sim[sample(nrow(list_sim)),]
    
    # initialize progress bar
    pb <- utils::txtProgressBar(min = 0, max = nrow(path_sim)*f*nrow(mInformedTI), style = 3)
    status <- 0
    # Looping simulation nodes
    for (simnod in 1:nrow(path_sim)){
      path_this_sim = path_sim[simnod,]
      
      # Finding the n closest compositions (hard or simulated) to build the data event
      tki = !is.na(SimGrid_ilr_list[[ii]][,,1])
      dataevesim_discode <- FNN::get.knnx(
          data=which(tki,arr.ind = TRUE), t(as.matrix(path_this_sim)), 
          k=n, 
          algorithm=c("kd_tree")
        )
      dataevesim_loc <- which(tki, arr.ind = TRUE)[dataevesim_discode$nn.index,]
      dataevesim <- mapply(function(i, j) SimGrid_ilr_list[[ii]][i, j, 1:ncol(TI_ilr)], dataevesim_loc[,1], dataevesim_loc[,2])
      dataevesim_vec <- dataevesim_loc - matrix(rep( t(as.matrix(path_this_sim)),each=n),nrow=n)
      
      # Scanning TI for a close pattern
      path_TI <- mInformedTI[sample(nrow(mInformedTI)),]
      
      # Initial best distance is set to inf. Update with every best distance encountered
      mindist <- Inf
      
      # Number of tries in the TI
      nb_of_tries <- ceiling(nrow(path_TI)*f)
      # Store best pattern encountered so far
      BestPoint <- matrix(data = NA,nrow = 1,ncol = 2)
      
      for (tinod in 1:nb_of_tries){
        # update progress bar
        status = status + 1
        utils::setTxtProgressBar(pb, status)
        
        # Building training pattern and measuring distance
        dataeveti_loc <- dataevesim_vec + matrix(rep( t(as.matrix(path_TI[tinod,])),each=n),nrow=n) 
        outwin <- dataeveti_loc[,1] <= nx_TI & dataeveti_loc[,2] <= ny_TI & dataeveti_loc[,1] > 0 & dataeveti_loc[,2] > 0 
        if(sum(outwin)==0){next}
        dataeveti <- mapply(function(i, j) TI_ilr_array[i, j, 1:ncol(TI_ilr)], dataeveti_loc[outwin,1], dataeveti_loc[outwin,2])
        if(sum(is.na(dataeveti[1,]))>=ncol(dataeveti)){next}
        mydist <- mean(sqrt(colSums((dataevesim[,outwin] - dataeveti)^2))/CRange,na.rm = TRUE)
        
        # Checking for the minimum distance found so far
        if (mydist < mindist){
          mindist <- mydist
          BestPoint <- t(as.matrix(path_TI[tinod,]))
        }
        # break the loop if the distance is less than t
        if (mindist <= t){break}
      }
      # update status bar
      status = simnod* f*nrow(mInformedTI)
      utils::setTxtProgressBar(pb, status)
      
  
      # pasting the whole composition
      SimGrid_ilr_list[[ii]][path_this_sim[1],path_this_sim[2],] <- TI_ilr_array[BestPoint[,1],BestPoint[,2],]
      #return(SimGrid_ilr_list[[ii]])
    }
  }
  
  #SimGrid_ilr_list = foreach(ii=1:n_realiz,.combine = list) %dopar% myfun(ii) 
  
    
  # Empty array to store the backtransfomed realizations
  SimGrid <- array(data = NA, dim = c(nrow(SimGrid_ilr),ncol(SimGrid_ilr)+1, n_realiz))
  
  # Backtransform to compositional space
  if(is.matrix(W)){
    for (i in 1:n_realiz){
      SimGrid_ilr[,,i] <- matrix(as.vector(SimGrid_ilr_list[[i]]),nrow = nrow(SimGrid_ilr),ncol = ncol(SimGrid_ilr))
      SimGrid[mask,,i] <- ilrInv(SimGrid_ilr[mask,,i], V=W)
      varnames_out = rownames(W)
    } 
  }else{
    for (i in 1:n_realiz){
      SimGrid[mask,,i] <- matrix(as.vector(SimGrid_ilr_list[[i]]),nrow = nrow(SimGrid_ilr),ncol = ncol(SimGrid_ilr))
      varnames_out = colnames(TI_input)[ivars_TI]
    }
  }
  

  # addition by Raimon 20200402
  ddd = dim(SimGrid)[2]
  if(length(varnames_out)!=ddd)   
    varnames_out[is.null(varnames_out)] = paste("v", 1:ddd, sep="")[is.null(varnames_out)]
  dimnames(SimGrid) = list(loc=1:nrow(SimGrid_ilr),
                           var=varnames_out,
                           sim=paste("sim", 1:n_realiz, sep="")
  )
  SimGrid=DataFrameStack(SimGrid, stackDim="sim")
  
  return(SimGrid)  
}





#' Workhorse function for direct sampling
#' 
#' This function implements in R the direct sampling algorithm
#'
#' @param n size of the conditioning data event (integer)
#' @param f fraction of the training image to scan (numeric between 0 and 1)
#' @param t maximal acceptable discrepance between conditioning data event and TI event (numeric between 0 and 1)
#' @param n_realiz number of simulations desired
#' @param dim_TI dimensions of the grid of the training image (ie. either \eqn{(n_x, n_y)} 
#' for dimension \eqn{k=2} or \eqn{(n_x, n_y, n_z)} for dimension \eqn{k=3})
#' @param dim_SimGrid  dimensions of the simulation grid (ie. either \eqn{(m_x, m_y)} or 
#' \eqn{(m_x, m_y, m_z)})
#' @param TI_input training image, as a matrix of \eqn{(n_x\cdot n_y\cdot n_z, k+D)} 
#' elements; WITH NAMED COLUMNS and including spatial coordinates
#' @param SimGrid_input simulation grid with conditioning data, as a matrix of 
#' \eqn{(m_x\cdot m_y\cdot m_z, k+D)} elements; with same columns as `TI_input`
#' @param ivars_TI which colnames of `TI_input` and `SimGrid_input` identify variables to consider in the data event
#' @param SimGrid_mask either a logical vector of length \eqn{m_x\cdot m_y\cdot m_z}, or else a column name of `SimGrid_input` 
#' giving a logical column
#' @param invertMask logical, does `SimGrid_mask` identify with TRUE the data OUTSIDE the simulation area?
#'
#' @return A [sp::SpatialPixelsDataFrame()] or  [sp::SpatialGridDataFrame()], depending on whether the whole
#' grid is simulated. The '@data' slot of these objects contains a [DataFrameStack()] with the stacking dimension
#' running through the realisations. It is safer to use this functionality through the interface
#' [make.gmCompositionalMPSSpatialModel()], then request a direct simulation with [DSpars()] and
#' finally run it with [predict_gmSpatialModel].
#' @export
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom stats complete.cases lm 
#' @author Hassan Talebi (copyright holder), Raimon Tolosana-Delgado
#' @examples
#' ## training image:
#' x = 1:10
#' y = 1:7
#' xy_TI = expand.grid(x=x, y=y)
#' TI_input = cbind(xy_TI, t(apply(xy_TI, 1, function(x) c(sum(x), abs(x[2]-x[1]))+rnorm(2, sd=0.01))))
#' colnames(TI_input) = c("x", "y", "V1", "V2")
#' o1 = image_cokriged(TI_input, ivar="V1")
#' o2 = image_cokriged(TI_input, ivar="V2")
#' ## simulation grid:
#' SimGrid = TI_input
#' SimGrid$mask = with(SimGrid, x==1 | x==10 | y==1 | y==7)
#' tk = SimGrid$mask
#' tk[sample(70, 50)] = TRUE 
#' SimGrid[tk,3:4]=NA
#' image_cokriged(SimGrid, ivar="V1", breaks=o1$breaks, col=o1$col)
#' image_cokriged(SimGrid, ivar="V2", breaks=o2$breaks, col=o2$col)
#' image_cokriged(SimGrid, ivar="mask", breaks=c(-0.0001, 0.5, 1.001))
#' ## res = gsi.DS(n=5, f=0.75, t=0.05, n_realiz=2, dim_TI=c(10,7),  dim_SimGrid=c(10,7), 
#' ##        TI_input=as.matrix(TI_input), SimGrid_input=as.matrix(SimGrid), 
#' ##        ivars_TI = c("V1", "V2"), SimGrid_mask="mask", invertMask=TRUE)
#' ## image_cokriged(cbind(xy_TI, getStackElement(res,1)), ivar="V1", breaks=o1$breaks, col=o1$col)
#' ## image_cokriged(cbind(xy_TI, getStackElement(res,2)), ivar="V1", breaks=o1$breaks, col=o1$col)
#' ## image_cokriged(cbind(xy_TI, getStackElement(res,1)), ivar="V2", breaks=o2$breaks, col=o2$col)
#' ## image_cokriged(cbind(xy_TI, getStackElement(res,2)), ivar="V2", breaks=o2$breaks, col=o2$col)
gsi.DS <- function(n, f, t, n_realiz, 
                   dim_TI, dim_SimGrid,  
                   TI_input, SimGrid_input,
                   ivars_TI = 3:ncol(TI_input), 
                   SimGrid_mask = ncol(SimGrid_input), invertMask = TRUE
){
  if(!requireNamespace("FNN", quietly = TRUE)) stop("direct sampling requires package 'KNN' installed")
  ## constants:
  nx_TI=dim_TI[1]
  ny_TI=dim_TI[2]
  if(length(dim_TI)>2){nz_TI=dim_TI[3]}else{nz_TI=1}
  nx_SimGrid=dim_SimGrid[1]
  ny_SimGrid=dim_SimGrid[2]
  if(length(dim_SimGrid)>2){nz_SimGrid=dim_SimGrid[3]}else{nz_SimGrid=1}

  ### extract elements out of TI and simgrid
  # mask
  if(length(SimGrid_mask)==1){
    mask = as.logical(SimGrid_input[, SimGrid_mask])
    if(is.numeric(SimGrid_mask)) SimGrid_input = SimGrid_input[,-SimGrid_mask]
    if(is.character(SimGrid_mask)) SimGrid_input = SimGrid_input[,setdiff(colnames(SimGrid_input), SimGrid_mask)]
  }else if(length(SimGrid_mask)==nrow(SimGrid_input)){
    mask = as.logical(SimGrid_mask)
  }else stop("gsi.DS4CoDa: SimGrid_mask not interpretable")
  if(invertMask) mask = !mask
  # full grid
  if(is.numeric(ivars_TI)) fullgrid = SimGrid_input[,-ivars_TI]
  if(is.character(ivars_TI)) fullgrid = SimGrid_input[,setdiff(colnames(SimGrid_input), ivars_TI)]
  # nr of variables
  D=length(ivars_TI)
  # TI
  TI <- matrix(data = NA, nrow = nrow(TI_input), ncol = D)
  TI[which(complete.cases(TI_input)),] <- TI_input[which(complete.cases(TI_input)), ivars_TI]
  # conditioning data  + SimGrid
  SimGrid <- matrix(data = NA, nrow = nrow(SimGrid_input), ncol = D)
  tk0 = complete.cases(SimGrid_input)
  if (sum(tk0)>=1)  SimGrid[tk0,] <- SimGrid_input[tk0, ivars_TI]
  
  # Array to store realizations
  SimGrid <- replicate(n_realiz, SimGrid)
  
  # If not sufficient conditioning data are provided, complete with some unconditional simulations   
  if (sum(tk0)<n){
    for (i in 1:n_realiz){
      nmissing = n-sum(tk0)
      stk = sample(x=which(mask &!tk0), size=nmissing)
      SimGrid[stk,,i] <- TI[sample(x=which(complete.cases(TI)), size = nmissing),]
    }
  }
  
  # Data range
  tkTI = complete.cases(TI)
  CRange <- max(dist(TI[tkTI,]))
  
  # Change TI to an array
  TI_array <- array(as.vector(TI),dim = c(nx_TI,ny_TI, nz_TI, ncol(TI)))
  
  # List to store realization
  SimGrid_list <- list()
  for (i in 1:n_realiz){
    SimGrid_list[[i]] <- array(as.vector(SimGrid[,,i]), dim = c(nx_SimGrid,ny_SimGrid,nz_SimGrid,D) )
  }
  
  # matrix of the informed nodes in the training image
  TIinformed_array  <- which(!is.na(TI_array[,,,1]), arr.ind = TRUE)
  
  # array of nodes to be simulated
  mask_array <- array(mask, dim = c(nx_SimGrid,ny_SimGrid,nz_SimGrid))
  
  #pb = list()
  #myfun = function(ii){
  for(ii in 1:n_realiz){
    cat(paste("\n Realization number #",ii, "\n"))
    
    # Defining a fully random path for simulation
    list_sim <- which(mask_array[,,,drop=T] & is.na(SimGrid_list[[ii]][,,,1]), arr.ind = TRUE)
    path_sim <- list_sim[sample(nrow(list_sim)),]
    
    # initialize progress bar
    pb <- utils::txtProgressBar(min = 0, max = nrow(path_sim)*f*nrow(TIinformed_array), style = 3)
    status <- 0
    # Looping simulation nodes
    for (simnod in 1:nrow(path_sim)){
      path_this_sim = path_sim[simnod,]
      
      # Finding the n closest compositions (hard or simulated) to build the data event
      tki = !is.na(SimGrid_list[[ii]][,,,1])
      dataevesim_discode <- FNN::get.knnx(
        data=which(tki,arr.ind = TRUE), t(as.matrix(path_this_sim)),  # why t(as.matrix(...))?
        k=n, 
        algorithm=c("kd_tree")
      )
      dataevesim_loc <- which(tki, arr.ind = TRUE)[c(dataevesim_discode$nn.index),]
      if(ncol(dataevesim_loc)==3){
        G = 3
        dataevesim <- mapply(function(i, j, k) SimGrid_list[[ii]][i, j, k, 1:D], dataevesim_loc[,1], dataevesim_loc[,2], dataevesim_loc[,3])        
      }else{
        G = 2
        dataevesim <- mapply(function(i, j) SimGrid_list[[ii]][i, j, 1, 1:D], dataevesim_loc[,1], dataevesim_loc[,2])
      }
      dataevesim_vec <- dataevesim_loc - matrix(rep( t(as.matrix(path_this_sim)),each=n),nrow=n) # compute lag constellation 
      
      # Scanning TI for a close pattern
      path_TI <- TIinformed_array[sample(nrow(TIinformed_array)),]
      
      # Initial best distance is set to inf. Update with every best distance encountered
      mindist <- Inf
      
      # Number of tries in the TI
      nb_of_tries <- ceiling(nrow(path_TI)*f)
      # Store best pattern encountered so far
      BestPoint <- matrix(data = NA, nrow = 1, ncol = G)
      
      for (tinod in 1:nb_of_tries){
        # update progress bar
        status = status + 1
        utils::setTxtProgressBar(pb, status)
        
        # Building training pattern and measuring distance
        dataeveti_loc <- dataevesim_vec + matrix(rep( t(as.matrix(path_TI[tinod,])),each=n),nrow=n)  # place the lag constellation on the training image
        outwin <- dataeveti_loc[,1] <= nx_TI & dataeveti_loc[,2] <= ny_TI &  dataeveti_loc[,1] > 0 & dataeveti_loc[,2] > 0  
        if(G==3) outwin <- outwin & dataeveti_loc[,3] <= nz_TI  & dataeveti_loc[,3] > 0 
        if(sum(outwin)==0){next}
        if(G==3){
          dataeveti <- mapply(function(i, j, k) TI_array[i, j, k, 1:ncol(TI)], dataeveti_loc[outwin,1], dataeveti_loc[outwin,2], dataeveti_loc[outwin,3])
        }else{
          dataeveti <- mapply(function(i, j) TI_array[i, j, 1, 1:ncol(TI)], dataeveti_loc[outwin,1], dataeveti_loc[outwin,2])
        }
        if(sum(is.na(dataeveti[1,]))>=ncol(dataeveti)){next}
        mydist <- mean(sqrt(colSums((dataevesim[,outwin] - dataeveti)^2))/CRange,na.rm = TRUE)
        
        # Checking for the minimum distance found so far
        if (mydist < mindist){
          mindist <- mydist
          BestPoint <- t(as.matrix(path_TI[tinod,]))
        }
        # break the loop if the distance is less than t
        if (mindist <= t){break}
      }
      # update status bar
      status = simnod* f*nrow(TIinformed_array)
      utils::setTxtProgressBar(pb, status)
      
      
      # pasting the whole composition
      if(G==3){
        SimGrid_list[[ii]][path_this_sim[1],path_this_sim[2],path_this_sim[3],] <- TI_array[BestPoint[,1],BestPoint[,2],BestPoint[,3],]        
      }else{
        SimGrid_list[[ii]][path_this_sim[1],path_this_sim[2],1,] <- TI_array[BestPoint[,1],BestPoint[,2],1,]
      }

      #return(SimGrid_ilr_list[[ii]])
    }
  }
  
  # set as DataFrameStack
  if(is.numeric(ivars_TI)){
    varnames_out = tryCatch(colnames(TI_input)[ivars_TI])
  }else if(is.character(ivars_TI)){
     varnames_out = ivars_TI
  }
  if(length(varnames_out)!=D | class(varnames_out)=="try-error") varnames_out = paste("v", 1:D, sep="")
  dm = list(loc=1:length(mask), var=varnames_out, sim=paste("sim", 1:n_realiz, sep="") )
  
  SimGrid_list = lapply(SimGrid_list, function(x){
    dim(x) = c(length(x)/D,D)
    rownames(x) = 1:length(mask)
    colnames(x) = varnames_out
    x
  })
  
  SimGrid=DataFrameStack(SimGrid_list, stackDimName="sim", Dimnames=dm)
  

  return(SimGrid)  
}

