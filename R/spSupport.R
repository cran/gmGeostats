par_remove_readonly = function(x, args=c("cin", "cra", "csi", "cxy", "din", "fig", "page")){
  x[args]<- NULL
  return(x)
}


oneimage = function(x,y,Z,isim,ivar,asp=1,...){
  aux = Z[,ivar,isim]
  dim(aux) = c(length(x),length(y))
  image(x,y,aux,asp=asp,...)
}


#' Write a regionalized data set in GSLIB format
#' 
#' Write a regionalized data set in plain text GSLIB format

#' @param x regionalized data set
#' @param file filename
#' @param header the first line of text for the file, defaults to filename  
#'
#' @return The status of closing the file, see \code{\link{close}} 
#' for details, although this is seldom problematic. This function is basically called
#' for its side-effect of writing a data set in the simplified Geo-EAS format that is 
#' used in GSLIB.
#' 
#' @seealso \url{http://www.gslib.com/gslib_help/format.html}
#' 
#' @export
#' @importFrom utils write.table 
#'
#' @examples
#' data("jura", package="gstat")
#' \dontrun{write.GSLib(jura.pred, file="jurapred.txt")}
write.GSLib <- function(x,file, header=basename(file)) {
  if(is(x, "Spatial")){
    if("data" %in% slotNames(x)) y = x@data else y=NULL
    x = cbind(data.frame(sp::coordinates(x)),y)
  } 
  x = as.data.frame(lapply(x,as.numeric))
  x[is.na(x)] = -1e99
  f<-file(file,"w")
  cat(file=f,header,"\n")
  cat(file=f,ncol(x),"\n")
  write(file=f,colnames(x),sep="\n",append=TRUE)
  write.table(file=f,x,append=TRUE,col.names=FALSE,row.names=FALSE)
  close(f)
}


ISATISrotation = function(Az, Ay, Ax){
	myfun = function(a, i, j){
		m = diag(3)
		a = pi*a/180
		m[c(i,j),c(i,j)] = matrix(c(cos(a), sin(a), -sin(a), cos(a)), ncol=2, nrow=2)
		return(m)
	}
	M = myfun(Ax, 2,3) %*% myfun(Ay, 3,1) %*% myfun(Az, 1,2)
	#M = myfun(Az, 1,2) %*%  myfun(Ay, 3,1) %*% myfun(Ax, 2,3) 
	return(M)
 }






#################################################
## function to generate a polarplot representing
## the number of pairs of data in each lag class
pointpairs2polargrid = function(loc, # 
                                maxdist=max(dist(loc))/2, nbins=10,  # h maximal value and nr of classes
                                dists = seq(0,maxdist,length.out=nbins+1), # h classes to consider
                                azimuths=(0:11)*30,   # azimuths to consider
                                plotit=TRUE){
  # calculate the distances between all points, in radius and in angle
  dr = as.matrix(dist(loc))
  da = outer(1:nrow(loc), 1:nrow(loc), function(i,j) atan2(x=loc[j,2]-loc[i,2], y=loc[j,1]-loc[i,1]) )
  # count the number of radial distances in each class
  distmids = dists[-1] - diff(dists)/2
  comparisonR = sapply(c(dr), function(dd){
    aux = abs(dd-distmids)
    which(aux==min(aux))[1]
  })
  # count the number of angular distances in each class
  angdist = function(a,b) gmApply(cbind(abs(a-b), abs(a-b+2*pi), abs(a-b-2*pi)), 1, min) # construct an angular distance function
  azimuths = azimuths * pi/180
  comparisonA = outer(c(da), azimuths, "-" ) %% (2*pi)  # angular distances are modulo 2pi
  comparisonA = gmApply(comparisonA, 1, function(x) which(x==min(x)) )
  # result
  tt = table(comparisonR, comparisonA)
  colnames(tt) = azimuths
  rownames(tt) = distmids
  if(plotit){
    dfaz = diff(azimuths)[1]
    aux = c(azimuths[1]-0.5*dfaz, azimuths+0.5*dfaz)
    aux = pi/2-aux
    image.polargrid(r = dists, phi = aux, tt, breaks = unique(quantile(tt[tt!=0], probs=seq(0,1,0.1))))
  }
  return(tt)
}
#################################################

#################################################
## ensure that azimuths are numeric
gsi.azimuth2angle = function(x) {
  if(is.numeric(x)) return(x)
  # if "N"
  x = sub("N","",x)
  return(as.numeric(x))
}
#################################################



#################################################
## internal alternative to paint one single pair

imageOneAnisVario = function(avg, v1=NULL, v2=NULL){
  # construct the polar grid
  r = attr(avg, "dists")
  rlim = range(r)
  phi = as.double(colnames(avg))
  dphi = mean(diff(phi))
  phi = phi - dphi/2
  phi = c(phi,phi[length(phi)]+dphi)*pi/180
  phi2 = pi/2-phi
  # extract the dimension of the composition
  D = 1
  # extract common z levels
  aux = c(unlist(sapply(1:ncol(avg), function(i) avg["vg",i][[1]][,v1,v2])))
  bks = quantile(aux[aux!=0], probs=seq(0,1,0.1), na.rm=TRUE)
  # set the matrix of figures
  plot(c(-1,1)*rlim[2], c(-1,1)*rlim[2],type="n", asp=1, xlab="", ylab="", bty="n", xaxt="n", yaxt="n", ann=FALSE, xaxs="i",yaxs="i")
  # extract the directional variogram for the logratio k vs j
  z = sapply(1:ncol(avg), function(i){
    avg["vg",i][[1]][,v1,v2]
  })
  z[is.nan(z)]=NA
  dim(z) = c(length(r)-1, length(phi)-1)
  # plot it
  image.polargrid(r,phi2,z,add=TRUE)
}
#################################################



#################################################
## Utility function to 
#### generate a colored polar sector diagram ---------
image.polargrid = function(
  r = seq(0, 1, length.out = nrow(z)), # radii
  phi = seq(0, 2*pi, length.out = ncol(z)),  #  angles 
  z,  # elevation
  zlim = range(z[is.finite(z)], na.rm=TRUE),
  rlim = range(r, na.rm=TRUE), philim = range(phi, na.rm=TRUE), 
  col = spectralcolors(length(breaks)-1),   # colors to use; by default blue to red
  add = FALSE, xaxs = "i", yaxs = "i",  
  probs = seq(0,1,0.1),
  breaks=quantile(z[is.finite(z)], probs=probs),  ...  # elevation cuts to use; by default, deciles
){
  #requireNamespace("DescTools", quietly = TRUE)
  # create a plotting region
  if(!add) plot(c(-1,1)*rlim[2], c(-1,1)*rlim[2],type="n", asp=1, xlab="", ylab="")
  # define the grid of polar coordinate indices
  polargrid = expand.grid(1:(length(r)-1), 1:(length(phi)-1))
  # choose the color of each sector
  mycol <- col[cut(c(z),breaks=breaks)]
  dim(mycol) = dim(z)
  # check which is the direction of the series of angles
  ss = sign(diff(phi))
  if(!all(ss>0) & !all(ss<0)) stop("neither clockwise nor counterclockwise angular series")
  # draw counterclockwise
  if(all(ss>0)){
    res <- gmApply( polargrid,1,
                  function(i) gsi.DrawAnnulusSector(x = 0, y = 0, radius.in = r[i[1]], 
                                                    radius.out = r[i[1]+1],
                                                    angle.beg = phi[i[2]], 
                                                    angle.end = phi[i[2]+1], 
                                                    col = mycol[i[1],i[2]],
                                                    border=NA)
    )}
  # draw clockwise  
  if(all(ss<0)){
    res <- gmApply( polargrid,1,
                  function(i) gsi.DrawAnnulusSector(x = 0, y = 0, radius.in = r[i[1]], radius.out = r[i[1]+1],
                                                    angle.end = phi[i[2]], angle.beg = phi[i[2]+1], col = mycol[i[1],i[2]],
                                                    border=NA)
    )}
  invisible(res)
}
#################################################

gsi.DrawAnnulusSector = function(x, y, radius.in, 
                                 radius.out,
                                 angle.beg, 
                                 angle.end, 
                                 ...){
  if(requireNamespace("DescTools", quietly=TRUE)){
    DescTools::DrawCircle(x=x, y=y, r.out = radius.out, r.in=radius.in, theta.1 = angle.beg,
                          theta.2 = angle.end, ...)
  # }else if(requireNamespace("circlize", quietly=TRUE)){
  #   circlize::draw.sector(
  #     start.degree = 0,
  #     end.degree = 360,
  #     rou1 = 1,
  #     rou2 = NULL,
  #     center = c(0, 0),
  #     clock.wise = TRUE,
  #     col = NA,
  #     border = "black",
  #     lwd = par("lwd"),
  #     lty = par("lty"))
  }else stop("DrawAnnulusSector: one of the packages 'Desctools' or 'circlize' is needed. Install one!")
  
}

gsi.DrawCircle = function(x, y, r.in, 
                          r.out,
                          theta.1, 
                          theta.2, 
                          ...){
  if(requireNamespace("DescTools", quietly=TRUE)){
    DescTools::DrawCircle(x=x, y=y, r.out = r.out, r.in=r.in, theta.1 = theta.1,
                          theta.2 = theta.2, ...)
  # }else if(requireNamespace("circlize", quietly=TRUE)){
  #   circlize::draw.sector(
  #     start.degree = theta.2,
  #     end.degree = theta.1,
  #     rou1 = r.out,
  #     rou2 = r.in,
  #     center = c(x, y),
  #     clock.wise = FALSE,
  #     ...)
  }else stop("DrawCircle: one of the packages 'Desctools' or 'circlize' is needed. Install one!")
}


#### general internal functions -------------------------------

gsi.powM <- function(A,alpha=1/2) {
  with(svd(A),u%*%diag(d^alpha)%*%t(v))
}



gsiDiag <- function(d) {
  if( length(d)>1)
    return(diag(d))
  return( structure(d,dim=c(1,1)))
}

gsiInv <- function(A,tol=1E-15) {
  with(svd(A),v%*% gsiDiag(ifelse(abs(d)/max(d)>tol,1/d,0))%*%t(u))
}

#### commons for C-R mixed code ---------------------------

checkInt <- function(x,n) {
  x <- as.integer(x)
  if( !missing(n) ) stopifnot(length(x)==n)
  x
}

checkDouble <- function(x,n) {
  x <- as.numeric(x)
  if( !missing(n) ) stopifnot(length(x)==n)
  x
}

checkLogical <- function(x,n) {
  x <- as.logical(x)
  if( !missing(n) ) stopifnot(length(x)==n)
  x
}
