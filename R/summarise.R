#' @title Plot LogLikelihood
#' @aliases plotlike
#' @description  A log likelihood trace plot of hapsampler output.
#' @param  x a \code{HS} object, obtained from running \code{\link{hapsampler}}.
#'  @param ...   additional plot options. (not yet implemented).
#' @details
#' Produces a plot of the change in log likelihood for each chain, over 
#' the run length.  It is useful for looking for agreeement between chains and also 
#' that the log likelihoods have stabilized over a run. 
#' @return
#'   a plotting window is opened.
#' @seealso  \code{\link{plotEffects}}
#' @export
plotLike <- function(x)
{
  UseMethod("plotLike")
}

#' @title Plot LogLikelihood
#' @aliases plotlike.HS
#' @description  A log likelihood trace plot of hapsampler output.
#' @param  x a \code{HS} object, obtained from running \code{\link{hapsampler}}.
#'  @param ...   additional plot options. (not yet implemented).
#' @details
#' Produces a plot of the change in log likelihood for each chain, over 
#' the run length.  It is useful for looking for agreeement between chains and also 
#' that the log likelihoods have stabilized over a run. 
#' @return
#'   a plotting window is opened.
#' @seealso  \code{\link{plotEffects}}
#' @export
plotLike.HS <- function(x ) {

# to be consistent with John's original naming of variable
canon <- x[["canonical.hap"]]


  # Plot variables
  w <- 1024
  h <- 1024
  q <- 100
  ho <- FALSE
  ps <- 24
  fm <- 2
  
  chain <- NULL
  for (repl in 1:x[["nchains"]]) {
    chain <- rbind(chain,
	data.frame(
		rep=repl,
		read.table(paste("Temp/chain_",repl,
					".txt",sep=""),header=T)))
  }

  #######################################
  # Log likelihood
  #######################################

#  jpeg(filename = "/loglik.jpg""), 
#  width = 1*w, height = 1*h,
#  pointsize = ps, quality = q, bg = "white", res = NA)

#  par(mfrow=c(1,1))

  loglik <- chain$loglik

  maxlik <- max(loglik) + 100
  minlik <- max((min(loglik) - 100),3*maxlik)

  loglik[loglik < minlik] <- minlik
  
  leg.txt <- NULL
  leg.col <- NULL

  for (repl in 1:x[["nchains"]]) {

    logfile <- loglik[chain$rep == repl]

    leg.txt <- c(leg.txt,paste("chain ",repl,sep="")) 
    leg.col <- c(leg.col,repl+1) 
  
    if (repl == 1) {
  
      plot(logfile,type="l",col=(repl+1),
	  #xlim = c(0,(nburn+nsamp)),
	  ylim = c(minlik,maxlik),
	  #main = "Chains",
	  xlab = "sample",
	  ylab = "log likelihood")
  
    } else {
      lines(logfile,type="l",col=(repl+1))
    }
  }

  lines(c(x[["nburn"]],x[["nburn"]]),c(-9999,0))
  
  legend("bottomright", inset=0.05, leg.txt,fill=leg.col,horiz=ho,cex=0.8)

#  dev.off()

 }  ## end plot.HS


##--------------------------
#   plot of mu and sd
##--------------------------
#' @title Plot Effects
#' @aliases ploteffects
#' @description  A trace plot of the mean effect sizes of AA, AB, and BB from 
#'       the running of \code{\link{hapsampler}}.
#' @param  x a \code{HS} object, obtained from running \code{\link{hapsampler}}.
#' @details
#' Produces a plot of the change in mean effect size  for each chain, over 
#' the run length.  It is useful for looking for agreeement between chains and also 
#' that the mean effect size estimates are behaving. 
#' @return
#'   a plotting window is opened.
#' @seealso  \code{\link{plotLike}}
#' @export
plotEffects <- function(x )
{
  UseMethod("plotEffects")
}



#' @title Plot Effects
#' @aliases ploteffects.HS
#' @description  A trace plot of the mean effect sizes of AA, AB, and BB from 
#'       the running of \code{\link{hapsampler}}.
#' @param  x a \code{HS} object, obtained from running \code{\link{hapsampler}}.
#' @details
#' Produces a plot of the change in mean effect size  for each chain, over 
#' the run length.  It is useful for looking for agreeement between chains and also 
#' that the mean effect size estimates are behaving. 
#' @return
#'   a plotting window is opened.
#' @seealso  \code{\link{plotLike}}
#' @export
plotEffects.HS <- function(x){

# construct results object - chain
chain <- NULL
  for (repl in 1:x[["nchains"]]) {
    chain <- rbind(chain,
        data.frame(
                rep=repl,
                read.table(paste("Temp/chain_",repl,
                                        ".txt",sep=""),header=T)))
  }



  #######################################
  # Mu and Stdev
  #######################################

  lab <- c("AA","AB","BB")

#  fn <- paste(x[["trait"]],"/effects.jpg",sep="") 
#  jpeg(filename = fn, 
#  width = 3*w, height = 2*h,
#  pointsize = ps, quality = q, bg = "white", res = NA)

  par(mfrow=c(2,3))


  mu <- chain[,4:6]
  miny <- min(mu)
  maxy <- max(mu)

  # Identify those that need to be flipped
  flip <- array(FALSE,x[["nchains"]])

  for (repl in 1:x[["nchains"]]) {
    cf <- read.table(paste("Temp/hapsols_",repl,
                                        ".txt",sep=""),header=T)$P
    if (cf[x[["canonical.hap"]]]  < 0.5) flip[repl] <- T
  }

 for (i in 1:3) {

    leg.txt <- NULL
    leg.col <- NULL

    for (repl in 1:x[["nchains"]]) {

      j <- i

      if (flip[repl]) {j <- 4-i} # flip

      mmu <- mu[chain$rep == repl,]

      leg.txt <- c(leg.txt,paste("chain ",repl,sep=""))
      leg.col <- c(leg.col,repl+1)


      if (repl == 1) {

        plot(mmu[,j],type="l",col=(repl+1),
          #xlim = c(0,(nburn+nsamp)),
          ylim = c(miny,maxy),
          main = lab[i],
          xlab = "sample",
          ylab = "effect")

      } else {
        lines(mmu[,j],type="l",col=(repl+1))
      }
      lines(c(x[["nburn"]],x[["nburn"]]),c(-9999,9999))
    }
  }


  sd <- chain[,7:9]
  miny <- 0
  maxy <- max(sd)

  for (i in 1:3) {

    leg.txt <- NULL
    leg.col <- NULL

    for (repl in 1:x[["nchains"]]) {

      j <- i

      if (flip[repl]) {j <- 4-i} # flip

      msd <- sd[chain$rep == repl,]

      leg.txt <- c(leg.txt,paste("chain ",repl,sep=""))
      leg.col <- c(leg.col,repl+1)

      lab <- c("AA","AB","BB")

      if (repl == 1) {

        plot(msd[,j],type="l",col=(repl+1),
          #xlim = c(0,(nburn+nsamp)),
          ylim = c(miny,maxy),
          main = lab[i],
          xlab = "stdev",
          ylab = "effect")

      } else {
        lines(msd[,j],type="l",col=(repl+1))
      }
      lines(c(x[["nburn"]],x[["nburn"]]),c(-9999,9999))
    }
  }


#  dev.off()

  return(flip)


}




##--------------------------
#   Scatter plots 
##--------------------------
#' @title Haplotype and Allele Scatter Plots
#' @description  Various scatter plots involving haplotypes and observed alleles (ordered by mean Prob(A))
#' verse Prob(A). 
#' @param  x a \code{HS} object, obtained from running \code{\link{hapsampler}}.
#' @details This function opens a graphic window that contains three plots. The first 
#' plot is of the sampled haplotypes (ordered by Prob(A)) against Prob(A).  The second plot 
#' is similar, except the scater plot is of the unique haplotypes, binned according to the levels 
#' given in the legend. The third plot is of the QTL alleles, ordered by mean P....
#' Produces a plot of the change in mean effect size  for each chain, over 
#' the run length.  It is useful for looking for agreeement between chains and also 
#' that the mean effect size estimates are behaving. 
#' @return
#'   a plotting window is opened.
#' @seealso  \code{\link{plotLike}}
#' @export
plotscatter  <- function(x )
{
  UseMethod("plotscatter")
}


#' @title Haplotype and Allele Scatter Plots
#' @description  Various scatter plots involving haplotypes and observed alleles (ordered by mean Prob(A))
#' verse Prob(A). 
#' @param  x a \code{HS} object, obtained from running \code{\link{hapsampler}}.
#' @details This function opens a graphic window that contains three plots. The first 
#' plot is of the sampled haplotypes (ordered by Prob(A)) against Prob(A).  The second plot 
#' is similar, except the scater plot is of the unique haplotypes, binned according to the levels 
#' given in the legend. The third plot is of the QTL alleles, ordered by mean P....
#' Produces a plot of the change in mean effect size  for each chain, over 
#' the run length.  It is useful for looking for agreeement between chains and also 
#' that the mean effect size estimates are behaving. 
#' @return
#'   a plotting window is opened.
#' @seealso  \code{\link{plotLike}}
#' @export
plotscatter.HS <- function(x) {
  # Identify those that need to be flipped
  flip <- array(FALSE,x[["nchains"]])

  for (repl in 1:x[["nchains"]]) {
    cf <- read.table(paste("Temp/hapsols_",repl,
                                        ".txt",sep=""),header=T)$P
    if (cf[x[["canonical.hap"]]]  < 0.5) flip[repl] <- T
  }



par(mfrow=c(2,2))


  # First need haplotype counts

  anis = read.table("Temp/anisols_1.txt",header=T)
  haps <- c(anis$hap1,anis$hap2)

  nhap <- max(haps)

  ct <- array(NA,nhap)
  for (i in 1:nhap) {
    ct[i] <- sum(haps == i)
  }

  repl <- 1
  csum <- data.frame(
    haplotype = read.table(paste("Temp/hapsols_",repl,
                                        ".txt",sep=""),header=T)$haplotype)
  chain <- NULL
  for (repl in 1:x[["nchains"]]) {
    cf <- read.table(paste("Temp/hapsols_",repl,
                                        ".txt",sep=""),header=T)$P
    if (flip[repl]) {cf <- 1 - cf}
    chain <- cbind(chain,cf)
  }

  csum <- cbind(csum,
    data.frame(
        count = ct,
        PA = rowMeans(chain)))


  # Plots
  # Only observed haplotypes
  ohp <- csum
  ohp <- cbind(ohp,chain)
  ohp <- subset(ohp,ohp$count > 0)
  nhap <- dim(ohp)[1]

  # Order probabilities for plotting
  ohp <- ohp[order(ohp$PA),]
  ohp <- cbind(ohp,xx = 1:nhap)

  pP <- ohp[,4:(x[["nchains"]]+3)]


  for (repl in 1:x[["nchains"]]) {
    if (repl == 1) {
      if (x[["nchains"]] > 1) {
        plot(pP[,repl],col=repl+1,pch=20,cex=0.6,
          xlim = c(0,nhap),
          ylim = c(0,1),
          xlab = "Haplotype (ordered by mean Prob(A))",
          ylab = "Prob(A)")
      } else {
        plot(pP,col=repl+1,pch=20,cex=0.6,
          xlim = c(0,nhap),
          ylim = c(0,1),
          xlab = "Haplotype (ordered by mean Prob(A))",
          ylab = "Prob(A)")
      }
    } else {
      points(pP[,repl],col=repl+1,pch=20,cex=0.6)
    }
  }


  #  Add line of means


  lines(ohp$xx,ohp$PA,col=2,lwd=2)


  #dev.off()


  # Now plot means, both raw and by count

#  jpeg(filename = paste(trait,"/chain_by_count.jpg",sep=""),
#        width = 2*w, height = 1*h,
#        pointsize = ps, quality = q, bg = "white", res = NA)

#  par(mfrow=c(1,2))

  # Left panel is by haplotype allele

  ss <- ohp[ohp$count == 1,]
  plot(ss$xx,ss$PA,col=2,pch=20,cex=0.7,
        main = "Unique Haplotypes",
        xlim = c(0,nhap),
        ylim = c(0,1),
        xlab = "Haplotype (ordered by mean Prob(A))",
        ylab = "Prob(A)")

  ss <- subset(ohp,ohp$count == 2)
  points(ss$xx,ss$PA,col=3,pch=20,cex=0.6)

  ss <- subset(ohp,ohp$count > 2 & ohp$count <= 5)
  points(ss$xx,ss$PA,col=4,pch=20,cex=0.6)

  ss <- subset(ohp,ohp$count > 5 & ohp$count <= 10)
  points(ss$xx,ss$PA,col=5,pch=20,cex=0.6)

  ss <- subset(ohp,ohp$count > 10 & ohp$count <= 20)
  points(ss$xx,ss$PA,col=6,pch=20,cex=0.6)

  ss <- subset(ohp,ohp$count > 20)
  points(ss$xx,ss$PA,col=7,pch=20,cex=0.6)

  leg.txt <- c(
        "1",
        "2",
        "3-5",
        "6-10",
        "11-20",
        ">20")
  leg.col <- c(2,3,4,5,6,7)
  legend("topleft", inset=0.05, leg.txt,fill=leg.col,horiz=FALSE,cex=0.8)

  # In right panel the x axis is by times seen

  full.ohp <- NULL
  for (i in 1:nhap) {
    ct <- round(ohp$count[i])
    for (j in 1:ct) {full.ohp <- rbind(full.ohp,ohp[i,])}
  }

  nfhap <- dim(full.ohp)[1]
  full.ohp$xx <- 1:nfhap

  plot(full.ohp$xx,full.ohp$PA,col=2,pch=20,cex=0.6,
        main = "Observed alleles",
        xlim = c(0,nfhap),
        ylim = c(0,1),
        xlab = "Allele (ordered by mean Prob(A))",
        ylab = "Prob(A)")



#  ohp <- ohp[,1:3]
#  ohp <- ohp[order(ohp$haplotype),]

#  ohp

}







