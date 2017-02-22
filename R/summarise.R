#' @title Plot LogLikelihood
#' @aliases plotlike
#' @description  A log likelihood trace plot of hapsampler output.
#' @param  x a \code{HS} object, obtained from running \code{\link{hapsampler}}.
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
#' @description  A trace plot of the mean qtl effect sizes of QQ, Qq, and qq from 
#'       the running of \code{\link{hapsampler}}.
#' @param  x a \code{HS} object, obtained from running \code{\link{hapsampler}}.
#' @details
#' Produces a plot of the change in mean qtl effect size  for each chain, over 
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



#' @title Mean and Standard Deviation Plots of QTL  Effects
#' @aliases ploteffects.HS
#' @description  A trace plot of the mean and standard deviation of 
#' qtl effect sizes of QQ, Qq, and qq from 
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

  lab <- c("QQ","Qq","qq")

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
          xlab = "Mean",
          ylab = "QTL genotype effect size")

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


      if (repl == 1) {

        plot(msd[,j],type="l",col=(repl+1),
          #xlim = c(0,(nburn+nsamp)),
          ylim = c(miny,maxy),
          main = lab[i],
          xlab = "Standard Deviation",
          ylab = "QTL genotype effect size")

      } else {
        lines(msd[,j],type="l",col=(repl+1))
      }
      lines(c(x[["nburn"]],x[["nburn"]]),c(-9999,9999))
    }
  }


}




##--------------------------
#   Scatter plots 
##--------------------------
#' @title Haplotype and QTL Allele Scatter Plots
#' @aliases plotScatter
#' @description  Various scatter plots involving the sampled haplotypes and sampled QTL alleles. 
#' @param  x a \code{HS} object, obtained from running \code{\link{hapsampler}}.
#' @param  nbins the number of bins for grouping the frequency each of the haploytpes are 
#' sampled. 
#' @details  A single graphic window is open, containing three plots. 
#'
#' The top left  plot is of the sampled haplotypes within a chain verse the probability of that 
#' haplotype carrying the Q  allele.  The sammpled haplotypes are ordered based on their average 
#' Q allele probability (Prob(Q)).  The results from each chain are coloured 
#' differently.  This plot is useful for identifying those marker-QTL haplotypes that have been sampled differently 
#' across the parallel chains  (i.e. we don't want points a long way away from the red fitted line). 
#'
#'
#' The top right plot is of the sampled haplotypes across chains verse the probability of that 
#' haplotype carrying QTL allele Q. The sampled haplotypes are colour coded, depending upon 
#' their frequency of being sampled across all n chains. The frequency ranges and colours, for the bins,
#' are given in the legend. The number of bins is adjusted by setting \code{nbins}. 
#' Bin ranges with 0 frequency are not listed in the legend. 
#' 
#' The bottom left plot is a scatter plot of the sampled QTL allele Q verse 
#' its probability.  
#' @return
#'   a plotting window is opened.
#' @seealso  \code{\link{plotLike}}
#' @export
plotScatter  <- function(x, nbins=10 )
{
  UseMethod("plotScatter")
}


#' @title Haplotype and QTL Allele Scatter Plots
#' @aliases plotScatter.HS
#' @description  Various scatter plots involving the sampled haplotypes and sampled QTL alleles. 
#' @param  x a \code{HS} object, obtained from running \code{\link{hapsampler}}.
#' @param  nbins the number of bins for grouping the frequency each of the haploytpes are 
#' sampled. 
#' @details  A single graphic window is open, containing three plots. 
#'
#' The top left  plot is of the sampled haplotypes within a chain verse the probability of that 
#' haplotype carrying the Q  allele.  The sammpled haplotypes are ordered based on their average 
#' Q allele probability (Prob(Q)).  The results from each chain are coloured 
#' differently.  This plot is useful for identifying those marker-QTL haplotypes that have been sampled differently 
#' across the parallel chains  (i.e. we don't want points a long way away from the red fitted line). 
#'
#'
#' The top right plot is of the sampled haplotypes across chains verse the probability of that 
#' haplotype carrying QTL allele Q. The sampled haplotypes are colour coded, depending upon 
#' their frequency of being sampled across all n chains. The frequency ranges and colours, for the bins,
#' are given in the legend. The number of bins is adjusted by setting \code{nbins}. 
#' Bin ranges with 0 frequency are not listed in the legend. 
#' 
#' The bottom left plot is a scatter plot of the sampled QTL allele Q verse 
#' its probability.  
#' @return
#'   a plotting window is opened.
#' @seealso  \code{\link{plotLike}}
#' @export
plotScatter.HS <- function(x, nbins=10) {
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
        main="Sampled QTL-marker haplotypes",
          xlim = c(0,nhap),
          ylim = c(0,1),
          xlab = "Haplotype (ordered by mean Prob(Q))",
          ylab = "Prob(Q)")
      } else {
        plot(pP,col=repl+1,pch=20,cex=0.6,
        main="Sampled QTL-marker haplotypes",
          xlim = c(0,nhap),
          ylim = c(0,1),
          xlab = "Haplotype (ordered by Prob(Q))",
          ylab = "Prob(Q)")
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


  h <- hist(ohp$count, plot=FALSE, breaks = nbins)
  breaks <- h$breaks
  counts <- h$counts
  breaks[counts>0]

  ## set up initial plot
  ss <- subset(ohp, ohp$count >= breaks[1] & ohp$count <= breaks[2])
  indx <- sample(2:657,nbins,FALSE)
  leg.col <- colors()[indx[1]]
  plot(ss$xx,ss$PA,col=leg.col,pch=20,cex=0.7,
        main="Sampled QTL-marker haplotypes",
        xlim = c(0,nhap),
        ylim = c(0,1),
        xlab = "Haplotype (ordered by mean Prob(Q))",
        ylab = "Prob(Q)")

  for(ii in 2:length(counts)){
     if(counts[ii] > 0){
        ss <- subset(ohp, ohp$count > breaks[ii] & ohp$count <= breaks[ii+1])
        points(ss$xx,ss$PA,col=colors()[indx[ii]],pch=20,cex=0.6)
        leg.col <- c(leg.col, colors()[indx[ii]])
     }
  }

  leg.txt <- vector("character", 0)
  for(ii in 1:length(counts)){
     if(counts[ii] > 0){
       leg.txt <- c(leg.txt, paste(breaks[ii], "-", breaks[ii+1], sep=" "))
     }
  }
  
  legend("topright", inset=0.05, leg.txt,fill=leg.col,horiz=FALSE,cex=0.8)





  # In right panel the x axis is by times seen

  full.ohp <- NULL
  for (i in 1:nhap) {
    ct <- round(ohp$count[i])
    for (j in 1:ct) {full.ohp <- rbind(full.ohp,ohp[i,])}
  }

  nfhap <- dim(full.ohp)[1]
  full.ohp$xx <- 1:nfhap

  plot(full.ohp$xx,full.ohp$PA,col=2,pch=20,cex=0.6,
        main = "Sampled alleles",
        xlim = c(0,nfhap),
        ylim = c(0,1),
        xlab = "Allele (ordered by mean Prob(Q))",
        ylab = "Prob(Q)")




}







