\name{HapSampler-package}
\alias{HapSampler-package}
\alias{HapSampler}
\docType{package}
\title{
A implementation  for constructing  joint marker-QTL haplotypes and 
their probability 
}
\description{
An MCMC implementation for constructing joint marker-QTL haplotypes and for 
approximating their probability. The input data consists of 
the trait data  and the observed marker haplotypes  on a set of animals. 
The trait data can be continuous or discrete.  Functions are available for 
reading in the data, for analysing the data via MCMC, and for plotting the results from the MCMC runs. 
}
\details{
\tabular{ll}{
Package: \tab HapSampler\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-11-2\cr
License: GPL (>2)
}
}
\author{
John Henshall (original author)
Andrew George (developer of package)

Maintainer: Andrew George <andrew.george@csiro.au>
}
\keyword{HapSampler}
\examples{
# a worked example

# read in data from dataexample.dat.
# The file is organised where each row is the data for a different animal.
# The data file must contain a column with the trait data, 
# two columns with the haplotype indexes and 
# a column for the probability of observing these haplotypes in the population. The 
# haplotypes  and probability are obtained by an external phasing package such as BEAGLE.
complete.name <- system.file("extdata", "dataexample.dat", package="HapSampler")
dt <- read_data(path=dirname(complete.name),
                file=basename(complete.name))
     
# print the first couple of rows of the data 
head(dt)  

# construct the joint QTL-marker haplotypes and their probabilitie via MCMC.
# Here, we will run 4 parallel chains, with each chain being of run length 5.
# In practice, you would have longer run lengths.  
hapres <- hapsampler(data=dt, trait="NAVEL", nchains=4, runlength=5)


# plot the log likelihoods of these four chains over their run lengths.  We are looking
# for a stabilization within and across chains of the log likelihoods. 
plotLike(hapres)

# plot the sampled haplotypes and QTL alleles against the mean haplotype probability
# The number of bins for the top left plot has been set to 20. See help(plotScater) for details 
plotScatter(hapres, nbins=20)


}
