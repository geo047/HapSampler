\name{HapSampler-package}
\alias{HapSampler-package}
\alias{HapSampler}
\docType{package}
\title{
An MCMC implementation for constructing joint QTL-marker haplotypes
}
\description{
An MCMC implementation for constructing joint QTL-marker haplotypes. It builds these 
haplotypes from phenotypic and marker haplotype data. The phenotypic data must include
a trait which may be either continuous or binary. Functions are available for reading 
in the data, analysing the data, and for summariing the results. 
}
\details{
\tabular{ll}{
Package: \tab HapSampler\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-11-2\cr
License: \tab What license is it under?\cr
}
Overview of what package does
}
\author{
John Henshall (original author)
Andrew George (developer of package)

Maintainer: Andrew George <andrew.george@csiro.au>
}
\keyword{HapSampler}
\examples{
# a worked example

# read in the phenotypic and haplotype  data
complete.name <- system.file("extdata", "dataexample.dat", package="HapSampler")
# read in data from file dataexample.dat. This file contains the phenotypic and 
# haplotype data. 
dt <- read_data(path=dirname(complete.name),
                file=basename(complete.name))
     
# to print out the first couple of rows of the data now read into R
head(dt)  

# construct the joint QTL marker haplotypes via MCMC methods.
# Here, we will run 4 parallel chains, with each chain being of run length 20
hapres <- hapsampler(data=dt, trait="NAVEL", id="Animal_id", breed="Breed",
                              nchains=4, runlength=20)


# plot the log likelihoods of these four chains over their run lengths.  We are looking
# for a stabilization with and across chains of the log likelihoods. 
plotLike(hapres)

# plot the sampled haplotypes and QTL alleles against the mean haplotype probability
# The number of bins for the top left plot has been set to 20.  
plotscatter(hapres, nbins=20)


}
