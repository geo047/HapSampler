
# Haplotype association sampler

# John.Henshall@csiro.au
# Nov 2013

# Oct 2014, continuous version added 

# Fast version

#####################################################################


.forfortran <- function(ani.haplo,pheno,cand.pen,samp.pen,
					nburn,nsamp,nchains)  {

# Discrete trait version, especially for polled.

  # cand.pen is the probability of a P given phenotype and other allele
  # Row names are phenotypes (ptpp means progeny tested PP, etc)
  # Columns are other allele

  # Probability of an H is 1 - cand.pen
  
  # log them
  lPP <- log(cand.pen)
  lPH <- log(1-cand.pen)


  # samp.pen is the diploid penetrance, should be P(phen|gen) but this is 
  # other way around so rows sum to one. 

  # Needn't be totally consistent with cand.pen, as that is only the 
  # distribution, for generating proposals.

  # Rows are phenotype
  # ptpp means progeny tested PP

  # Columns are  number of copies of H, but this didn't store properly so
  # overwrite

  colnames(samp.pen) <- 0:2

  # log it
  ldpen <- log(samp.pen)

  haps <- c(ani.haplo$h1,ani.haplo$h2)
  nhap <- max(haps)

  # Keep only highly likely haplotypes
  anidat <- subset(ani.haplo, ani.haplo$prob > 0.99)

  # Add in phenotype
  anidat <- merge(anidat,pheno,by.x="sample_name",by.y="Sample_Name",all.x=T)

  # Angus, assumed PP so make make phenotype ptpp
  anidat$Phenotype <- as.vector(anidat$Phenotype)
  anidat$Phenotype[anidat$Breed=="Angus"] <- "ptpp"

  # Update phenotypes for progeny test sires
  anidat$Phenotype[anidat$Prog_Test=="PH"] <- "ptph"
  anidat$Phenotype[anidat$Prog_Test=="HH"] <- "pthh"
  anidat$Phenotype[anidat$Prog_Test=="PP"] <- "ptpp"
  

  anirec  <- anidat[,c(3,4,6)]

  nanis <- dim(anirec)[1]

  # Flag homozygous animals
  homoz <- anirec$h1 == anirec$h2
  anirec <- cbind(anirec,homoz)

  h1 <- anirec$h1
  h2 <- anirec$h2

  # Build penetrance arrays

  # Polled sampling array
  Psamp.array <- array(NA,c(nanis,4))

  # Horned sampling array
  Hsamp.array <- array(NA,c(nanis,4))

  # Likelihood penetrance array
  lik.array <- array(NA,c(nanis,3))

  for (ani in 1:nanis) {

    if (anirec$Phenotype[ani] == "horned") {
      Psamp.array[ani,] <- as.matrix(lPP[1,])
      Hsamp.array[ani,] <- as.matrix(lPH[1,])
      lik.array[ani,] <- as.matrix(ldpen[1,])
    }
    if (anirec$Phenotype[ani] == "polled") {
      Psamp.array[ani,] <- as.matrix(lPP[2,])
      Hsamp.array[ani,] <- as.matrix(lPH[2,])
      lik.array[ani,] <- as.matrix(ldpen[2,])
    }
    if (anirec$Phenotype[ani] == "scurred") {
      Psamp.array[ani,] <- as.matrix(lPP[3,])
      Hsamp.array[ani,] <- as.matrix(lPH[3,])
      lik.array[ani,] <- as.matrix(ldpen[3,])
    }
    if (anirec$Phenotype[ani] == "ptpp") {
      Psamp.array[ani,] <- as.matrix(lPP[4,])
      Hsamp.array[ani,] <- as.matrix(lPH[4,])
      lik.array[ani,] <- as.matrix(ldpen[4,])
    }
    if (anirec$Phenotype[ani] == "ptph") {
      Psamp.array[ani,] <- as.matrix(lPP[5,])
      Hsamp.array[ani,] <- as.matrix(lPH[5,])
      lik.array[ani,] <- as.matrix(ldpen[5,])
    }
    if (anirec$Phenotype[ani] == "pthh") {
      Psamp.array[ani,] <- as.matrix(lPP[6,])
      Hsamp.array[ani,] <- as.matrix(lPH[6,])
      lik.array[ani,] <- as.matrix(ldpen[6,])
    }
    if (anirec$Phenotype[ani] == "unknown") {
      Psamp.array[ani,] <- 0
      Hsamp.array[ani,] <- 0
      lik.array[ani,] <- 0
    }

  }

  Psamp.array[anirec$homoz,3] <- Psamp.array[anirec$homoz,4] 
  Psamp.array <- Psamp.array[,1:3] 
  
  Hsamp.array[anirec$homoz,3] <- Hsamp.array[anirec$homoz,4] 
  Hsamp.array <- Hsamp.array[,1:3] 

  # Write out input file

  ofn <- "TmpData/forfortran.txt"

  write(nchains,file=ofn, append = F)
  write(nburn,file=ofn, append = T)
  write(nsamp,file=ofn, append = T)
  write(nanis,file=ofn, append = T)

  h <- cbind(h1,h2)
  for (i in 1:nanis) {
    write(anidat$sample_name[i],file=ofn, append = T)
    write(as.character(anidat$Breed[i]),file=ofn, append = T)
    write(anidat$Phenotype[i],file=ofn, append = T)
    write(h[i,],file=ofn, append = T)
    write(Psamp.array[i,],file=ofn, append = T)
    write(Hsamp.array[i,],file=ofn, append = T)
    write(lik.array[i,],file=ofn, append = T)
  }


} 

##########################################################################


.forfortran.cont <- function(id,cphen,h1,h2,samp.pen.mu,samp.pen.sd,
			    hap.assign,nsamp)  {

# Continuous trait version

  nanis <- length(h1)

  # Flag homozygous animals
  homoz <- h1 == h2

  # Build penetrance arrays

  # Polled sampling array
  Psamp.array <- array(NA,c(nanis,4))

  # Horned sampling array
  Hsamp.array <- array(NA,c(nanis,4))

  # Likelihood penetrance array
  lik.array <- array(NA,c(nanis,3))

  likpp <- dnorm(cphen,samp.pen.mu[1],samp.pen.sd[1])
  likph <- dnorm(cphen,samp.pen.mu[2],samp.pen.sd[2])
  likhh <- dnorm(cphen,samp.pen.mu[3],samp.pen.sd[3])

  tot <- likpp+likph+likhh
  likpp <- likpp/tot
  likph <- likph/tot
  likhh <- likhh/tot

  likpp[is.na(likpp)] <- 1/3
  likph[is.na(likph)] <- 1/3
  likhh[is.na(likhh)] <- 1/3

  lik.array <- log(cbind(likpp,likph,likhh))

  # sampling arrays are functions of the likelihood array
  Psamp.array[,1] <- likpp/(likpp+likph)
  Psamp.array[,2] <- likph/(likph + likhh)
  Psamp.array[,3] <- likpp + likph/2
  Psamp.array[,4] <- likpp/(likpp + likhh)

  Hsamp.array <- 1 - Psamp.array
  Psamp.array <- log(Psamp.array)
  Hsamp.array <- log(Hsamp.array)

  Psamp.array[homoz,3] <- Psamp.array[homoz,4] 
  Psamp.array <- Psamp.array[,1:3] 
  
  Hsamp.array[homoz,3] <- Hsamp.array[homoz,4] 
  Hsamp.array <- Hsamp.array[,1:3] 

  # Write out input file

  ofn <- "TmpData/forfortran.txt"

  write(nsamp,file=ofn, append = F)
  write(nanis,file=ofn, append = T)

  h <- cbind(h1,h2)
  for (i in 1:nanis) {
    write(id[i],file=ofn, append = T)
    write(cphen[i],file=ofn, append = T)
    write(h[i,],file=ofn, append = T)
    write(Psamp.array[i,],file=ofn, append = T)
    write(Hsamp.array[i,],file=ofn, append = T)
    write(lik.array[i,],file=ofn, append = T)
  }

  # Write out starting values

  ofn <- "TmpData/starting.txt"

  l <- length(hap.assign)
  write(l,file=ofn, append = F)
  for (i in 1:l) {
    write(hap.assign[i],file=ofn, append = T)
  }

} 

