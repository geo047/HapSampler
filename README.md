# HapSampler

An R  package for building  joint marker and qtl haplotypes via Markov chain Monte Carlo methods.
The package requires phenotypic and haplotype data.  The phenotypic data must contain a 
trait (either continuous or binary), and an animal id.  The haplotype data are the maternal 
and paternal haplotypes. These haploytpes can be  obtained from unphased genotype data via
phasing programs such as BEAGLE and fastPAHSE. 

Once the package has been installed (see below), a listing of its functions can be 
obtained from within your R session via the command:

```R
library(, HapSampler)
```

and documentation on the running of any command can be seen with the help command. For 
example, to obtain help on the hapsampler command, use

```R
help(hapsampler)
```


## Installation

This package is best installed from within R.

To install the package, directly from github, you need the "devtools" package which can be installed with

```R
install.packages("devtools")
```

Then, to install the HapSampler package, run the commands

```R
library("devtools")
devtools::install_github("geo047/HapSampler")
```

Suppose, you would like to install "HapSampler" into a local directory, say /home/id/RLib 

Then, you just need to add a ".libPaths" command as follows

```R
library("devtools")
.libPaths(c(.libPaths(), "/home/id/RLib"))
devtools::install_github("geo047/HapSampler")
```

 

