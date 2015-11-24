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

To install the package, directly from github, use the following commands:

```R
install.packages("devtools")
devtools::install_github("geo047/HapSampler")
```



 

