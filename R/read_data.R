#' @title read_data  
#' @aliases read_data
#' @description Read in data file.
#' @param path  a character vector containing the absolute path for where the
#'          phenotype file is located.  The default is that the data file
#'          is contained in the directory from which R is being run.
#'
#' @param    file a character vector containing the name of the data file.
#' @param     csv a logical value. When \code{TRUE}, a csv file format is assumed.
#'          When \code{FALSE}, a space separated format is assumed.  The
#'          default is for the file to be space separaed.
#' @details Here, \code{\link{read_data}} reads in the ASCII data into the package. A
#'     space separated ASCII file with column headings is assumed. A csv
#'     file can also be read if \code{csv} is set to \code{TRUE}. Each row of the
#'     data file is the data for a seperate animal.  The data file is
#'     allowed to have any number of columns. 
#'
#'     Three columns are mandatory in order for ‘hapsampler’ to function
#'     correctly; the two haplotype indexes and a haplotype probability.
#'     These three columns must be named hap1, hap2, and prob where prob
#'     is the haplotype probability of the first haplotype (hap1). The
#'     function will stop if the column names hap1, hap2, and prob are
#'     not found.
#'
#'     The columns hap1 and hap2 contain the haplotype indexes. These
#'     indexes can be any integer number.
#'
#'     NOTE WHAT about missing data???
#'
#' @return 
#'     a data frame is returned of the data. The names of the columns
#'     will be as specified by the first row of the data file.
#'
#'@examples 
#'     # Read in  example phenotypic data from ./inst/extdata/
#'     
#'     # finding the full location of the phenotypic data
#'     complete.name <- system.file("extdata", "dataexample.dat", package="HapSampler")
#'     
#'     # read in phenotypic data which is in csv format
#'     dt <- read_data(path=dirname(complete.name),
#'                                  file=basename(complete.name))
#'     
#'      ## print a couple of lines of the data file
#'      head(dt)
#' @export
 read_data <- function (path = getwd(), file = NULL, csv = FALSE) 
{
    if (.Platform$OS.type == "unix") {
        dir_path <- paste(path, "/", sep = "")
    } else {
        dir_path <- paste(path, "\\", sep = "")
    }
    file_name <- paste(dir_path, file, sep = "")
    if (!file.exists(file_name)) {
        cat(" File ", file_name, " does not exist. \n")
        stop(" Please modify path and/or name of data file. \n")
    }
    if (is.null(file)) 
        stop(" No name of the data file has been supplied.")
    cf <- count.fields(file_name)
    if (length(unique(cf)) > 1) {
        cat(" File ", file, " has  differing numbers of fields per line ... \n")
        cat(" The differing number of fields  per line are ", 
            unique(cf), "\n")
        stop(" Please modify the genotype file to have the same number of genotypes per row. \n")
    }
    sep <- ""
    if (csv) 
        sep <- ","
    phdata <- read.table(file_name, sep = sep, header=TRUE)
    indx <- which(names(phdata) == "prob")
    if (length(indx) == 0) 
        stop(" Data file must contain a column heading called prob.")
    indx <- which(names(phdata) == "hap1")
    if (length(indx) == 0) 
        stop(" Data file must contain a column heading called hap1.")
    indx <- which(names(phdata) == "hap2")
    if (length(indx) == 0) 
        stop(" Data file must contain a column heading called hap2.")
    indx <- with(phdata, which(is.na(prob)))
    if (length(indx) > 0) 
        stop(" Data file must not contain missing haplotype probabilities. ")
    indx <- with(phdata, any(prob < 0 || prob > 1))
    if (indx) 
        stop(" Invalid haplotype probabilities have been found (values < 0 or > 1).")
    cat(" \n\n Loading data file ...  \n\n")
    cat(" SUMMARY OF DATA FILE  \n")
    cat(" file location(path):         ", dirname(file_name), 
        "\n")
    cat(" file name:                   ", basename(file_name), 
        "\n")
    cat(" number of rows:              ", nrow(phdata), "\n")
    cat(" number of columns:           ", ncol(phdata), "\n")
    cat("\n                    Column classes  \n")
    for (ii in 1:ncol(phdata)) cat(c(sprintf("%20s   %15s", names(phdata)[ii], 
        class(phdata[[ii]])), "\n"))
    cat("\n Warning: if the column classes are incorrect, these will need to be changed by the user.\n\n\n")
    return(phdata)
}



