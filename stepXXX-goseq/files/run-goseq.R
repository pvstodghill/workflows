#! /usr/bin/env Rscript
# -*- r -*-

# $ setup-bioc goseq

options(width=10000) # Pring wide lines

## ------------------------------------------------------------------------
## Parse command line
## ------------------------------------------------------------------------

## Skip the cruft at the beginning of the arg list
lst <- commandArgs()
while (length(lst) > 0) {
    if (substr(lst[1],1,7) == "--file=") {
        programName <- substr(lst[1],8,nchar(lst[1]))
    } else if (lst[1] == "--args") {
        lst <- lst[-1] # shift
        break
    }
    lst <- lst[-1] # shift
}

## parse into `args`
args <- list()
args$help <- FALSE
args$pvalue <- 1
args$qvalue <- 1
args$.error <- FALSE
args$.rest <- c()
while (length(lst)) {
    if (substr(lst[1],1,1) != '-') {
        args$.rest <- lst
        break
    } else if (lst[1] == "--") {
        lst <- lst[-1] # shift
        args$.rest <- lst
        break
    } else if (lst[1] == "-h" || lst[1] == "--help") {
        lst <- lst[-1] # shift
        args$help = TRUE;
    } else if (lst[1] == "-p" || lst[1] == "--pvalue") {
        if ( length(lst) == 1) {
            cat("Missing argument to",lst[1],"\n\n")
            args$.error <- TRUE
            break
        }
        args$pvalue <- as.numeric(lst[2])
        lst <- lst[-1] # shift
        lst <- lst[-1] # shift
    } else if (lst[1] == "-q" || lst[1] == "--qvalue") {
        if ( length(lst) == 1) {
            cat("Missing argument to",lst[1],"\n\n")
            args$.error <- TRUE
            break
        }
        args$qvalue <- as.numeric(lst[2])
        lst <- lst[-1] # shift
        lst <- lst[-1] # shift
    } else {
        cat("Unknown option:",lst[1],"\n\n")
        args$.error <- TRUE
        break
    }
}
remove(lst)

if (length(args$.rest) != 4) {
    cat("Expected exactly 4 positional arguments\n");
    args$.error <- TRUE
}

if (args$help || args$.error) {
    cat("Usage:",programName,"[options] arg assayed.txt len.txt regulon.txt funcs.txt\n\n")
    cat("-h, --help             this message\n")
    cat("-p NUM, --pvalue NUM   use NUM as p-value cutoff (default: 1.0)\n")
    cat("-q NUM, --qvalue NUM   use NUM as q-value cutoff (default: 1.0)\n")
    q(save="no")
}

assayed.txt <- args$.rest[1]
len.txt <- args$.rest[2]
regulon.txt <- args$.rest[3]
funcs.txt <- args$.rest[4]

cutoff.qval <- args$qvalue
cutoff.pval <- args$pvalue

## ------------------------------------------------------------------------
## Set up the data frame for the analysis
## ------------------------------------------------------------------------

assayed.genes <- scan(assayed.txt, what=character(), quiet=TRUE)
assayed.lens <-  scan(len.txt, what=numeric(), quiet=TRUE)

regulon.genes <- scan(regulon.txt, what=character(), quiet=TRUE)

genes.vector <- as.integer(assayed.genes %in% regulon.genes)
names(genes.vector) <- assayed.genes

## ------------------------------------------------------------------------
## Do the analysis
## ------------------------------------------------------------------------

suppressPackageStartupMessages(library("goseq"))

pwf <- nullp(genes.vector,bias.data=assayed.lens,plot.fit=FALSE)

terms <- read.table(funcs.txt,header=FALSE,sep="\t")

## use_genes_without_cat=TRUE is a design decision that might have to be revisited.
## suppressMessages(
results <- goseq(pwf,gene2cat=terms,use_genes_without_cat=TRUE)
## )

## ------------------------------------------------------------------------
## Filter, reorder, and munge the output
## ------------------------------------------------------------------------

results$over_represented_qvalue <- p.adjust(results$over_represented_pvalue,method="fdr")
if (cutoff.qval < 1 ) {
    results=results[results$over_represented_qvalue<cutoff.qval,]
}
if (cutoff.pval < 1 ) {
    results=results[results$over_represented_pvalue<cutoff.pval,]
}

output.columns <- c("category", "numDEInCat", "numInCat",
                    "over_represented_pvalue", "over_represented_qvalue")
if ( "term" %in% colnames(results) ) {
     output.columns <- c(output.columns,"term")
}
write.table(results[order(results$over_represented_pvalue),output.columns],
            sep="\t",row.names=FALSE)
