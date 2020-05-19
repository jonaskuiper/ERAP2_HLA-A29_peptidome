###
### Sample code to reproduce the NMDS plots for the 948 9-mers from the HLA-A29 peptidome
### This code is adapted from  Sarkizova et al Nature Biotechnology volume 38, pages199–209(2020) 

## Setup dependencies ----
library(HDMD)
library(ecodist)
### END: Setup dependencies ----

### Load  data ----
rndseed = 12311
load("HLA_A29_948_9mers.RData")
peps<-HLA_A29_9_mers
npeps <- length(peps)
len=9
### Calculate Molecular Entropy for amino acid positions ----
molecularEntropy <- MolecularEntropy(peps, type='AA')

#NMDS 
### NMDS plot ----
### Compute peptide distances using amino acid substitutionmatrix "distPMBEC file" 
distPMBEC <- read.table('./distPMBEC.txt', header=TRUE, stringsAsFactors=FALSE)

# functions
getDistances <- function(peps, pos_weights=NULL, len=9) {
        if (is.null(pos_weights)) { pos_weights = rep(1, len) }
        n <- length(peps)
        dists <- matrix(nrow=n, ncol=n)
        . <<- lapply(1:n, function(i) { dists[i,i:n] <<- unlist(lapply(i:n, function(j) {
                pepDistPMBEC(strsplit(peps[i], '')[[1]], strsplit(peps[j], '')[[1]], len, pos_weights)
        }))})
        dists[lower.tri(dists)] <- t(dists)[lower.tri(dists)]
        colnames(dists) <- rownames(dists) <- peps
        return(dists)
}

pepDistPMBEC <- function(pepA, pepB, len, pos_weights) {
        return ((sum(unlist(lapply(1:len, function(i) {
                distPMBEC[pepA[i], pepB[i]]*pos_weights[i]
        })))) / len)
}

pos_weights <- 1 - molecularEntropy$H
dists <- getDistances(peps, pos_weights)

### Reduce dimensions (NMDS)
### Note: This step is computationally intensive, it will run relatively
###       quickly for a few hundred peptides but with a few thousand
###       peptides it will take a while.
###       nits=10 is used to generate plots with maxit=500. The trace will help you track progress
set.seed(rndseed) #this is important. Keep it the same
nmds.tmp <- nmds(as.dist(dists), mindim=2, maxdim=2, nits=10, maxit=500, trace=TRUE)
nmds.dat <- nmds.tmp$conf[[which.min(nmds.tmp$stress)]]

### Plot
pdf('./nmds.pdf', width=4.5, height=5)
plot(nmds.dat, xlab='nmds1', ylab='nmds2', col="grey")
mtext(side=3, line=0, text="A2902_NMDS")
dev.off()
### END: NMDS plot ----
















