#' Calculate Jost's D
#'
#' This function calculates Jost's D from a genind object
#'
#' Takes a genind object with population information and calculates Jost's D
#' Returns a list with values for each locus as well as two global estimates.
#' 'global.het' uses the averages of Hs and Ht across all loci while
#' 'global.harm_mean' takes the harmonic mean of all loci.
#'
#' Because estimators of Hs and Ht are used, its possible to have negative
#' estimates of D. You should treat these as numbers close to zero.
#'
#' @return per.locus values for each D for each locus in the dataset
#' @return global estimtes for D based on overall heterozygosity or the harmonic
#' mean of values for each locus
#' @param hsht_mean The type of mean to use to calculate values of Hs and Ht
#' for a global estimate. (Default is teh airthmetic mean, can also be set to
#' the harmonic mean).
#' @param x freqency table in dataframe format. First column must contain population names, second number of samples per population.
#' Rest of the columns contain alle frequencies and must have names formatted like "locusname_allelename".
#' @export
#' @examples
#' @references
#'  Jost, L. (2008), GST and its relatives do not measure differentiation. Molecular Ecology, 17: 4015-4026.
#' @family diffstat
#' @family D


D_Jost.Allo <- function(x, hsht_mean = "arithmetic"){
  freqs<-Allo.format(x)
  nbp<-x[,2]
  names(nbp)<-x[,1]
  mean_type <- match.arg(hsht_mean, c("arithmetic", "harmonic"))
  mean_f <- if(mean_type == "arithmetic") mean else harmonic_mean
  gn <- length(nbp)
  loci <- t(sapply(freqs, D.per.locus.Allo, nbp=nbp))
  global_Hs <- mean_f(loci[,1], na.rm=T)
  global_Ht <- mean_f(loci[,2], na.rm=T)
  global_D <-  (global_Ht - global_Hs)/(1 - global_Hs ) * (gn/(gn-1))
  harm_D <- harmonic_mean(loci[,3])
  return(list("per.locus"=loci[,3],
              "global.het"=global_D,
              "global.harm_mean" = harm_D
  ))

}

Allo.format<-function(x){
  locs<-colnames(x)
  locs<-locs[3:length(locs)]
  locs<-gsub("_.*","",locs)
  locs%<>%unique()
  locsep<-list()
  for(i in 1:length(locs)){
    currl<-select(x, starts_with(locs[i]))
    #currl<-select(currl, !contains("null"))
    rownames(currl)<-freqs[,1]
    locsep[[i]]<-currl
  }
  names(locsep)<-locs
  return(locsep)
}

HtHsAllo<-function(x, nbp){
  HpS <- mean(1 - rowSums(x^2))
  Hs_est <- (2*harmN/(2*harmN-1))*HpS
  HpT <- 1 - sum(colMeans(x)^2)
  Ht_est <- HpT + Hs_est/(2*harmN*n)
  return(c(Ht_est, Hs_est))
}


D.per.locus.Allo <- function(g,nbp) {
  hets <- HtHsAllo(g, nbp)
  if(all(is.na(hets))){
    return(hets)
  }
  Ht_est <- hets[[1]]
  Hs_est <- hets[[2]]
  n <- length(nbp)
  D <- (Ht_est-Hs_est)/(1-Hs_est) * (n/(n-1))
  return(c(Hs_est, Ht_est, D))
}


