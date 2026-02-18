require(dplyr)
require(magrittr)
require(pegas)
require(PopGenReport)
require(hierfstat)
require(poppr)

#' A convinience function to check the usability of loci for population genetics analysis.
#'
#' Function checks the deviation from Hardy-Weinberg equilibrium, the significance of Fis deviation from zero and the presence of null-alleles for each supplied locus and population.
#' HW check and null-alleles are done using functions from packages pegas and PopGenREport. The Fis deviations is checked with Xi squared on the basis of Fis calculated with poppr.
#' If the locus shows deviation from HW, deviation of Fis from zero and null-alleles in half or more populations it's recommended to exclude it from population genetics analysis.
#'
#' @param genpop
#' Genpop file or a genind object containing genotypes
#' @param file
#' Set to TRUE if supllying genotypes in a file or to FALSE if in an already loaded genind object. Defaults to TRUE
#' @param OutputHW
#' Set to TRUE to output results of Hardy-Weinberg analysis in a file. Defaults to TRUE
#' @param OutputNull
#' Set to TRUE to output results of Null allele analysis in a file. Defaults to TRUE
#' @param OutputFis
#' Set to TRUE to output results of Fis analysis in a file. Defaults to TRUE
#' @param OutputPerLoc
#' Set to TRUE to output summary for each locus in a file. Defaults to TRUE
#' @param NullTreshold
#' Integer. What null allele frequency count as significant. Defaults to 0.05
#'
#' @returns
#' A list containting raw results for three analysis and two tables with summary statistics.
#' Quicksummary shows the number and names of populations for each locus in which the deviation by all three checks occurs.
#' Qualitytest shows extended summary
#' HW, NuAl and Fis contain the results for individual checks.
#' Optionally parts of the analysis are outputed in the .txt files
#' @export
#'
#' @examples
#' TBD
LocQual<-function(genpop, file=T, OutputHW=T, OutputNull=T, OutputFis=T, OutputPerLoc=T, NullTreshold=0.05){


if(file==T){
gens<-adegenet::read.genepop(genpop)
}else{gens<-genpop}
pops<-levels(gens@pop)
npops<-length(pops)
genssep<-adegenet::seppop(gens)
loci<-levels(gens@loc.fac)
nloc<-length(loci)

qualitytest<-list()
locqual<-as.data.frame(matrix(ncol = 5, nrow = npops))
colnames(locqual)<-c("Pop", "Hw_deviation", "Null_allelle", "Fis_significance", "overall")
locqual$Pop<-pops


for(i in 1:nloc){
  qualitytest[[i]]<-locqual
}
names(qualitytest)<-loci

message("Calculating HW")
HW<-list()
for(i in 1:npops){
  tempHW<-pegas::hw.test(genssep[[i]], B=1000)
  for(k in 1:nloc){
    qualitytest[[k]][i, "Hw_deviation"]<-tempHW[k,"Pr.exact"]
  }
  HW[[i]]<-tempHW
  message(paste0("Calculated HW for pop ", i, " out of ", npops,  sep=""))
}
names(HW)<-pops
message("Calculated HW for all pops")

if(OutputHW==T){
  cat("Deviation from Hardy-Weinber Equilibrium for each locus for each pop \n", file = "HW_raw.txt", append=F)
  for(i in 1:length(HW)){
    cat(paste(pops[i], "\n", sep=""), file = "HW_raw.txt", append=T)
    suppressWarnings(write.table(HW[[i]], "HW_raw.txt", sep="\t", fileEncoding = "UTF-8", row.names = F, append = T))
  }
  message("Outputed HW")
}

message("Calculating nulls")
NuAl<-list()
for(i in 1:npops){
  NuTemp<-PopGenReport::null.all(genssep[[i]])
  NuTemp<-NuTemp$null.allele.freq$summary2
  for(k in 1:nloc){
    qualitytest[[k]][i, "Null_allelle"]<-NuTemp["Median frequency",k]
  }
  NuAl[[i]]<-NuTemp
  message(paste0("Calculated nulls for pop ", i, " out of ", npops,  sep=""))
}
names(NuAl)<-pops
message("Calculated nulls for all pops")

if(OutputNull==T){
  cat("Null allele frequency for each locus for each pop \n", file = "Null_raw.txt", append=F)
  for(i in 1:length(NuAl)){
    cat(paste(pops[i], "\n", sep=""), file = "Null_raw.txt", append=T)
    suppressWarnings(write.table(NuAl[[i]], "Null_raw.txt", sep="\t", fileEncoding = "UTF-8", row.names = T, append = T))
  }
  message("Outputed nulls")
}

message("Calculating significance of Fis deviation from zero")

Fis<-list()
for(i in 1:npops){
  tempstat<-hierfstat::basic.stats(genssep[[i]], diploid = TRUE, digits = 2)
  tempstat<-cbind(tempstat$perloc, tempstat$n.ind.samp)
  tempstat<-tempstat[,c("1", "Ho", "Hs", "Fis")]
  colnames(tempstat)<-c("n", "Ho", "Hs", "Fis")
  tempstat<-tempstat%>%dplyr::mutate(Xi=n*Fis^2)
  tempstat<-tempstat%>%dplyr::mutate(p=dplyr::case_when(Xi>3.84&Xi<=6.63 ~ "*", Xi>6.63&Xi<=10.83 ~ "**", Xi>10.83 ~ "***", .default = "Non Significant"))
  tempstat<-tempstat%>%dplyr::mutate(p=dplyr::case_when(p!="Non Significant"&Fis<0 ~ "Significantly negative!", .default = p))
  for(k in 1:nloc){
    qualitytest[[k]][i, "Fis_significance"]<-tempstat[k, "p"]
  }
  Fis[[i]]<-tempstat
  message(paste0("Calculated Fis", i, " out of ", npops,  sep=""))
}
names(Fis)<-pops
message("Calculated significance of Fis deviation from zero")

if(OutputFis==T){
  cat("Fis each locus for each pop \n", file = "Fis_raw.txt", append=F)
  for(i in 1:length(NuAl)){
    cat(paste(pops[i], "\n", sep=""), file = "Fis_raw.txt", append=T)
    suppressWarnings(write.table(Fis[[i]], "Fis_raw.txt", sep="\t", fileEncoding = "UTF-8", row.names = T, append = T))
  }
  message("Outputed Fis")
}

message("Summarising data for each locus")

for(i in 1:nloc){
  qualitytest[[i]]%<>%dplyr::rowwise()%<>%dplyr::mutate(overall=dplyr::case_when(Hw_deviation<0.05&Null_allelle>NullTreshold&Fis_significance!="Non Significant"&Fis_significance!="Significantly negative!"~"*", .default = ""))
  tempdf<-qualitytest[[i]]

  message(paste0("Locus ", loci[i], " deviates by all three parameters in ", nrow(tempdf[tempdf$overall=="*",]), " out of ", npops))
}

if(OutputPerLoc==T){
  cat("Three quality parameters for each locus for each pop \n", file = "PerLoc.txt", append=F)
  for(i in 1:length(qualitytest)){
    cat(paste(loci[i], "\n", sep=""), file = "PerLoc.txt", append=T)
    suppressWarnings(write.table(qualitytest[[i]], "PerLoc.txt", sep="\t", fileEncoding = "UTF-8", row.names = F, append = T))
  }
  message("Outputed quality parameters")
}

quicksummary<-as.data.frame(matrix(nrow=nloc, ncol=3))
colnames(quicksummary)<-c("Locus", "N_significant_pops", "Significant_pops")
quicksummary$Locus<-loci

for(i in 1:nloc){
  tempdf<-qualitytest[[i]]
  quicksummary[i, 2]<-nrow(tempdf[tempdf$overall=="*",])
  if(quicksummary[i, 2]>0){
  spops<-tempdf[tempdf$overall=="*","Pop"][[1]]
  quicksummary[i, 3]<-paste(spops, collapse = ", ")
  }else{quicksummary[i, 3]<-""}
}
write.table(quicksummary, "Quicksummary.txt", row.names = F, fileEncoding = "UTF-8")

Res<-list(quicksummary, qualitytest, HW, NuAl, Fis)
names(Res)<-c("quicksummary", "qualitytest", "HW", "NuAl", "Fis")

return(Res)
}


#' Calculate the Garza-Williamson index (M-ratio) for microsatellite loci
#'
#'Function calculates Mratio for microsatellite loci suitable for the amplification fragment length analysis. The allelic range is calculated by automatically determining the motif on the basis of smallest difference between alleles.
#'
#' @param genal
#' Genpop object or a genslex file to be read with poppr::read.genalex()
#' @param file
#' Set to TRUE if supllying genotypes in a file or to FALSE if in an already loaded genind object. Defaults to TRUE
#' @param output
#' Set to TRUE to output results into .txt file
#' @param popnames
#' Names for the populations to use. If not supplied sequential numbers will be used.
#'
#' @returns
#' Named list with per locus and mean MRatio for for each pop
#'
#' @export
#'
#' @examples
#' TBD
MRatio<-function(genal, file=T, output=T, popnames=NULL){
  if(file==T){
    gens<-poppr::read.genalex(genal, genclone = F)
  }else{gens<-genal}
  if(!is.null(popnames)){
    pops<-popnames
  }else{
    pops<-levels(gens@pop)
  }
  locs<-levels(gens@loc.fac)
  genssep<-seppop(gens, drop=T)
  meansbypop<-c()
  Results<-list()
  for(i in 1:length(pops)){
    curgen<-genssep[[i]]
    GW<-c()
    for(l in 1:length(locs)){
      curloc<-locs[l]
      loc<-curgen@all.names[[curloc]]
      loc<-as.numeric(loc)
      loc<-sort(loc)
      k<-curgen@loc.n.all[[curloc]]
      r<-diff(range(loc))/min(abs(diff(loc)))
      M<-k/(r+1)
      GW[l]<-M
    }
    names(GW)<-locs
    GW["Mean"]<-mean(GW)
    meansbypop[i]<-GW["Mean"]
    Results[[i]]<-GW
  }
  names(Results)<-pops
  names(meansbypop)<-pops
  Results[["MeansByPop"]]<-meansbypop
  if(output==T){
    cat("Garza-Williamson index for each pop \n --- \n",  file = "GW.txt", append = F)
    cat("Means for each pop \n",  file = "GW.txt", append = T)
    for(i in 1:length(pops)){
      cat(c(pops[i], ": ", Results[["MeansByPop"]][i], "\n"),  file = "GW.txt", append = T, sep="")
    }
    cat("--- \n Values for each pop for each locus \n",  file = "GW.txt", append = T)
    for(i in 1:length(pops)){
      cat(c(pops[i], "\n"),  file = "GW.txt", append = T)
      for(l in 1:length(locs)){
        cat(c(locs[l], ": ", Results[[i]][l], "\n"),  file = "GW.txt", append = T, sep="")
      }
    }
  }
  return(Results)
}


