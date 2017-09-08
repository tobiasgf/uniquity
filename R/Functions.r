#' Estimating biotic uniqueness of biological samples
#' @description This algorithm \code{uniquity} estimates the uniqueness of a set of communities
#' @param sstable a data.frame with with a species/site matrix (or an OTU table)
#' that has sites/samples as columns and species/OTUs as rows, and observations as
#' presence/absence data or abundance data.
#' @param classes a list (dataframe) assigning samples to a set of unique classes (e.g. habitat-types). Each class
#'  represents an exclusive proportion af the investigated total area). (data.frame should contain these column names: site, class).
#'  If no classes are supplied, each samples is assigned to a unique class, and these will receive equal weights.
#' @param weights a list (dataframe) with the weight (proportion of total) of each
#' class. (data.frame should contain these column names: class, weight).
#' If no weights are supplied, classes are assigned equal weights.
#' @param presabs TRUE/FALSE. reduce species/site matrix to presence absence (0/1). Default is TRUE
#' @param background_sites any positive number. The number of background sites to construct
#' for counting the number of unique species/OTUs of the sample/site in focus (i.e the uniquity of the site/sample).
#' each each compare. If \code{background_sites} < 1, the algorithm uses that proportion
#'   of the sites investigated as the number of background sites. If \code{background_sites} =>1, that absolute number
#'   of background sites wll be constructed. Default is 0.66 (=66\% of samples/sites)
#' @param rep any positive integer. The number of replicates pr site evaluation.
#' (i.e. the number of times to construct a pool of background sites for evaluating each site/sample). Default is 100.
#' @param non_correspondence TRUE/FALSE. If set to TRUE non-matching entries (classes, weights, sites) will be removed from the dataet.
#' If set to FALSE, non non-matching entries will throw an error. Can be usefull for subsetting data, 
#' as subsetting will only be necessary in one part of the data (species-site table)
#' @param total_coverage any number from 0 to 1. The fraction of total area investigated. If set to 0 the supplied
#' weights will be used, unless the total weight exceeds 1, in wich case weights are scaled for a total of 1. If set to any
#' number above 0 (including 1), weights are scaled for a total weight of \code{total_coverage}.
#' @return Function \code{uniquity} returns a list of results based on the input
#' table and parameters.
#'   \enumerate{
#'   \item \code{simple_site_score} the sum of inverse species weights, as a fast approximation of the estimated uniquity.
#'   \item \code{uniquity_score} the simulated uniquity contribution of each site, relative to species richness.
#'   \item \code{adjusted_uniquity} - number of curated (parent) OTUs
#'   \item \code{uniquity_table} the final table with number of unique observations pr species per site
#'   \item \code{uniquity_species_scores} sum of contribution per species.
#'   \item \code{Unique_species_avg} average of contribution per species for sites where it occurs.
#'   \item \code{Summed_species_weights} sum of weight pr species (used for drawing species for background sites).
#'   \item \code{Average_species_weight_pr_site} average species weight pr site (with a weight).
#'   \item \code{Top20_unique_sites} top 20 unique sites.
#'   \item \code{Top20_unique_species} top 20 species with highest uniquity impact.
#'   \item \code{Accumulated_uniqueness} uniquity uniquity metric per replicate.
#'   \item \code{site_richness} species richness per site.
#'   \item \code{replicates} The number of replicates used for each site uniquity estimation.
#'   \item \code{background_sites} The number of background sites sampled.
#'   \item \code{Removed_classes} classes removed due to non correspondence.
#'   \item \code{Removed_sites} sites removed due to non correspondence.
#'   \item \code{total_coverage} sum of site weights used in weighting of species.}
#' @examples
#' uniquity(my_table)
#' uniquity(my_table, my_classes, my_weights, presabs = FALSE, rep = 50, size = 100, non_correspondence = TRUE, total_coverage = 0.5)
#' @author Tobias Guldberg Frøslev
#' @export

#Variables
# ft: final table
# sn: site number
# lr: loop rounds
# rep: repetitions. Number of resampling events pr plot evaluation
# css: current site species. IDs of the species present in the currently evaluated plot
# rfp: number of background plots used to evaluate each site
# crsp : current background site pool. Pooled species IDs of full background site pool.
# csu : current site unique. Unique species IDs of current site in current repetition
# csut : current site unique total. vector of all unique species IDs pr current site.
# nsi : number of sites in study
# nsp : number of species in study
# spr : species richness pr site
# siw : site weights. percentage of area in the study represented by each site.
# ssw : species-site weights. calculated weight pr species pr site (site). siw applied to relevant species in each site
# spw : species weights. summed weight of each species across all sites
# claw : class weight. he figure needed to multiply each weight with. (equal to 1/(number of sites referred to class))
# classw : temporary datafra to calculate adjusted site weights
# rep: number of replicates pr site evaluation
# sv: vector of species IDs (numbers)
# refspecies: input number of species to draw for each background site. (0: use random)
# ssrsp: the actual number of species being drawn for each background site.


uniquity2 <- function(sstable, classes = NULL, weights = NULL, presabs = TRUE,  background_sites = 50, rep = 50, non_correspondence = TRUE, total_coverage = 0){
 #Exclude singletons if selected
 
 require(vegan)
 
 #Check if a species/site table is present 
 if (missing(sstable)) stop("Need to specify a valid species/site matrix as a data.frame with sites as rownames and species/OTUs as column names")
 
 sstable_names <- row.names(sstable)
 sstable <- as.matrix(sstable)
 
 removed_classes <- NULL
 removed_sites <- NULL
 
 #Check if valid combination of classes and weights are present.
 if (is.null(classes)) {
  if (is.null(weights)) {
   message("no classes supplied. Assigning separate class to each site!")
   classes <- data.frame(site=as.character(row.names(sstable)),class=rep(1,length(sstable_names)))
  } else {
   stop("Need to specify a valid classes, if supplying weights")
  }
 }
 
 #Apply equal weights when no weights supplied
 if (is.null(weights)) {
  message("no weights supplied. Assigning equal weights to all classes!")
  weights <- data.frame(class=names(table(classes$class)),weight=1)
 }
 
 sorted_sstable_names <- as.character(sort(sstable_names))
 sorted_class_site_names <- as.character(sort(classes$site))
 
 sstable_uniq <- paste(setdiff(sorted_sstable_names,sorted_class_site_names),collapse=" ")
 class_uniq <- paste(setdiff(sorted_class_site_names,sorted_sstable_names),collapse=" ")
 
 if(!identical(sorted_sstable_names, sorted_class_site_names)){
  if(non_correspondence){
   common_site_names <- intersect(sorted_sstable_names,sorted_class_site_names)
   rem_s_uniq <- paste(setdiff(sorted_sstable_names,common_site_names),collapse=" ")
   rem_c_uniq <- paste(setdiff(sorted_class_site_names,common_site_names),collapse=" ")
   message("Non-identical classes between classes and weights")
   message(paste("removing the following from the species_site_table:",rem_s_uniq))
   message(paste("removing the following from the class file:",rem_c_uniq))
   removed_sites <- paste(rem_s_uniq,rem_c_uniq,collapse=" ")
   sstable <- sstable[row.names(sstable) %in% common_site_names,]
   classes <- classes[classes$site %in% common_site_names,]
  } else {
   message("Site/sample names are not identical between site/species matrix and class-data")
   message(paste("unique to site/species table:",sstable_uniq))
   message(paste("unique to classes:",class_uniq))
   stop("stopped!")
  }
 }
 
 sorted_class_site_classes <- sort(names(table(classes$class)))
 sorted_weight_classes <- sort(weights$class)
 
 class_uniq <- paste(setdiff(sorted_class_site_classes,sorted_weight_classes),collapse=" ")
 weight_uniq <- paste(setdiff(sorted_weight_classes,sorted_class_site_classes),collapse=" ")
 
 if(!identical(sorted_class_site_classes, sorted_weight_classes)){
  if(non_correspondence){
   common_class_names <- intersect(sorted_class_site_classes, sorted_weight_classes)
   rem_c_uniq <- paste(setdiff(sorted_class_site_classes,common_class_names),collapse=" ")
   rem_w_uniq <- paste(setdiff(sorted_weight_classes,common_class_names),collapse=" ")
   message("Non-identical classes between classes and weights")
   message(paste("removing the following from the class file:",rem_c_uniq))
   message(paste("removing the following from the weight file:",rem_w_uniq))
   classes <- classes[classes$class %in% common_class_names,]
   weights <- weights[weights$class %in% common_class_names,]
   removed_classes <- paste(rem_w_uniq,rem_c_uniq,collapse=" ")
  } else {
   message("Names of classes  are not identical between class data and weight data")
   message(paste("class names unique to classes:",class_uniq))
   message(paste("class names unique to weights:",weight_uniq))
   stop("stopped!")
  }
 }
 
 sorted_sstable_names <- as.character(sort(sstable_names))
 sorted_class_site_names <- as.character(sort(classes$site))
 
 sstable_uniq <- paste(setdiff(sorted_sstable_names,sorted_class_site_names),collapse=" ")
 class_uniq <- paste(setdiff(sorted_class_site_names,sorted_sstable_names),collapse=" ")
 
 if(!identical(sorted_sstable_names, sorted_class_site_names)){
  if(non_correspondence){
   common_site_names <- intersect(sorted_sstable_names,sorted_class_site_names)
   rem_s_uniq <- paste(setdiff(sorted_sstable_names,common_site_names),collapse=" ")
   rem_c_uniq <- paste(setdiff(sorted_class_site_names,common_site_names),collapse=" ")
   message("Non-identical classes between classes and weights")
   message(paste("removing the following from the species_site_table:",rem_s_uniq))
   message(paste("removing the following from the class file:",rem_c_uniq))
   removed_sites <- paste(removed_sites,rem_s_uniq,rem_c_uniq,collapse=" ")
   sstable <- sstable[row.names(sstable) %in% common_site_names,]
   classes <- classes[classes$site %in% common_site_names,]
  }
 }
 
 sstable_names <- row.names(sstable)
 
 tab <- sstable
 #distribute class weights equally among sites assigned to each class
 claw <- as.data.frame(1/table(classes$class)) #  proportion of weight to assign to each site within each class
 classw <- merge(weights,claw, by.x = "class", by.y="Var1")
 classw$adj <- classw$weight * classw$Freq #calculate adjusted weight pr site within each class
 siw <- merge(classes,classw, by="class")[,c("site","adj")]
 row.names(siw) <- siw$site
 siw <- siw[sstable_names,]$adj #sort site-weights vector
 siw_sum <- sum(siw)
 if (siw_sum > 1 & total_coverage == 0) {
  message("total site-weights are above 1, scaling weight for a total of 1")
  siw <- siw/siw_sum
 }
 if (siw_sum <= 1 & total_coverage == 0) {
  message(paste0("Using weights as supplied, for total site weights of ",siw_sum))
 }
 if (total_coverage == 1){
  siw <- siw/siw_sum
  message(paste0("Scaling site weight for a total of :",sum(siw)))
 }
 if (siw_sum < 1 & total_coverage > 0 & total_coverage < 1) {
  siw <- (siw/siw_sum)*total_coverage
  message(paste0("Scaling site weight for a total of ", sum(siw)))
 }
 
 if (presabs){ tab[tab>0] <- 1 }# Convert to pres/abs table
 nsi <- nrow(tab) # number of sites
 nsp <- ncol(tab) # number of species/otus
 spr <- rowSums(tab>0) # Get number of species
 
 #Setting resampling parameters

 
 lr <- nsi*rep*size # calculate loop rounds (to be able to make "progress bar")
 sv <- seq(1:nsp) # make a vector of species-numbers
 
 #Normalized_table1C <- decostand(tab,"total") # Normalize Table1 to get relative share pr species (this choice has been deactivated)
 ssw <- apply(tab,2,function(x)x*siw) # multiply relative share by siteprobs to get species_weight relative to plot
 #ssw <- apply(Normalized_table1C,2,function(x)x*siw) # multiply relative share by siteprobs to get species_weight relative to plot
 spw <- colSums(ssw) ## Add all species probs pr site from all sites to get total weight
 spw2 <- spw
 names(spw2) <- 1:length(spw2)
 sspl <- length(spw2)

 norm_score <- max(spw) - spw
 score_tab <- apply(tab,1,function(x)x*norm_score)
 site_cores <- colSums(score_tab)
 
 csut <- vector(mode="list", length=nsi) # prepare list for total unique species for each site
 ft <- tab # Prepare a table for results
 ft[ft>0] = 0  # resetting the results table
 total_table <- list()
 site_table <- list()
 site_rep_uniq <- ft#[,0]
 
 for (sn in 1:nsi){      # looping through the sites one by one
  print(paste0("evaluating sample ", sn ," of ", nsi)) # make a progressline
  for (rn in 1:rep){   # looping through the replicates one by one
   css <- which(tab[sn,]>0) # which species (ID) are present in the current site?
   crsp = NULL #resetting the current background set
   for (jj in 1:size){
    crsp <- unique(c(crsp,which(runif(sspl, 0, 100) < spw2)))
   }
   csu <- css[!(css %in% crsp)] # find species IDs unique to currently investigated site compared to total background pool
   site_table[[rn]] <- csu
   site_rep_uniq[sn,rn] <- length(csu)
   csut[[sn]] <- c(csut[[sn]],csu) # add unique IDs of current replicate to growing vector of unique IDs for current plot ()
  }
  total_table[[sn]] <- site_table
 }
 
 # process list of vectors of unique species 
 for (sn in 1:nsi){
  ft[sn,as.numeric(names(table(csut[[sn]])))] <- table(csut[[sn]])
 }
 Uniqueness <- rowSums(ft)/rep #Calculating average uniqueness pr plot
 Adj_uniquity <- Uniqueness/spr
 Uniquenss_table <- ft
 Unique_species <- colSums(ft)
 av_un_sp <- (Unique_species/colSums(ft != 0))/rep
 nzmean <- function(x) {
  zvals <- x==0
  if (all(zvals)) 0 else mean(x[!zvals])
 }
 average_ssw <- apply(ssw,1,nzmean)
 top20 <- sort(Uniqueness,decreasing = TRUE)[1:20]
 top20sp <- sort(av_un_sp,decreasing = TRUE)[1:20]
 result <- list(simple_site_score=site_cores,
                uniquity_score=Uniqueness, 
                adjusted_uniquity=Adj_uniquity, 
                uniquity_table=Uniquenss_table, 
                uniquity_species_scores=Unique_species, 
                Unique_species_avg=av_un_sp, 
                Summed_species_weights=spw, 
                Average_species_weight_pr_site=average_ssw, 
                Top20_unique_sites=top20, 
                Top20_unique_species=top20sp, 
                Accumulated_uniqueness=total_table, 
                Site_richness=spr,
                replicates=rep,
                background_sites=background_sites,
                Removed_classes = removed_classes,
                Removed_sites = removed_sites,
                ft, total_coverage=sum(siw)
 )
 return(result)
}
