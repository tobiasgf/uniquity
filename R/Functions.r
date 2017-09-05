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
#' @param size any positive number. The number of reference sites to construct
#' for counting the number of unique species/OTUs of the sample/site in focus (i.e the uniquity of the site/sample).
#' each each compare. If \code{size} < 1, the algorithm uses that proportion
#'   of the sites investigated as the number of reference sites. If \code{size} =>1, that absolute number
#'   of reference sites wll be constructed. Default is 0.66 (=66\% of samples/sites)
#' @param spec_acc_size any number between 0 and 1. The number of reference sites to construct
#' will be based a species accumulation curve, so that the number of sites where the fraction \code{spec_acc_size}
#'  of species are sampled.
#' This option will be used when \code{spec_acc_size} >0, otherwise the settings of \code{size} will be used.
#' @param rep any positive integer. The number of replicates pr site evaluation.
#' (i.e. the number of times to construct a pool of reference sites for evaluating each site/sample). Default is 100.
#' @param presabs TRUE/FALSE. reduce species/site matrix to presence absence (0/1). Default is TRUE
#' @param refspecies any positive integer or 0. The number of species/OTUs to draw to construct a 
#' reference site. If \code{refspecies} = 0,  the function to draw a random figure based on
#' the species richness and weight of investigated sites. Any other figure will 
#' result in a fixed number being drawn for each reference site.
#' @param abundance_cutoff any integer. Remove species that occur in \code{size} sites or less
#' before simulation. Default = FALSE
#' @param replace_sp TRUE/FALSE. sample species with replacement.
#' Default = FALSE
#' @return Function \code{uniquity} returns a list of results based on the input
#' table and parameters.
#'   \enumerate{
#'   \item \code{uniquity_score} the simulated uniquity contribution of each site, relative to species richness.
#'   \item \code{adjusted_uniquity} - number of curated (parent) OTUs
#'   \item \code{uniquity_table} the final table with number of unique observations pr species per site
#'   \item \code{uniquity_species_scores} sum of contribution per species.
#'   \item \code{Summed_species_weights} sum of weight pr species (used for drawing species for reference sites).
#'   \item \code{Average_species_weight_pr_site} average species weight pr site (with a weight).
#'   \item \code{Top20_unique_sites} top 20 unique sites.
#'   \item \code{Top20_unique_species} top 20 species with highest uniquity impact.
#'   \item \code{removed_singleton_species} which species were removed (if rem-single is et to TRUE).
#'   \item \code{Accumulated_uniqueness} uniquity uniquity metric per replicate.
#'   \item \code{site_richness} species richness per site.}
#' @examples
#' uniquity(my_table)
#' uniquity(my_table, my_classes, my_weights, presabs = FALSE, rep = 50, size = 0, spec_acc_size = 0.75, refspecies = 50, abundance_cutoff = 1, replace_sp = FALSE, non_correspondence = FALSE)
#' @author Tobias Guldberg Frøslev
#' @export

#Valriables
# ft: final table
# sn: site number
# lr: loop rounds
# rep: repetitions. Number of resampling events pr plot evaluation
# css: current site species. IDs of the species present in the currently evaluated plot
# rfp: number of reference plots used to evaluate each site
# crsp : current reference site pool. Pooled species IDs of full reference site pool.
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
# refspecies: input number of species to draw for each reference site. (0: use random)
# ssrsp: the actual number of species being drawn for each reference site.

uniquity <- function(sstable, classes = NULL, weights = NULL, presabs = TRUE, rep = 100, size = 0, spec_acc_size = 0.75, refspecies = 0, abundance_cutoff = 0, replace_sp = FALSE, non_correspondence = FALSE){
  #Exclude singletons if selected

require(vegan)
 
#Check if a species/site table is present 
if (missing(sstable)) stop("Need to specify a valid species/site matrix as a data.frame with sites as rownames and species/OTUs as column names")

sstable_names <- row.names(sstable)

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

removed_species = NULL
  if (abundance_cutoff > 1){
   removed_species <- names(sstable[,colSums(sstable) < (abundance_cutoff+1)])
   sstable <- sstable[,colSums(sstable) > abundance_cutoff]
  }

tab <- sstable
  #distribute class weights equally among sites assigned to each class
  claw <- as.data.frame(1/table(classes$class)) #  proportion of weight to assign to each site within each class
  classw <- merge(weights,claw, by.x = "class", by.y="Var1")
  classw$adj <- classw$weight * classw$Freq #calculate adjusted weight pr site within each class
  siw <- merge(classes,classw, by="class")[,c("site","adj")]
  row.names(siw) <- siw$site
  siw <- siw[sstable_names,]$adj #sort site-weights vector

  if (presabs){ tab[tab>0] <- 1 }# Convert to pres/abs table
  nsi <- nrow(tab) # number of sites
  nsp <- ncol(tab) # number of species/otus
  spr <- rowSums(tab>0) # Get number of species

  #Setting resampling parameters
  if (size<1) {rfp <- size*nsi} else {rfp <- size} # are we using absolute sample size or relative (how many reference sites to construct)
  if (spec_acc_size>0) {
   specac <- specaccum(sstable)
   rfp <- which(abs(specac$richness-max(specac$richness)*spec_acc_size) == min(abs((specac$richness-max(specac$richness)*spec_acc_size)))) # set nsites to the number of sites where 75% of the species have been sampled. 
  }
  
  lr <- nsi*rep*rfp # calculate loop rounds (to be able to make "progress bar")
  sv <- seq(1:nsp) # make a vector of species-numbers
  
  #Normalized_table1C <- decostand(tab,"total") # Normalize Table1 to get relative share pr species (this choice has been deactivated)
  ssw <- apply(tab,2,function(x)x*siw) # multiply relative share by siteprobs to get species_weight relative to plot
  #ssw <- apply(Normalized_table1C,2,function(x)x*siw) # multiply relative share by siteprobs to get species_weight relative to plot
  spw <- colSums(ssw) ## Add all species probs pr site from all sites to get total weight
  if (refspecies > 0){ssrsp <- refspecies} # Are we using an absolute number of reference species pr site, then set that number.
  
  norm_score <- max(spw) - spw
  score_tab <- apply(tab,1,function(x)x*norm_score)
  site_cores <- colSums(score_tab)
  
  csut <- vector(mode="list", length=nsi) # prepare list for total unique species for each site
  ii=1 # counter used for progress status
  ft <- tab # Prepare a table for results
  ft[ft>0] = 0  # resetting the results table
  total_table <- list()
  site_table <- list()
  site_rep_uniq <- ft#[,0]
  for (sn in 1:nsi){      # looping through the sites one by one
    print(paste0("progress: ", round(((ii/lr) * 100),0) ,"%")) # make a progressline
    for (rn in 1:rep){   # looping through the replicates one by one
      css <- which(tab[sn,]>0) # which species (ID) are present in the current site?
      crsp = NULL #resetting the current reference set
      for (jj in 1:rfp){  # loop through and construct the number of reference sites
        if (refspecies == 0){ssrsp <- sample(spr, 1,replace = FALSE, prob = siw)} # set number of reference species to draw, if that option is selected
        crs <- sample(sv, ssrsp ,replace = replace_sp, prob = spw) # draw a number of species IDs from pool. Number is chosen randomly among real sites according to their weights
        crsp <- c(crsp,crs) # add current selection to total vector of species IDs drawn.
        ii=ii+1
      }
      csu <- css[!(css %in% crsp)] # find species IDs unique to currently investigated site compared to total reference pool
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
  top10 <- sort(Uniqueness,decreasing = TRUE)[1:20]
  top10sp <- sort(Unique_species,decreasing = TRUE)[1:20]
  result <- list(simple_site_score=site_cores,
                 uniquity_score=Uniqueness, 
                 adjusted_uniquity=Adj_uniquity, 
                 uniquity_table=Uniquenss_table, 
                 uniquity_species_scores=Unique_species, 
                 Unique_species_avg=av_un_sp, 
                 Summed_species_weights=spw, 
                 Average_species_weight_pr_site=average_ssw, 
                 Top20_unique_sites=top10, 
                 Top20_unique_species=top10sp, 
                 removed_singleton_species=removed_species, 
                 Accumulated_uniqueness=total_table, 
                 Site_richness=spr,
                 replicates=rep,
                 reference_sites=rfp,
                 Removed_classes = removed_classes,
                 Removed_sites = removed_sites,
                 ft
                 )
  return(result)
}