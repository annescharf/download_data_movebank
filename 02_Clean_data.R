#### in this script - data cleaning -
## remove empty locations 
## remove duplicated timestamps
## remve outliers based on speed -- could not be done here - pc crashed due to RAM -- will filter after reducing to one loc per hour

library(move2)
library(units)
# library(doParallel)
# library(plyr)
# mycores <- detectCores()-1
# registerDoParallel(mycores)
library(dplyr)

pathTOfolder <- "path_to_folder_where_all_data_will_be_saved/"
pthDownld <- paste0(pathTOfolder,"1.Indiv_data/")
dir.create(paste0(pathTOfolder,"2.Indiv_data_clean_empty_duply"))
pthClean <- paste0(pathTOfolder,"2.Indiv_data_clean_empty_duply/")

flsMV <- list.files(pthDownld, full.names = F)

## remove empty locs
## remove duplicated ts
start_time <- Sys.time()
lapply(flsMV, function(indPth){
  # llply(flsMV, function(indPth){
  mv2 <- readRDS(paste0(pthDownld,indPth))
  #### remove empty locs
  mv2 <- mv2[!sf::st_is_empty(mv2),]
  ## sometimes only lat or long are NA
  crds <- sf::st_coordinates(mv2)
  rem <- unique(c(which(is.na(crds[,1])),which(is.na(crds[,2]))))
  if(length(rem)>0){mv2 <- mv2[-rem,]}
  ## remove 0,0 coordinates
  rem0 <- which(crds[,1]==0 & crds[,2]==0)
  if(length(rem0)>0){mv2 <- mv2[-rem0,]}
  
  #### retain the duplicate timestamp entry which contains the least number of columns with NA values
  mv2 <- mv2 %>%
    mutate(n_na = rowSums(is.na(pick(everything())))) %>%
    arrange(n_na) %>%
    mt_filter_unique(criterion='first') %>% # this always needs to be "first" because the duplicates get ordered according to the number of columns with NA. 
    dplyr::arrange(mt_track_id()) %>%
    dplyr::arrange(mt_track_id(),mt_time())
   saveRDS(mv2, file=paste0(pthClean,indPth))
} )
# } ,.parallel = T)
end_time <- Sys.time()
end_time - start_time

## excluding outliers by speed should be done here, but the large individulas crash the pc due to RAM-> doing filtering on 1h tracks
