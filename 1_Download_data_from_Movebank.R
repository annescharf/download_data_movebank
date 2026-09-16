## Title: Download data of multiple studies from Movebank
## Author: Anne K Scharf, MPI of Animal Behavior
## Date: September 2026
## Description: this script is developed to download hundreds of studies, but of 
##              course it also works just for one or a few studies
##  This script has 3 steps:
##    1. creates table of studies to be downloaded. In this example studies 
##        from a given movebank user and all studies with license "CC_0","CC_BY" 
##        & "CC_BY_NC" are downloaded. This table can also be filtered by taxon 
##        or any other column of interest. A "metadata table" is created containing 
##        the study name, contact person and license type and terms among others.
##    2. data are downloaded by individual. The reason for not downloading by 
##        study, is that some studies contain an enormous amount of individuals, 
##        or individuals with a huge amount of data, and managing these huge 
##        datasets in one table crashes R, in the best case makes it very slow. 
##        Working on an individual base makes it easier to manage. 
##        License terms are accepted in the code if need be.
##    3. basic cleaning of the data: ensuring that timestamps are ordered, 
##        empty locations removed, duplicated timestemps removed.
##
## ATTENTION: please go though the license terms and types of all the 
##    studies that you are downloading and act accordingly. If you have any 
##    doubt if and how you can use the data for your own work, please 
##    get in touch with the contact person of the study.

## OUTPUT:
## Step 1: - `full_table_all_studies.rds`: table containing all studies to download (one row per study)
##         - `metadata_all_studies.rds`: table containing the metadata, including license and contact per study
## Step 2: - one folder containing one .rds file per individual downloaded
## Step 3: - one folder containing one .rds file per individual downloaded. File names are the same as in step 2, but they are placed              in a new folder with an intuitive name to identify the content


#######-----------------------------########
## 1. create table of studies to download ##
#######-----------------------------########
### gathers all studies to download, 
### those shared with the specific movebank user 
### and/or those that are publicly available
### metadata table is created and saved including study name, 
### owner, license terms, download date, etc

library(move2)
library(units)
library(dplyr)

# specify account to use in the R session
# Follow this link to find out how to setup the credentials via keyring: https://bartk.gitlab.io/move2/articles/movebank.html
# keyring::key_list()
options("move2_movebank_key_name" = "movebank")

dir.create("MBdata")
pathTOfolder <- "./MBdata/"

#### downloading studies to which the user "XXXX" has been added as collaborator or manager
# download list of studies available through this account
all_shared <- movebank_download_study_info(study_permission=c("data_manager","collaborator"))
## here you can also filter the table as below

### searching for public studies
all <- movebank_download_study_info() # some studies have years in weird formats, just ignore this warning message
all <- all[grep("GPS", all$sensor_type_ids),] # studies can have multiple sensors, making sure gps is included, adjust as needed
all <- all[all$number_of_deployed_locations > units::set_units(0,"count"),] # removing those with 0 locations
all <- all[!is.na(all$number_of_deployed_locations),] # removing those with no deployed locations
all <- all[!is.na(all$taxon_ids),] ## removing those with NO taxon
# all <- all[grep("Ciconia ciconia", all$taxon_ids),] ## selecting only studies that contain species of interest
all <- all[!all$is_test==T,] ## removing studies marked as tests
all_open <- all[which(all$license_type %in% c("CC_0","CC_BY","CC_BY_NC")),] 
## - CC_O: can use the data, do not need to mention names
## - CC_BY: can use the data, but names of owners should appear somewhere, eg acknowledgments
## - CC_BY_NC: can use the data, but names of owners should appear somewhere, eg acknowledgments

### making one large table and removing duplicated studies
allstudies <- rbind(all_shared,all_open)
allstudies <- allstudies[!duplicated(allstudies$id),] ## when duplicated, entry from all shared will be kept
allstudies$download_date <- Sys.Date()
saveRDS(allstudies, paste0(pathTOfolder,"full_table_all_studies.rds")) ## saving all columns just in case they need to be revisited

## creating table to record which studies and license terms have been downloaded
metadata_studies <- allstudies[,c(
  "id",
  "name",
  "taxon_ids",
  "number_of_individuals",
  "timestamp_first_deployed_location",
  "timestamp_last_deployed_location",
  "number_of_deployed_locations",
  "principal_investigator_name",
  "principal_investigator_email",
  "contact_person_name",
  "citation",
  "license_terms",
  "license_type",
  "download_date"
)]
saveRDS(metadata_studies, paste0(pathTOfolder,"metadata_all_studies.rds"))
## ensure to inspect this table before using the data to comply with all license terms

######-----------------------#######
## 2. download data by individual ##
######-----------------------#######
### data is downloaded from movebank, each individual is downloaded separately. 
### license agreements are accepted.
### in case the download gets interrupted because of internet issues or 
### connection to movebank, the script can be run again and it will check 
### which data has already been downloded and only download the missing individuals.


library(move2)
library(bit64)
library(units)
library(R.utils)

# specify account to use in the R session
# keyring::key_list()
options("move2_movebank_key_name" = "movebank")

pathTOfolder <- "./MBdata/"
dir.create(paste0(pathTOfolder,"01_MB_indv_mv2"))
pthDownld <- paste0(pathTOfolder,"01_MB_indv_mv2/")

### studies to download
allstudies <- readRDS(paste0(pathTOfolder,"full_table_all_studies.rds"))
Ids_toDo <- allstudies$id 

#### download by individual. object "result" only contains the error messages ######
start_time <- Sys.time()
# studyId <- Ids_toDo[1]
results <- lapply(Ids_toDo, function(studyId)try({
  ## create table with individuals per study, to be able to download per individual
  class(studyId) <- "integer64" ## lapply changes the class when looping though it
  print(studyId)
  
  ## if license terms have to be accepted, this is done here. These license terms are recorded in the metadata_table
  reftb <-  tryCatch({
    movebank_download_deployment(study_id=studyId, omit_derived_data=F)
  }, error = function(e) {
    movebank_download_deployment(study_id=studyId, omit_derived_data=F,
                                 'license-md5'= sub('...Alternat.*','',sub('.*se-md5.=.','',as.character(rlang::catch_cnd(movebank_download_study(studyId))))))
  })
  
  reftb <- reftb[reftb$number_of_events > units::set_units(0,"count"),]
  ## when you are interested in the data of one specific species, you might want to filter again here for the species. Some studies on Movebank contain multiple studies. This way you can ensure to only download your species of interest.
  
  ## intuitively one would use "individual_local_identifier" problem is: it sometimes does not exist, names often contains symbols that mess with R like e.g. "/". The "individual_id" is an internal number not visible on the webpage, it is also consistent unless the study gets deleted and uploaded again. When downloading many individuals, this is the safest option
  reftb$pthName <- paste0(studyId,"_",reftb$individual_id,".rds")
  
  ## here it is checked which indiv have already been downloaded and which are missing
  doneIndv <- list.files(pthDownld)[grep(studyId,list.files(pthDownld))] 
  allStInd <- reftb$pthName ## individuals in study
  missInd <- allStInd[!allStInd%in%doneIndv] ## missing indiv
  print(paste0("done:",length(doneIndv),"-todo:",length(missInd)))
  todoIndv <-  reftb$individual_id[reftb$pthName%in% missInd] 
  
  ## download each individual separatly
  # ind <- todoIndv[2]
  results2 <- lapply(todoIndv, function(ind)try({
    class(ind) <- "integer64" ## lapply changes the class when looping though it
    print(paste0(studyId,"_",ind))
    mv2 <- movebank_download_study(studyId,
                                   sensor_type_id=c("gps"),                                                                         
                                   individual_id=ind, 
                                   attributes = c("individual_local_identifier","deployment_id"), ## here only lat, lon, time, and the stated columns are downloaded. 
                                   # attributes="all", if all attributes should be downloaded, use this argument instead
                                   timestamp_end=as.POSIXct(Sys.time(), tz="UTC")) # to avoid locations in the future
    saveRDS(mv2, file=paste0(pthDownld,studyId,"_",ind,".rds"))
  }))
}))
end_time <- Sys.time()
end_time-start_time # 

is.error <- function(x) inherits(x, "try-error")
table(vapply(results, is.error, logical(1)))
names(results) <- seq_along(results)
results[vapply(results, is.error, logical(1))]
# Check studies that returned errors:
giveError <- Ids_toDo[vapply(results, is.error, logical(1))]



#######-----------------#######
## 3. basic cleaning of data ##
#######-----------------#######
## data are cleaned: empty locations, "0,0" coordinates and duplicated timestamps are removed

library(move2)
library(units)
library(dplyr)

## in case in parallel is an option
# library(doParallel)
# library(plyr)
# mycores <- detectCores()-1
# registerDoParallel(mycores)
# library(dplyr)

pathTOfolder <- "./MBdata/"
pthDownld <- paste0(pathTOfolder,"01_MB_indv_mv2/")
dir.create(paste0(pathTOfolder,"02_MB_indv_mv2_clean"))
pthClean <- paste0(pathTOfolder,"02_MB_indv_mv2_clean/")

flsMV <- list.files(pthDownld, full.names = F)
done <- list.files(pthClean, full.names = F) #checking which have been already done in case an error occurs and script stops
flsMV <- flsMV[!flsMV%in%done]

## remove empty locs, 0,0 corrds and duplicated ts
start_time <- Sys.time()
# indPth <- flsMV[10]
lapply(flsMV, function(indPth){
  # llply(flsMV, function(indPth){
  mv2 <- readRDS(paste0(pthDownld,indPth))
  if(!mt_is_track_id_cleaved(mv2)){mv2 <- mv2 |> dplyr::arrange(mt_track_id(mv2))} ## order by tracks
  if(!mt_is_time_ordered(mv2)){mv2 <- mv2 |> dplyr::arrange(mt_track_id(mv2),mt_time(mv2))} # order time within tracks
  if(!mt_has_no_empty_points(mv2)){mv2 <- mv2[!sf::st_is_empty(mv2),]} ## remove empty locs
  
  ## sometimes only lat or long are NA
  crds <- sf::st_coordinates(mv2)
  rem <- unique(c(which(is.na(crds[,1])),which(is.na(crds[,2]))))
  if(length(rem)>0){mv2 <- mv2[-rem,]}
  
  ## remove 0,0 coordinates
  rem0 <- which(crds[,1]==0 & crds[,2]==0)
  if(length(rem0)>0){mv2 <- mv2[-rem0,]}
  
  ## retain the duplicate entry which contains the least number of columns with NA values
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
end_time - start_time # 

