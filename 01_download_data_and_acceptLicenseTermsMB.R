#### in this scrip -  data download -
## search for studies to download
## download data per individual/entire study
## accept license terms "automatically"


library(move2)
library(units)
library(dplyr)
library(bit64)

## for instructions on providing credentials and more details see vignette: https://bartk.gitlab.io/move2/articles/movebank.html

# I would recommend to download and save each study separately, It makes debugging easier. And also putting all studies in one large object can be tricky as some of them can be enormously large
pathTOfolder <- "path_to_folder_where_all_data_will_be_saved/"
dir.create(paste0(pathTOfolder,"1.Indiv_data"))
pthDownld <- paste0(pathTOfolder,"1.Indiv_data/")

######################################
## find list of studies to download 
#######################################

## get a table with metadata of all available studies on movebank. This table than can be filtered to get the studies that one is looking for.
all <- movebank_download_study_info() # some studies have years in weird formats, just ignore this warning message
all <- all[all$i_have_download_access==T,] # excluding studies for which one does not have download rights
all <- all[grep("Ciconia ciconia", all$taxon_ids),] ## selecting studies that contain stork data
all <- all[grep("GPS", all$sensor_type_ids),] # studies can have multiple sensors, making sure gps is included
all <- all[all$number_of_deployed_locations > units::set_units(0,"count"),] # removing those with 0 locations
all <- all[!is.na(all$number_of_deployed_locations),] # removing those with no deployed locations
all <- all[which(all$license_type == "CC_0"),] ## !there are more license types

###############################
### save a table with license type and license terms to have them all in one place
##################################
license_table <- all[,c("id","name","taxon_ids","principal_investigator_name","contact_person_name","license_type","license_terms")]
write.csv(license_table, paste0(pathTOfolder,"license_terms_table.csv"))

#################################################################
## download data, one file per individual (for entire studies see below)
#################################################################

## download all studies, and catching those where license agreement is needed
## downloading entire studies at once can become challenging as size of the data increases. Personally I always download each individual separatly. It is easier for me to work with 1000s of smallish files than a few very very heavy ones
allIDs <- all$id
# allIDs <- as.integer64(c(374990463, 3791354435))
# studyId <- allIDs[1]
start_time <- Sys.time()
results <- lapply(allIDs, function(studyId)try({
  class(studyId) <- "integer64" ## lapply transforms the class into double for some reason
  
  reftb <- movebank_download_deployment(study_id=studyId, omit_derived_data=F )
  # vroom::problems(reftb)
  reftb <- reftb[reftb$number_of_events > units::set_units(0,"count"),]

  ### sometimes connection to movebank is interrupted and indvidual does not get downloaded. When running the code again, here it checks which files exist, and which are missing
  indiv <- reftb$individual_local_identifier
  if(any(grepl("/", indiv)==T)){indiv <- gsub("/","-",indiv)}
  reftb$individual_local_identifierNObslsh <- indiv
  reftb$pthName <- paste0(studyId,"_",reftb$individual_local_identifierNObslsh,".rds")
  doneIndv <- list.files(pthDownld)[grep(studyId,list.files(pthDownld))] ## indiv from study X downloaded
  
  todoIndv <-  reftb$individual_id[!reftb$pthName%in%doneIndv]  ## using movebank internal id ("individual_id") as sometimes "individual_local_identifier" gives error (because of symbols in the name) or does not exist
  
  ## download each individual separately
  # ind <- todoIndv[1]
  lapply(todoIndv, function(ind)try({
    print(paste0("MBid: ",studyId," Indv: ",ind))
    mv2 <- movebank_download_study(studyId,
                                   sensor_type_id=c("gps","argos-doppler-shift","radio-transmitter"),
                                   # individual_local_identifier= ind,
                                   individual_id=ind,
                                   timestamp_end=as.POSIXct(Sys.time(), tz="UTC")) # to avoid locations in the future
    
    ## if one individual has several deployments, move2 automatically assigned each deployment to a different track. If you want to work on a individual base, this is how one makes sure that there is only one track per individual, independent of the number of deployments:
    # if(mt_track_id_column(mv2)=="individual_local_identifier"){mv2 <- mv2}else{
    #   mv2 <- mt_set_track_id(mv2, "individual_local_identifier")}
    
    if("individual_local_identifier" %in% names(mt_track_data(mv2))){
      indiv <- unique(mt_track_data(mv2)$individual_local_identifier)
      if(any(grepl("/", indiv)==T)){indiv <- gsub("/","-",indiv)} ## on movebank often names have "/" which messes with R
      print(indiv)
    }else{
      indiv <- mt_track_data(mv2)$individual_id
      print(indiv)
    }
    saveRDS(mv2, file=paste0(pthDownld,studyId,"_",indiv,".rds"))
  }))
}))
end_time <- Sys.time()
end_time-start_time

is.error <- function(x) inherits(x, "try-error")
table(vapply(results, is.error, logical(1)))
names(results) <- allIDs
giveError <- allIDs[vapply(results, is.error, logical(1))]



#############################
## downloading those again that need license agreement, getting the 'license-md5' for each of the studies (it is different for each study).

start_time <- Sys.time()
resultsEr <- lapply(giveError, function(studyId)try({
  class(studyId) <- "integer64" ## lapply transforms the class into double for some reason
  
  reftb <- movebank_download_deployment(study_id=studyId, 
                                        omit_derived_data=F,'license-md5'= sub('...Alternat.*','',sub('.*se-md5.=.','',as.character(rlang::catch_cnd(movebank_download_study(studyId))))))
  # vroom::problems(reftb)
  reftb <- reftb[reftb$number_of_events > units::set_units(0,"count"),]
  allindv <- unique(reftb$individual_id) ## using movebank internal id ("individual_id") as sometimes "individual_local_identifier" gives error (because of symbols in the name) or does not exist
  
  ## download each individual separatly
  # ind <- allIndv[1]
  lapply(allIndv, function(ind)try({
    print(paste0("MBid: ",studyId," Indv: ",ind))
    mv2 <- movebank_download_study(studyId,
                                   sensor_type_id=c("gps","argos-doppler-shift","radio-transmitter"),
                                   # individual_local_identifier= ind,
                                   individual_id=ind,
                                   timestamp_end=as.POSIXct(Sys.time(), tz="UTC")) # to avoid locations in the future
    
    ## if one individual has several deployments, move2 automatically assignes each deployment to a different track. If you want to work on a individual base, this is how one makes sure that there is only one track per indiviudal, independant of the number of deployments:
    # if(mt_track_id_column(mv2)=="individual_local_identifier"){mv2 <- mv2}else{
    #   mv2 <- mt_set_track_id(mv2, "individual_local_identifier")}
    
    
    if("individual_local_identifier" %in% names(mt_track_data(mv2))){
      indiv <- unique(mt_track_data(mv2)$individual_local_identifier)
      if(any(grepl("/", indiv)==T)){indiv <- gsub("/","-",indiv)} ## on movebank often names have "/" which messes with R
      print(indiv)
    }else{
      indiv <- mt_track_data(mv2)$individual_id
      print(indiv)
    }
    saveRDS(mv2, file=paste0(pthDownld,studyId,"_",indiv,".rds"))
  }))
}))
end_time <- Sys.time()
end_time-start_time

table(vapply(resultsEr, is.error, logical(1)))
names(resultsEr) <- giveError
results[vapply(resultsEr, is.error, logical(1))]
giveError2 <- giveError[vapply(results, is.error, logical(1))]



###################################
### to download one file per study 
## -> be aware if studies are very large the connection might break, or R will not be able to handle the RAM needed
####################################

start_time <- Sys.time()
results <- lapply(allIDs, function(studyId)try({
  class(studyId) <- "integer64" ## lapply transforms the class into double for some reason
  mv2 <- movebank_download_study(studyId,
                                 sensor_type_id="gps",
                                 # attributes="all", # if all attributes want to be downloaded
                                 attributes = "individual_local_identifier", # to save space, here only individual name is downloaded
                                 timestamp_end=as.POSIXct(Sys.time(), tz="UTC")) # to avoid locations in the future
  
  ## if one individual has several deployments, move2 automatically assignes each deployment to a different track. If you want to work on a individual base, this is how one makes sure that there is only one track per indiviudal, independant of the number of deployments:
  # if(mt_track_id_column(mv2)=="individual_local_identifier"){mv2 <- mv2}else{
  #   mv2 <- mt_set_track_id(mv2, "individual_local_identifier")}
  
  saveRDS(mv2, file=paste0(pthDownld,studyId,".rds"))
}))

is.error <- function(x) inherits(x, "try-error")
table(vapply(results, is.error, logical(1)))
names(results) <- allIDs
giveError <- allIDs[vapply(results, is.error, logical(1))]

## downloading those again that need license agreement, getting the 'license-md5' for each of the studies (it is different for each study).
results2 <- lapply(giveError, function(studyId)try({
  class(studyId) <- "integer64" ## lapply transforms the class into double for some reason
  mv2 <- movebank_download_study(studyId,
                                 sensor_type_id="gps",
                                 attributes = "individual_local_identifier",
                                 'license-md5'= sub('...Alternat.*','',sub('.*se-md5.=.','',as.character(rlang::catch_cnd(movebank_download_study(studyId))))),
                                 timestamp_end=as.POSIXct(Sys.time(), tz="UTC")) # to avoid locations in the future
  
  ## if one individual has several deployments, move2 automatically assignes each deployment to a different track. If you want to work on a individual base, this is how one makes sure that there is only one track per indiviudal, independant of the number of deployments:
  # if(mt_track_id_column(mv2)=="individual_local_identifier"){mv2 <- mv2}else{
  #   mv2 <- mt_set_track_id(mv2, "individual_local_identifier")}
  
  saveRDS(mv2, file=paste0(pthDownld,studyId,".rds"))
}))

table(vapply(results2, is.error, logical(1)))
names(results2) <- giveError
results[vapply(results2, is.error, logical(1))]
giveError2 <- giveError[vapply(results, is.error, logical(1))]