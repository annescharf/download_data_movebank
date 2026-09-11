#### in this script - subset to 1h and remove outliers -
## subsetting tracks to 1 fix per hour
## removing outliers based on speed
## removing outliers based on distance

#ToDo: remove 1st location as his will not be 1h apart, or adjust to the new mt_filter_lag() or similar when it exists

pathTOfolder <- "path_to_folder_where_all_data_will_be_saved/"
referenceTableStudies <- readRDS(paste0(pathTOfolder,"/referenceTableStudies_ALL_excludedColumn.rds"))
referenceTableStudiesUsed <- referenceTableStudies[referenceTableStudies$excluded=="no",]
head(referenceTableStudiesUsed)
summary(referenceTableStudiesUsed)

referenceTableStudiesUsed <- referenceTableStudiesUsed[!is.na(referenceTableStudiesUsed$species),]

pthClean <- paste0(pathTOfolder,"2.Indiv_data_clean_empty_duply/")
dir.create(paste0(pathTOfolder,"3.Indiv_data_1h"))
savePath <- paste0(pathTOfolder,"3.Indiv_data_1h/")
dir.create(paste0(pathTOfolder,"4.Indiv_data_1h_outlspeed"))
savePathOutl <- paste0(pathTOfolder,"4.Indiv_data_1h_outlspeed/")
dir.create(paste0(pathTOfolder,"5.Indiv_data_1h_outlspeed_dist"))
savePathOutlDist <- paste0(pathTOfolder,"5.Indiv_data_1h_outlspeed_dist/")


flsMV2 <- as.character(referenceTableStudiesUsed$fileName)
flsMV2[1]

library(move2)
library(units)

# ind <- flsMV2[1]

start_time<- Sys.time()
results <- lapply(flsMV2, function(ind)try({
  mv2 <- readRDS(paste0(pthClean,ind))
  mv2_thinned <- mt_filter_per_interval(mv2,criterion = "first",unit="hour")
  mv2_thinned <- mv2_thinned[-1,] ## just in case
  
  saveRDS(mv2_thinned, file=paste0(savePath,ind))
}))
end_time <- Sys.time() 
end_time-start_time # ~41mins

is.error <- function(x) inherits(x, "try-error")
table(vapply(results, is.error, logical(1)))
names(results) <- seq_along(results)
results[vapply(results, is.error, logical(1))]

library(units)
## check distribution of speeds
flsMVs <- list.files(savePath, full.names = T)
# indPth <- flsMVs[1]
start_time <- Sys.time()
speedL <- lapply(flsMVs, function(indPth){
  mv2 <- readRDS(indPth)
  mv2_speed <- mt_speed(mv2, units="m/s")
  return(mv2_speed)
})
end_time <- Sys.time()
end_time - start_time # ~10mins

saveRDS(speedL, file=paste0(pathTOfolder,"speed_all_list",".rds"))


speedL <- readRDS(paste0(pathTOfolder,"speed_all_list",".rds"))
speedAll <- unlist(speedL)
speedAll <- speedAll[!is.na(speedAll)]
hist(speedAll)
hist(speedAll[speedAll<19])
round(quantile(speedAll, seq(0.95,1,0.001)),2)

## remove speeds higher than threshold XX -- remove top XX%
library(dplyr)
library(move2)
flsMVs <- list.files(savePath, full.names = F)
# indPth <- flsMVs[1]
start_time <- Sys.time()
maxspeed <- 20
results <- lapply(flsMVs, function(indPth)try({
  print(indPth)
  mv2 <- readRDS(paste0(savePath,indPth))
  while(any(mt_speed(mv2, units="m/s")>set_units(maxspeed, m/s), na.rm = TRUE)){
    mv2 <- mv2 %>% filter(mt_speed(., units="m/s")<=set_units(maxspeed, m/s) | is.na(mt_speed(., units="m/s")))
  }
  saveRDS(mv2, file=paste0(savePathOutl,indPth))
}))
end_time <- Sys.time()
end_time - start_time # 40min

is.error <- function(x) inherits(x, "try-error")
table(vapply(results, is.error, logical(1)))
names(results) <- seq_along(results)
results[vapply(results, is.error, logical(1))]

### remove outliers based on distance
## check distribution of speeds
flsMVs <- list.files(savePathOutl, full.names = T)
# indPth <- flsMVs[1]
start_time <- Sys.time()
distL <- lapply(flsMVs, function(indPth){
  mv2 <- readRDS(indPth)
  mv2_dist <- mt_distance(mv2, units="m")
  return(mv2_dist)
})
end_time <- Sys.time()
end_time - start_time #10min

distAll <- unlist(distL)
distAll <- distAll[!is.na(distAll)]
hist(distAll)
round(quantile(distAll, seq(0.9,1,0.01)),2)
hist(distAll[distAll<50000])
round(quantile(distAll, seq(0.95,1,0.001)))

dl <- distAll[distAll>1000000]


## remove distances higher than threshold 1000K km -- remove top 0.0001%
library(dplyr)
library(move2)
flsMVs <- list.files(savePathOutl, full.names = F)
# indPth <- flsMVs[1]
start_time <- Sys.time()
maxdist <- 1000000
results <- lapply(flsMVs, function(indPth)try({
  print(indPth)
  mv2 <- readRDS(paste0(savePathOutl,indPth))
  mv2 <- mv2 %>% filter(mt_distance(., units="m")<=set_units(maxdist, m) | is.na(mt_distance(., units="m")))
  saveRDS(mv2, file=paste0(savePathOutlDist,indPth))
}))
end_time <- Sys.time()
end_time - start_time # 25min

is.error <- function(x) inherits(x, "try-error")
table(vapply(results, is.error, logical(1)))
names(results) <- seq_along(results)
results[vapply(results, is.error, logical(1))]
