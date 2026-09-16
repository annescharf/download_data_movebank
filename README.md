# Download data from multiple studies from Movebank

This repository contains scripts to download data from Movebank and several basic cleaning steps.
What each script does (in each scrip yo can find more detailed comments):

- `1_Download_data_from_Movebank.R`: it contains 3 steps:      
   1. creates table of studies to be downloaded. A "metadata table" is created containing the study name, contact person and license type and terms among others.      
   2. data are downloaded by individual.    
   3. basic cleaning of the data         
           
- `2_Functions_Duplicated_tracks_across_studies.R`: contains functions that are sourced in `3_Duplicated_tracks_across_studies.R`. Currently in Movebank the same track can be found in multiple studies. The duplications of tracks can be exact, but often other studies contain only part of the tracking period, or be thinned to a coarser fix frequency, or both. In this scrip these duplicated tracks are identified from the timestamps and coordinates. Only files of the same species are compared with each other; files not associated to a species are excluded. The track with either the longest duration or the highest number of gps locations is retained.

- `3_Find_duplicated_tracks_across_studies.R`: it sources the previous scrip and runs the functions. Here the values of arguments can be adjusted if needed. It will produce a table containing all individuals, in the column "Keep" all duplicated individuals to exclude are marked with `FALSE`, those individuals that do not have any individuals are also marked with `TRUE`, and those individuals without taxon information will be `NA`.
Additionally it produces plots of all duplicated tracks for visual inspection.

**File organization across all scripts**: the files for each track are always named the same, this makes it easier to run the code independently of which folder is used or if the order changes along the way. For each step a folder is created that contains in its name the main change to the data, probably wise to enumerate them, e.g. "01_MB_ind_mv2", "02_MB_ind_mv2_basic_clean", "03_MB_ind_mv2_remove_outlier", "04_MB_ind_mv2_flying_segments", etc

For **removing outliers**, see our newly developed R library `move2utils`: https://github.com/move2universe

If you have any improvements, suggestions or questions, do not hesitate to contact me (ascharf@ab.mpg.de) and also feel free to make a PR.


