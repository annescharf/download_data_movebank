## Title: Identify duplicated tracks using only locations and timestamps. 
## Author: Anne K Scharf, MPI of Animal Behavior, with the help of Claude Code
## Date: September 2026
## Description: This script sources the script `2_Functions_Duplicated_tracks_across_studies.R`
##              Here the arguments of the functions are adjusted and run. For details 
##              on the functions see head of the script `2_Functions_Duplicated_tracks_across_studies.R`
##
## INPUT: The folder you obtained in Step 3 of the script ̀`1_Download_data_from_Movebank.R,
##         i.e. one folder containing one .rds file per individual
##
## OUTPUT: (see more details below) 
##  - `reference_table_all_studies.rds`: one row per individual, containing a set 
##             of metadata information, this is the inventory of files before 
##             duplicated tracks are identified. Individuals without taxon 
##             information will be flagged here, the column `species` will be NA 
##             and these individuals will be ignored in the duplication detection. 
##            Advice: add the taxon name to the Movebank study or ask the owner of the study to do so. 
##  - `duplicated_tracks.csv`: same table as above plus the result of the duplicate 
##            detection. Used to identify which individuals to keep. 
##            keep == FALSE are the files to exclude
##  - `duplicated_tracks_pairs.csv`: one row per pair of files that was compared 
##            location by location, use the grouping to check borderline cases.
##  - `duplicated_tracks.rds`: the complete result object of find_duplicated_tracks() 
##            that will be used for plot_duplicate_groups(). Can be used to get the 
##            duplicated plots without needing to rerun all the code. 
##  - `duplicatePlots`: contains one plot per duplicated group, these are organized 
##            in one folder per species
##
##
## Details on output files:
##  reference_table_all_studies.rds  One row per rds file, written before any 
##                                    duplicate detection. The table contains these columns: 
##                                    fileName, path, MBid, individual_local_identifier,
##                                    tag_local_identifier, species, animal_life_stage, manipulation_type,
##                                    n_tracks_in_file, tracking_duration_days, tracking_start/end_date,
##                                    GPSpts_total, GPSpts_used, median/min_timelag_mins, coord_decimals,
##                                    timestamp_ms, xmin/xmax/ymin/ymax, read_error, no_species. 
##                                    The inventory of all files, independent of the duplicate analysis; 
##                                    on disk even if a later step fails.
##  duplicated_tracks.csv            The same table, one row per file, plus the result of the duplicate
##                                    detection: dup_group (id of the files belonging to the same
##                                    individual, NA = no duplicate found), n_in_group, keep (TRUE = use
##                                    this file: best of its group or no duplicate; FALSE = duplicate to
##                                    drop; NA = no species, not checked), kept_fileName (kept file of the
##                                    group), frac_in_kept (share of a dropped file's locations that are
##                                    also in the kept file; NA when only linked via a third file). Rows
##                                    ordered by species, group, kept file first. This is the operational
##                                    table: keep == FALSE are the files to exclude.
##  duplicated_tracks_pairs.csv     One row per pair of files that was compared location by location
##                                    (candidate pairs after the pre-screening; pairs never sharing a
##                                    day/cell do not appear): fileName_a, fileName_b, species,
##                                    is_duplicate, n_a_window, n_b_window (locations of each file in the
##                                    common time window), n_coincident (fixes coinciding in time),
##                                    n_matched (of those, also within tol_dist_m),
##                                    frac_matched_of_coincident (the decisive statistic: ~1 same data,
##                                    ~0 different animals), frac_a_matched, frac_b_matched (matched fixes
##                                    as share of ALL locations of a / b), dist_med/q90/max_m (distances
##                                    of the time-coincident fixes). The evidence behind the grouping;
##                                    use it to check borderline cases.
##  duplicated_tracks.rds            The complete result object of find_duplicated_tracks(): a list
##                                    with $tracks (= duplicated_tracks.csv as data.table with proper
##                                    date/number types), $pairs (= duplicated_tracks_pairs.csv),
##                                    $groups (list of the duplicate groups, vector of fileNames each,
##                                    kept file first) and $settings (the thresholds used): 
##                                        $tracks  reference table, one row per file (columns of extract_track_info()), plus the columns:
##                                                  dup_group      id of the group of files belonging to the same
##                                                                  individual (NA = no duplicate found)
##                                                  n_in_group     number of files in that group
##                                                  keep           TRUE for the file with the longest tracking duration
##                                                                  of its group (durations within keep_duration_tol_s
##                                                                  count as equal; ties: most locations, then most
##                                                                  coordinate decimals, then timestamps with ms); TRUE
##                                                                  for all files without duplicates; NA for files
##                                                                  without species
##                                                  no_species     TRUE for files not associated to a species; these are
##                                                                  excluded from the comparison
##                                                  kept_fileName  the kept file of the group
##                                                  frac_in_kept   share of the locations of a non-kept file that are
##                                                                 also in the kept file (NA when the two were not
##                                                                compared directly). A low value means the file holds
##                                                                  locations that the kept file lacks.
##                                         $pairs   one row per compared candidate pair (see compare_track_pair()):
##                                                   fileName_a, fileName_b, species, is_duplicate, n_a_window,
##                                                   n_b_window, n_coincident, n_matched, frac_matched_of_coincident,
##                                                   frac_a_matched, frac_b_matched, dist_med_m, dist_q90_m, dist_max_m
##                                        $groups  list of duplicate groups (vector of fileNames each, kept file first)
##                                        $settings the thresholds used
##                                    Input of plot_duplicate_groups(); reload it to continue without re-running.
##
##  duplicatePlots/<species>/group_<id>_<kept file>.jpg   one plot per duplicate group (sanity check)



#____________________________________________________________________________


pthDownld <- paste0(pathTOfolder,"01_MB_indv_mv2/")
dir.create(paste0(pathTOfolder,"02_MB_indv_mv2_clean"))
pthClean <- paste0(pathTOfolder,"02_MB_indv_mv2_clean/")
pathTOfolder <- "./MBdata/"
path_to_rds_folder <- paste0(pathTOfolder,"02_MB_indv_mv2_clean/")   # one move2 object (one individual) per rds file
path_to_results    <- paste0(pathTOfolder,"03_identifying_duplicate_tracks/")      # tables are written here 
n_cores            <- 4                                   # works on Linux, macOS and Windows
path_to_plots      <- file.path(path_to_results, "duplicatePlots")  # sanity check plots: one folder per species
path_to_cache      <- file.path(path_to_results, "cache")           # compact copies of the locations (~24 bytes per location), can be deleted afterwards

source("2_Functions_Duplicated_tracks_across_studies.R")
for (p in c(path_to_results, path_to_plots, path_to_cache)) dir.create(p, showWarnings = FALSE, recursive = TRUE)

## Further arguments of find_duplicated_tracks() and their defaults (details at the function):
##   tol_time_s = 1              fixes are time-coincident when their timestamps differ by <= 1 s
##   tol_dist_m = 2              time-coincident fixes are matched when <= 2 m apart; raise it
##                               (e.g. 15) when coordinates are stored with only 4 decimals
##   min_coincident = 20         a pair needs >= 20 time-coincident fixes to be judged
##   min_frac_matched = 0.9      a pair is duplicated when >= 90% of its time-coincident fixes are matched
##   prescreen_cell_deg = 0.01   grid cell (~1 km) of the (day, cell) pre-screening of candidate pairs
##   keep_criterion = "duration" file to keep per group: longest duration ("duration") or most locations ("locations")
##   keep_duration_tol_s = 60    durations within 60 s count as equal for keep_criterion = "duration"
##   max_cache_locs = 2e7        locations held in memory per worker during the pair comparison (~24 bytes each)
##   big_file_mb = 50            files above 50 MB on disk (~2-3 million locations, several GB of RAM
##   n_cores_big = 2             each while being read) are read by at most 2 workers at a time
##   stop_on_error = FALSE       unreadable files are skipped with a warning (column read_error)
##   parallel_backend = "auto"   forking on Linux/macOS, PSOCK cluster on Windows ("fork" / "psock" to force)
##   verbose = TRUE              progress messages

res <- find_duplicated_tracks(files = path_to_rds_folder, 
                              n_cores = n_cores, 
                              cache_dir = path_to_cache,
                              reference_table_file = file.path(path_to_results, "reference_table_all_studies.rds")
                              )
fwrite(res$tracks, file = file.path(path_to_results, "duplicated_tracks.csv"))        # reference table: dup_group, keep, ...
fwrite(res$pairs,  file = file.path(path_to_results, "duplicated_tracks_pairs.csv"))  # match statistics of all compared pairs
saveRDS(res, file = file.path(path_to_results, "duplicated_tracks.rds"))

## Arguments of plot_duplicate_groups() and their defaults:
##   species = NULL              plot only these species (NULL = all)
##   max_groups_per_species = Inf
##   width = 20, height = 12     jpg size in inches
##   dpi = 100                   jpg resolution
plot_duplicate_groups(result = res, 
                      out_dir = path_to_plots
                      )                          # one jpg per duplicate group
unlink(path_to_cache, recursive = TRUE)                                               # remove the cache

## Optional
res$tracks[keep == FALSE]                              # files to exclude
res$tracks[keep == FALSE & frac_in_kept < 0.95]        # excluded files holding locations the kept file lacks
res$pairs[n_coincident > 0 & is_duplicate == FALSE]    # borderline pairs worth a look

