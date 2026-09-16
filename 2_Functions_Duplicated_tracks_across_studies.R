## Title: Function to identify duplicated tracks using only locations and timestamps. 
## Author: Anne K Scharf, MPI of Animal Behavior, with the help of Claude Code
## Date: September 2026
## Description: This script contains the functions that are sourced in the the script
##             `3_Find_duplicated_tracks_across_studies.R`
## Background: Currently in Movebank the same track can be found in multiple studies. 
##             The duplications of tracks can be exact, but often other studies 
##             contain only part of the tracking period, or be thinned to a coarser 
##             fix frequency, or both. In this scrip these duplicated tracks are 
##             identified from the timestamps and coordinates. Only files of the 
##             same species are compared with each other; files not associated to 
##             a species are excluded. The track with either the longest duration 
##             or the highest number of gps locations is retained. 
##
## INPUT: The folder you obtained in Step 3 of the script ̀`1_Download_data_from_Movebank.R,
##         i.e. one folder containing one .rds file per individual
##
## OUTPUT: (see more details in the script `3_Find_duplicated_tracks_across_studies.R`) 
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
##            duplicated plots without needing to rerun all the code
##  - `duplicatePlots`: contains one plot per duplicated group, these are organized 
##            in one folder per species
##
## Details on the script:
## Workflow of find_duplicated_tracks():
##  1. Reference table, one row per file (study, individual, tag, species, life
##     stage, manipulation, tracking period, number of locations, time lags,
##     precision, bounding box; see extract_track_info()). Meanwhile the
##     locations are reduced to (time, lon, lat) and cached on disk.
##  2. For every species: candidate pairs = files that share at least one
##     (day, ~1 km cell) and whose tracking periods overlap.
##  3. Candidate pairs are compared location by location (see below).
##  4. Groups of files belonging to the same individual = connected components
##     of the duplicated pairs; per group the file with the longest tracking
##     duration (default, see keep_criterion) is marked to be kept.
##
## Rationale of the location comparison:
##  - Copies of the same data have (near) identical timestamps and coordinates.
##    "Near" because timestamps may be rounded differently (+-1 s from dropped
##    milliseconds) and coordinates stored with a different number of decimals.
##  - Two copies with different fix frequencies (e.g. 5 vs 20 min) only share a
##    fraction of the fixes in time. Therefore the decisive statistic is:
##        among the fixes of two tracks that coincide in time (+- tol_time_s),
##        which fraction also coincide in space (<= tol_dist_m)?
##    This is ~1 for copies of the same data and ~0 for different animals,
##    independent of the fix frequency of either copy.
##  - Duplication is transitive but not necessarily pairwise detectable: two
##    subsets of the same original can cover disjoint months (both are then
##    duplicates of the full track, but share nothing with each other). Groups
##    are therefore built as connected components of the pairwise matches.
##
## Scalability: tracks can have a few hundred to several million locations,
## the compact cache holds 24 bytes per location, comparisons are rolling
## nearest-time joins (data.table, O(n log n)), steps 1 and 3 can run on several
## cores (forking on Linux/macOS, a PSOCK cluster on Windows; see parallel_apply()).

#____________________________________________________________________________

library(move2)
library(sf)
library(data.table)

#____________________________________________________________________________
## Step 1: read one move2 file -> reference table row + compact locations ####

# Reduce a move2 file to a data.table of (t, x, y): numeric seconds since epoch
# (time-zone independent), lon, lat; sorted by time, exact duplicate rows
# removed. All tracks of the file (an individual with several tags/deployments)
# are pooled. Returns the table with the move2 object's track data, number of
# rows and number of tracks as attributes.
read_track_locations <- function(path) {
  m <- readRDS(path)
  if (!inherits(m, "move2")) stop("Not a move2 object: ", path)
  # Work on the coordinate matrix, never subset the move2 object (very slow for
  # millions of rows): empty points give NA rows in st_coordinates, so rows stay
  # aligned with the timestamps and are dropped afterwards.
  tt <- mt_time(m)
  if (inherits(tt, "Date")) tt <- as.POSIXct(tt, tz = "UTC")
  xy <- st_coordinates(st_geometry(m))[, 1:2, drop = FALSE]
  ok <- is.finite(xy[, 1]) & is.finite(xy[, 2]) & !is.na(tt)
  # lon/lat coordinates; sf_project on the matrix is much faster than st_transform.
  # Objects without CRS are assumed to be lon/lat.
  if (any(ok) && isFALSE(sf::st_is_longlat(m))) {
    xy[ok, ] <- sf::sf_project(from = st_crs(m), to = st_crs(4326), pts = xy[ok, , drop = FALSE])
  }
  loc <- unique(data.table(t = as.numeric(tt)[ok], x = xy[ok, 1], y = xy[ok, 2]))
  setkey(loc, t)
  setattr(loc, "track_data", mt_track_data(m))
  setattr(loc, "n_locs", nrow(m))
  setattr(loc, "n_tracks", mt_n_tracks(m))
  loc
}

# Reference table row of one file. Writes the compact (t, x, y) table to
# cache_dir/<fileName without .rds>.rds (uncompressed for fast re-reading).
# Track attributes that differ between the tracks of a file are concatenated with ";".
# Columns: fileName, path, track_id (fileName without .rds, internal), MBid
# (study_id), individual_local_identifier, tag_local_identifier (tag_id if
# missing), species (taxon_canonical_name), animal_life_stage, manipulation_type
# (all three NA when missing), n_tracks_in_file,
# tracking_duration_days, tracking_start_date, tracking_end_date (UTC),
# GPSpts_total (rows of the move2 object), GPSpts_used (rows with location and
# timestamp, exact duplicate rows removed), median_timelag_mins, min_timelag_mins,
# coord_decimals, timestamp_ms (precision of the stored coordinates/timestamps,
# from a systematic sample of <= 1000 locations), xmin, xmax, ymin, ymax
# (bounding box, lon/lat), keys (internal, for the pre-screening).
extract_track_info <- function(path, cache_dir, prescreen_cell_deg = 0.01) {
  fileName <- basename(path)
  track_id <- sub("\\.rds$", "", fileName, ignore.case = TRUE)
  loc <- read_track_locations(path)
  n_locs <- attr(loc, "n_locs")
  n_tracks <- attr(loc, "n_tracks")
  td <- attr(loc, "track_data")
  # track attribute as character: NA when the column is missing or holds only
  # NA/empty values; several distinct values (deployments) joined with ";"
  attr_chr <- function(col) {
    if (!col %in% names(td)) return(NA_character_)
    v <- unique(as.character(td[[col]]))
    v <- v[!is.na(v) & v != ""]
    if (length(v) == 0) NA_character_ else paste(v, collapse = ";")
  }
  # tag: fall back to tag_id when tag_local_identifier is missing
  tag <- attr_chr("tag_local_identifier")
  if (is.na(tag)) tag <- attr_chr("tag_id")
  saveRDS(loc, file = file.path(cache_dir, paste0(track_id, ".rds")), compress = FALSE)

  n_used <- nrow(loc)
  lags <- if (n_used > 1) diff(loc$t) / 60 else numeric(0)
  # (day, grid cell) keys for the pre-screening of candidate pairs
  keys <- unique(data.table(day = floor(loc$t / 86400), cx = floor(loc$x / prescreen_cell_deg), cy = floor(loc$y / prescreen_cell_deg)))
  # precision of the stored data (tie-breakers when choosing which copy to keep)
  smp <- loc[unique(round(seq(1, n_used, length.out = 1000)))]   # systematic sample of <= 1000 locations
  n_decimals <- function(v) { s <- sub("0+$", "", sub("^[^.]*\\.?", "", formatC(v, digits = 9, format = "f"))); nchar(s) }
  coord_decimals <- if (n_used) max(n_decimals(smp$x), n_decimals(smp$y)) else NA_integer_
  timestamp_ms <- if (n_used) round(mean(smp$t %% 1 != 0), 2) else NA_real_
  data.table(fileName = fileName, 
             path = path, 
             track_id = track_id,
             MBid = attr_chr("study_id"),
             individual_local_identifier = attr_chr("individual_local_identifier"),
             tag_local_identifier = tag,
             species = attr_chr("taxon_canonical_name"),
             animal_life_stage = attr_chr("animal_life_stage"),   # not used for filtering, needed afterwards
             manipulation_type = attr_chr("manipulation_type"),   # "none" or NA are not manipulated
             n_tracks_in_file = n_tracks,
             tracking_duration_days = if (n_used) round((loc$t[n_used] - loc$t[1]) / 86400, 2) else NA_real_,
             tracking_start_date = if (n_used) as.POSIXct(loc$t[1], origin = "1970-01-01", tz = "UTC") else as.POSIXct(NA, tz = "UTC"),
             tracking_end_date = if (n_used) as.POSIXct(loc$t[n_used], origin = "1970-01-01", tz = "UTC") else as.POSIXct(NA, tz = "UTC"),
             GPSpts_total = n_locs,
             GPSpts_used = n_used,
             median_timelag_mins = if (length(lags)) round(median(lags)) else NA_real_,
             min_timelag_mins = if (length(lags)) round(min(lags), 2) else NA_real_,
             coord_decimals = coord_decimals,
             timestamp_ms = timestamp_ms,
             xmin = if (n_used) min(loc$x) else NA_real_, xmax = if (n_used) max(loc$x) else NA_real_,
             ymin = if (n_used) min(loc$y) else NA_real_, ymax = if (n_used) max(loc$y) else NA_real_,
             keys = list(keys))
}

#____________________________________________________________________________
## Helpers ####

# approximate distance in metres between lon/lat points (equirectangular; exact
# enough for the metre-scale tolerances used here), robust to the date line
geo_dist_m <- function(x1, y1, x2, y2) {
  r <- 6371008.8
  dlon <- ((x1 - x2 + 180) %% 360 - 180) * pi / 180
  dlat <- (y1 - y2) * pi / 180
  sqrt((dlon * cos((y1 + y2) * pi / 360))^2 + dlat^2) * r
}

# in-memory cache of compact tracks, so the same file is not re-read for every
# pair it is involved in; emptied when it holds more than max_locs locations
make_track_loader <- function(cache_dir, max_locs = 2e7) {
  store <- new.env(parent = emptyenv())
  held <- 0
  function(track_id) {
    if (!exists(track_id, envir = store, inherits = FALSE)) {
      loc <- readRDS(file.path(cache_dir, paste0(track_id, ".rds")))
      if (held + nrow(loc) > max_locs) { rm(list = ls(store), envir = store); held <<- 0 }
      assign(track_id, loc, envir = store)
      held <<- held + nrow(loc)
    }
    get(track_id, envir = store, inherits = FALSE)
  }
}

# union-find: connected components of an undirected graph given as edge list
connected_components <- function(nodes, from, to) {
  parent <- setNames(seq_along(nodes), nodes)
  find <- function(i) { while (parent[i] != i) { parent[i] <<- parent[parent[i]]; i <- parent[i] }; i }
  for (k in seq_along(from)) {
    a <- find(match(from[k], nodes)); b <- find(match(to[k], nodes))
    if (a != b) parent[b] <- a
  }
  roots <- vapply(seq_along(nodes), find, integer(1))
  as.integer(factor(roots, levels = unique(roots)))
}

#____________________________________________________________________________
## Step 2: candidate pairs, cheap pre-screening ####
# A pair of tracks is a candidate when both have at least one location in the
# same (day, grid cell) - keys computed in step 1 - and their tracking periods
# overlap. Pairs are found with a self join of the keys, i.e. without looping
# over all pairs of tracks.
find_candidate_pairs <- function(tracks, tol_time_s, max_join_rows = 5e6) {
  empty <- data.table(track_a = character(0), track_b = character(0))
  tr <- tracks[GPSpts_used > 0]
  if (nrow(tr) < 2) return(empty)
  keys <- rbindlist(Map(function(k, id) data.table(k, track_id = id), tr$keys, tr$track_id))
  keys[, n := .N, by = .(day, cx, cy)]
  keys <- keys[n > 1]                    # keys of a single track cannot give pairs
  if (nrow(keys) == 0) return(empty)
  # A key shared by n tracks (e.g. a roost used by hundreds of animals every
  # day) gives n^2 join rows, so the join is done in chunks of bounded size.
  # Over the rows of one key, cumsum(n) grows by n^2; a key never spans chunks.
  setkey(keys, day, cx, cy)
  keys[, chunk := floor(cumsum(as.numeric(n)) / max_join_rows)]
  keys[, chunk := chunk[1], by = .(day, cx, cy)]
  shared <- rbindlist(lapply(split(keys, keys$chunk), function(kl) {
    unique(kl[kl, on = .(day, cx, cy), allow.cartesian = TRUE,
              .(track_a = x.track_id, track_b = i.track_id)][track_a < track_b])
  }))
  shared <- unique(shared)
  # exact overlap of the tracking periods (sharing a day is coarser than that)
  period <- tr[, .(track_id, t0 = as.numeric(tracking_start_date) - tol_time_s, t1 = as.numeric(tracking_end_date) + tol_time_s)]
  shared[period, `:=`(t0_a = i.t0, t1_a = i.t1), on = .(track_a = track_id)]
  shared[period, `:=`(t0_b = i.t0, t1_b = i.t1), on = .(track_b = track_id)]
  shared[t0_a <= t1_b & t1_a >= t0_b, .(track_a, track_b)][order(track_a, track_b)]
}

#____________________________________________________________________________
## Step 3: compare one candidate pair location by location ####
# For every fix of the sparser track (within the common time window) the nearest
# fix in time of the denser track is found. Fixes closer than tol_time_s are
# "time-coincident"; those additionally closer than tol_dist_m are "matched".
# Returns one row: n_a_window, n_b_window (locations of a and b in the common
# time window), n_coincident, n_matched, frac_matched_of_coincident (the decisive
# statistic), frac_a_matched, frac_b_matched (share of ALL locations of a / b
# that are matched), dist_med_m, dist_q90_m, dist_max_m (distance between the
# time-coincident fixes).
compare_track_pair <- function(a, b, tol_time_s = 1, tol_dist_m = 2) {
  na_tot <- nrow(a); nb_tot <- nrow(b)
  t0 <- max(a$t[1], b$t[1]) - tol_time_s
  t1 <- min(a$t[na_tot], b$t[nb_tot]) + tol_time_s
  aw <- a[t >= t0 & t <= t1]
  bw <- b[t >= t0 & t <= t1]
  na_w <- nrow(aw); nb_w <- nrow(bw)
  out <- data.table(n_a_window = na_w, n_b_window = nb_w, n_coincident = 0L, n_matched = 0L,
                    frac_matched_of_coincident = NA_real_,
                    frac_a_matched = 0, frac_b_matched = 0,
                    dist_med_m = NA_real_, dist_q90_m = NA_real_, dist_max_m = NA_real_)
  if (na_w == 0 || nb_w == 0) return(out)

  # query with the sparser track, index the denser one. The denser track can
  # hold several fixes with the same timestamp (different coordinates), in which
  # case the join returns all of them: a query fix counts as coincident/matched
  # when any of them is.
  if (na_w <= nb_w) { S <- aw; D <- bw } else { S <- bw; D <- aw }
  S[, iS := seq_len(.N)]
  j <- D[S, on = "t", roll = "nearest", nomatch = NULL,
         .(iS = i.iS, dt = abs(i.t - x.t), d = geo_dist_m(i.x, i.y, x.x, x.y))]
  j <- j[dt <= tol_time_s]
  n_co <- uniqueN(j$iS)
  out[, n_coincident := n_co]
  if (n_co == 0) return(out)
  d <- j[, .(d = min(d)), by = iS]$d          # per coincident query fix: closest same-time fix
  n_m <- sum(d <= tol_dist_m)
  out[, `:=`(n_matched = n_m,
             frac_matched_of_coincident = n_m / n_co,
             frac_a_matched = n_m / na_tot,   # share of ALL locations of a that are also in b
             frac_b_matched = n_m / nb_tot,
             dist_med_m = as.numeric(median(d)), dist_q90_m = as.numeric(quantile(d, 0.9)),
             dist_max_m = max(d))]
  out
}

#____________________________________________________________________________
## Worker functions of the parallel steps ####
# Top-level functions with explicit arguments, so that they can be sent to
# parallel workers on all operating systems (see parallel_apply()).

# step 1: reference table row of one file; an unreadable file gives a row with
# read_error filled in (or an error when stop_on_error)
read_track_info_safe <- function(path, cache_dir, prescreen_cell_deg, stop_on_error) {
  r <- try(extract_track_info(path = path, cache_dir = cache_dir, prescreen_cell_deg = prescreen_cell_deg), silent = TRUE)
  if (inherits(r, "try-error")) {
    if (stop_on_error) stop("Error reading ", path, ": ", r)
    return(data.table(fileName = basename(path), path = path, track_id = sub("\\.rds$", "", basename(path), ignore.case = TRUE),
                      GPSpts_used = 0L, read_error = trimws(as.character(r))))
  }
  r[, read_error := NA_character_]
  r
}

# step 3: compare the candidate pairs in rows of cand (all with the same
# track_a). load_track: a track loader from make_track_loader(); on a PSOCK
# worker it is the loader created in the worker's global environment
compare_pair_rows <- function(rows, cand, tol_time_s, tol_dist_m, load_track = get("load_track", envir = globalenv())) {
  rbindlist(lapply(rows, function(k) {
    compare_track_pair(a = load_track(cand$track_a[k]), b = load_track(cand$track_b[k]),
                       tol_time_s = tol_time_s, tol_dist_m = tol_dist_m)
  }))
}

# Parallel lapply on all operating systems (plain lapply when n_cores is 1).
# backend "fork" (parallel::mclapply, Linux/macOS): workers share the memory of
# the main process; preschedule = FALSE hands out the elements one by one (for
# few, unequal tasks). backend "psock" (parallel::parLapplyLB, all systems,
# needed on Windows): separate R processes that load the packages and receive
# the functions of this script; the elements are always handed out one by one
# (the tasks are unequal). Extra arguments to FUN are sent to the workers with
# every element (keep them small); if they include cache_dir and max_cache_locs,
# these are not passed on but used to create one track cache (load_track) in the
# global environment of every worker, which compare_pair_rows() picks up.
parallel_apply <- function(X, FUN, ..., n_cores, backend, preschedule = TRUE) {
  if (n_cores <= 1 || length(X) == 0) return(lapply(X, FUN, ...))
  if (backend == "fork") return(parallel::mclapply(X, FUN, ..., mc.cores = n_cores, mc.preschedule = preschedule))
  cl <- parallel::makeCluster(min(n_cores, length(X)))
  on.exit(parallel::stopCluster(cl))
  # one data.table thread per worker, otherwise the workers oversubscribe the cores
  parallel::clusterEvalQ(cl, { suppressPackageStartupMessages({ library(move2); library(sf); library(data.table) }); setDTthreads(1) })
  parallel::clusterExport(cl, c("read_track_locations", "extract_track_info", "read_track_info_safe", "geo_dist_m",
                                "make_track_loader", "compare_track_pair", "compare_pair_rows"),
                          envir = environment(find_duplicated_tracks))
  extra <- list(...)
  if (all(c("cache_dir", "max_cache_locs") %in% names(extra))) {   # step 3: one track cache per worker
    parallel::clusterCall(cl, function(cache_dir, max_locs) {
      assign("load_track", make_track_loader(cache_dir = cache_dir, max_locs = max_locs), envir = globalenv()); NULL
    }, extra$cache_dir, extra$max_cache_locs)
    extra$cache_dir <- NULL; extra$max_cache_locs <- NULL
  }
  do.call(parallel::parLapplyLB, c(list(cl, X, FUN), extra, list(chunk.size = 1)))
}

#____________________________________________________________________________
## Main function ####
# Only files of the same species (taxon_canonical_name) are compared; files
# without species are excluded (they stay in the table with keep = NA).
# files             character vector of rds files (one move2 object, one individual each)
#                   or a single directory (all *.rds files in it, not recursive)
# tol_time_s        two fixes are time-coincident when their timestamps differ by <= this
# tol_dist_m        time-coincident fixes are matched when they are <= this far apart (m)
# min_coincident    minimum number of time-coincident fixes needed to judge a pair
# min_frac_matched  pair is duplicated when >= this fraction of the time-coincident
#                   fixes are matched in space
# prescreen_cell_deg grid cell size (degrees) of the (day, cell) pre-screening
# keep_criterion    which file of a duplicate group to mark as kept:
#                   "duration" (longest tracking duration, ties -> most locations)
#                   or "locations" (most locations, ties -> longest duration);
#                   remaining ties -> most coordinate decimals, then timestamps
#                   with ms, then file name
# keep_duration_tol_s  durations differing by less than this (seconds) count as
#                   equal when choosing the file to keep
# cache_dir         where compact (t, x, y) tables are written (24 bytes per
#                   location); a folder in tempdir() by default
# max_cache_locs    max locations kept in memory per worker during step 3 (~24 bytes each)
# n_cores           cores for steps 1 and 3
# parallel_backend  "auto": forking (mclapply) on Linux/macOS, a PSOCK cluster
#                   (parLapplyLB) on Windows; can be forced to "fork" or "psock"
# big_file_mb, n_cores_big  files larger than big_file_mb (MB on disk, ~50 MB is
#                   roughly 2-3 million locations) are read with at most
#                   n_cores_big workers at a time, as each needs several GB of RAM
# stop_on_error     FALSE: a file that cannot be read gets a row with read_error
#                   filled in and is skipped (with a warning); TRUE: stop
# reference_table_file  if given, the reference table of step 1 (all rows of
#                   extract_track_info() joined, before any duplicate columns)
#                   is saved there as rds right after step 1
# verbose           print progress messages
find_duplicated_tracks <- function(files,
                                   tol_time_s = 1,
                                   tol_dist_m = 2,
                                   min_coincident = 20,
                                   min_frac_matched = 0.9,
                                   prescreen_cell_deg = 0.01,
                                   keep_criterion = c("duration", "locations"),
                                   keep_duration_tol_s = 60,
                                   cache_dir = file.path(tempdir(), "dup_tracks_cache"),
                                   max_cache_locs = 2e7,
                                   n_cores = 1,
                                   parallel_backend = c("auto", "fork", "psock"),
                                   big_file_mb = 50,
                                   n_cores_big = 2,
                                   stop_on_error = FALSE,
                                   reference_table_file = NULL,
                                   verbose = TRUE) {
  keep_criterion <- match.arg(keep_criterion)
  if (length(files) == 1 && dir.exists(files)) files <- list.files(files, pattern = "\\.rds$", full.names = TRUE, ignore.case = TRUE)
  if (length(files) == 0) stop("No rds files found")
  if (!all(file.exists(files))) stop("File(s) not found: ", paste(head(files[!file.exists(files)]), collapse = ", "))
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  if (anyDuplicated(basename(files))) stop("File names must be unique: ",
                                           paste(unique(basename(files)[duplicated(basename(files))]), collapse = ", "))
  msg <- function(...) if (verbose) message(format(Sys.time(), "%H:%M:%S"), " ", ...)
  n_cores <- max(1L, as.integer(n_cores))
  parallel_backend <- match.arg(parallel_backend)
  if (parallel_backend == "auto") parallel_backend <- if (.Platform$OS.type == "windows") "psock" else "fork"
  if (parallel_backend == "fork" && .Platform$OS.type == "windows") stop("Forking is not available on Windows, use parallel_backend = \"psock\"")

  ## step 1: reference table + compact locations ------------------------------
  msg("Step 1: reading ", length(files), " files (", parallel_backend, ", ", n_cores, " cores)")
  read_args <- list(cache_dir = cache_dir, prescreen_cell_deg = prescreen_cell_deg, stop_on_error = stop_on_error)
  read_files <- function(f, n, preschedule = TRUE) {
    do.call(parallel_apply, c(list(f, read_track_info_safe), read_args,
                              list(n_cores = n, backend = parallel_backend, preschedule = preschedule)))
  }
  # big files (millions of locations) need several GB of RAM each while being
  # read, so at most n_cores_big workers read them (biggest first)
  files <- files[order(-file.size(files))]
  big <- file.size(files) > big_file_mb * 1e6
  n_big <- min(n_cores, n_cores_big)
  if (any(big) && n_cores > 1 && parallel_backend == "fork") {
    # forking: the big files are read in a background job while the remaining
    # cores read the small files
    job_big <- parallel::mcparallel(read_files(files[big], n_big, preschedule = FALSE))
    tracks_small <- read_files(files[!big], max(1L, n_cores - n_big))
    tracks_big <- parallel::mccollect(job_big)[[1]]
    if (inherits(tracks_big, "try-error")) stop("Error reading big files: ", tracks_big)
    if (length(tracks_big) != sum(big)) stop("Reading the ", sum(big), " big files failed (worker died, probably out of memory): ",
                                            "reduce n_cores_big or n_cores")
  } else {
    # PSOCK or one core: big files first, then the small ones
    tracks_big <- read_files(files[big], n_big, preschedule = FALSE)
    tracks_small <- read_files(files[!big], n_cores)
  }
  tracks <- rbindlist(c(tracks_big, tracks_small), fill = TRUE)
  if (any(!is.na(tracks$read_error))) warning(sum(!is.na(tracks$read_error)), " file(s) could not be read, see column read_error")
  if (!"species" %in% names(tracks)) tracks[, species := NA_character_]
  msg("  ", sum(tracks$GPSpts_used), " locations in total; ",
      length(unique(na.omit(tracks$species))), " species, ", sum(is.na(tracks$species)), " files without species")
  # exclude tracks not associated to a species
  tracks[, no_species := is.na(species)]
  if (!is.null(reference_table_file)) {
    saveRDS(tracks[, setdiff(names(tracks), c("track_id", "keys")), with = FALSE], file = reference_table_file)
    msg("  reference table saved to ", reference_table_file)
  }
  if (any(tracks$no_species)) warning(sum(tracks$no_species), " file(s) without species (taxon_canonical_name) are excluded from the comparison, see column no_species")

  ## step 2: candidate pairs (per species) ---------------------------------------
  cand <- rbindlist(lapply(unique(na.omit(tracks$species)), function(sp) {
    find_candidate_pairs(tracks = tracks[species == sp], tol_time_s = tol_time_s)
  }))
  msg("Step 3: comparing ", nrow(cand), " candidate pairs")

  ## step 3: compare candidate pairs ------------------------------------------
  if (nrow(cand) > 0) {
    # chunk by track_a so a worker mostly re-uses tracks it already loaded
    chunks <- split(seq_len(nrow(cand)), cand$track_a)
    load_track <- make_track_loader(cache_dir = cache_dir, max_locs = max_cache_locs)   # used by fork/serial workers
    stats <- if (parallel_backend == "psock" && n_cores > 1) {
      parallel_apply(chunks, compare_pair_rows, cand = cand, tol_time_s = tol_time_s, tol_dist_m = tol_dist_m,
                     cache_dir = cache_dir, max_cache_locs = max_cache_locs, n_cores = n_cores, backend = "psock")
    } else {
      parallel_apply(chunks, compare_pair_rows, cand = cand, tol_time_s = tol_time_s, tol_dist_m = tol_dist_m,
                     load_track = load_track, n_cores = n_cores, backend = parallel_backend)
    }
    err <- vapply(stats, inherits, logical(1), what = "try-error")
    if (any(err)) stop("Error comparing pairs of ", paste(head(names(chunks)[err]), collapse = ", "), ": ", stats[[which(err)[1]]])
    stats <- rbindlist(stats)
    pairs <- cbind(cand[unlist(chunks)], stats)
    pairs[, is_duplicate := n_coincident >= min_coincident & frac_matched_of_coincident >= min_frac_matched]
    pairs[is.na(is_duplicate), is_duplicate := FALSE]
  } else {
    pairs <- data.table(track_a = character(0), track_b = character(0), n_a_window = integer(0), n_b_window = integer(0),
                        n_coincident = integer(0), n_matched = integer(0), frac_matched_of_coincident = numeric(0),
                        frac_a_matched = numeric(0), frac_b_matched = numeric(0), dist_med_m = numeric(0),
                        dist_q90_m = numeric(0), dist_max_m = numeric(0), is_duplicate = logical(0))
  }
  # report pairs by fileName
  pairs[, `:=`(fileName_a = tracks$fileName[match(track_a, tracks$track_id)],
               fileName_b = tracks$fileName[match(track_b, tracks$track_id)],
               species = tracks$species[match(track_a, tracks$track_id)])]
  pairs[, c("track_a", "track_b") := NULL]
  setcolorder(pairs, c("fileName_a", "fileName_b", "species", "is_duplicate"))
  setorder(pairs, fileName_a, fileName_b)

  ## step 4: groups (connected components of duplicated pairs) -----------------
  dup <- pairs[is_duplicate == TRUE]
  tracks[, dup_group := NA_integer_]
  if (nrow(dup) > 0) {
    nodes <- unique(c(dup$fileName_a, dup$fileName_b))
    comp <- connected_components(nodes = nodes, from = dup$fileName_a, to = dup$fileName_b)
    tracks[match(nodes, fileName), dup_group := comp]
  }
  tracks[, n_in_group := if (is.na(dup_group[1])) 1L else .N, by = dup_group]

  # which file of a group to keep
  tracks[, `:=`(keep = TRUE, kept_fileName = NA_character_, frac_in_kept = NA_real_)]
  if (nrow(dup) > 0) {
    # durations within keep_duration_tol_s count as equal (timestamps rounded to
    # the second can make a copy 1 s "longer" than the original); ties: more
    # locations, then more coordinate decimals, then timestamps with ms
    pick_keep <- function(dur_s, n, dec, ms, fn) {
      if (keep_criterion == "duration") {
        cand <- which(dur_s >= max(dur_s) - keep_duration_tol_s)
        o <- order(-n[cand], -dec[cand], -ms[cand], fn[cand])
      } else {
        cand <- which(n == max(n))
        o <- order(-dur_s[cand], -dec[cand], -ms[cand], fn[cand])
      }
      seq_along(dur_s) == cand[o][1]
    }
    tracks[!is.na(dup_group), keep := pick_keep(dur_s = as.numeric(tracking_end_date) - as.numeric(tracking_start_date),
                                               n = GPSpts_used, dec = coord_decimals, ms = timestamp_ms, fn = fileName),
           by = dup_group]
    tracks[!is.na(dup_group), kept_fileName := fileName[keep], by = dup_group]
    # share of a discarded file's locations that are in the kept file (from the direct comparison)
    frac <- rbind(pairs[, .(fileName = fileName_a, kept_fileName = fileName_b, frac = frac_a_matched)],
                  pairs[, .(fileName = fileName_b, kept_fileName = fileName_a, frac = frac_b_matched)])
    tracks[frac, frac_in_kept := i.frac, on = .(fileName, kept_fileName)]
  }
  tracks[no_species == TRUE, keep := NA]
  setorder(tracks, species, dup_group, -keep, fileName, na.last = TRUE)
  tracks[, c("track_id", "keys") := NULL]
  setcolorder(tracks, c("fileName", "dup_group", "n_in_group", "keep", "kept_fileName", "frac_in_kept", "no_species",
                        "MBid", "individual_local_identifier", "tag_local_identifier", "species", "animal_life_stage",
                        "manipulation_type",
                        "tracking_duration_days", "tracking_start_date", "tracking_end_date",
                        "GPSpts_total", "GPSpts_used", "median_timelag_mins", "min_timelag_mins",
                        "coord_decimals", "timestamp_ms"))

  groups <- if (nrow(dup) > 0) split(tracks[!is.na(dup_group)]$fileName, tracks[!is.na(dup_group)]$dup_group) else list()
  msg("Done: ", length(groups), " duplicate groups involving ", sum(!is.na(tracks$dup_group)), " of ", nrow(tracks), " files")

  list(tracks = tracks[], pairs = pairs[], groups = groups,
       settings = list(tol_time_s = tol_time_s, tol_dist_m = tol_dist_m, min_coincident = min_coincident,
                       min_frac_matched = min_frac_matched, prescreen_cell_deg = prescreen_cell_deg,
                       keep_criterion = keep_criterion))
}

#____________________________________________________________________________
## Sanity check: plot the duplicate groups ####
# One folder per species in out_dir, one jpg per duplicate group
# (group_<id>_<kept file>.jpg): the tracks of the group side by side (one facet
# per file), coloured by whether the file is kept or a duplicate, plus a footer
# with the key numbers per file. Tracks are thinned to one location per day to
# keep the plots light.
# result                  output of find_duplicated_tracks()
# out_dir                 folder for the species folders (created if needed)
# species                 plot only these species (NULL = all)
# max_groups_per_species  plot at most this many groups per species
# width, height, dpi      size (inches) and resolution of the jpgs
# Returns the paths of the jpgs (invisibly).
plot_duplicate_groups <- function(result, out_dir = ".", species = NULL, max_groups_per_species = Inf,
                                  width = 20, height = 12, dpi = 100) {
  library(ggplot2)
  tr <- result$tracks[!is.na(dup_group)]
  if (!is.null(species)) tr <- tr[species %in% ..species]
  if (nrow(tr) == 0) { message("No duplicate groups to plot"); return(invisible(character(0))) }
  # one location per day, as an sf line (or point) per file
  daily_line <- function(path) {
    loc <- read_track_locations(path)
    d <- loc[, .SD[1], by = .(day = floor(t / 86400))]
    g <- if (nrow(d) > 1) st_linestring(as.matrix(d[, .(x, y)])) else st_point(c(d$x, d$y))
    st_sfc(g, crs = 4326)
  }
  files <- character(0)
  for (sp in unique(tr$species)) {
    sp_dir <- file.path(out_dir, gsub("[^A-Za-z0-9]+", "_", sp))
    dir.create(sp_dir, showWarnings = FALSE, recursive = TRUE)
    groups <- head(unique(tr[species == sp]$dup_group), max_groups_per_species)
    for (g in groups) {
      grp <- tr[dup_group == g][order(-keep, fileName)]
      jpg_file <- file.path(sp_dir, sprintf("group_%d_%s.jpg", g, sub("\\.rds$", "", grp$fileName[grp$keep], ignore.case = TRUE)))
      res <- try({
        lines <- st_sf(fileName = factor(grp$fileName, levels = grp$fileName),
                       keep = ifelse(grp$keep, "keep", "duplicate"),
                       geometry = do.call(c, lapply(grp$path, daily_line)))
        info <- grp[, sprintf("%s: %s | %s | %.0f days | %d locations | %s min | in kept: %s", fileName,
                              individual_local_identifier, tag_local_identifier, tracking_duration_days, GPSpts_used,
                              median_timelag_mins, ifelse(is.na(frac_in_kept), "-", sprintf("%.0f%%", 100 * frac_in_kept)))]
        pl <- ggplot() + geom_sf(data = lines, aes(color = keep)) +
          facet_wrap(~fileName, nrow = 1) +
          scale_color_manual(values = c(keep = "#1b9e77", duplicate = "#d95f02")) +
          labs(title = sprintf("%s - duplicate group %d - %d files (kept: %s)", sp, g, nrow(grp), grp$fileName[grp$keep]),
               caption = paste(info, collapse = "\n"), color = NULL) +
          theme(plot.caption = element_text(hjust = 0))
        ggsave(filename = jpg_file, plot = pl, width = width, height = height, dpi = dpi, device = "jpeg")
        files <- c(files, jpg_file)
      }, silent = TRUE)
      if (inherits(res, "try-error")) warning("Group ", g, " could not be plotted: ", res)
    }
  }
  invisible(files)
}
