# Copyright (C) 2026 ACTRIS ARES DC Unit services@actris-ares.eu
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

library(isotone)
library(ncdf4)
library(radiant.data)
library(dplyr)
library(tidyr)


# Function for closing all connection and files 
closeAllConnections <- function() {
  open_files <- showConnections(all = TRUE)
  for (i in seq_len(nrow(open_files))) {
    if (open_files[i, "description"] != "stdin" && open_files[i, "description"] != "stdout" && open_files[i, "description"] != "stderr") {
      close(getConnection(as.integer(open_files[i, "class"])))
    }
  }
  cat("Tutti i file aperti sono stati chiusi.\n")
}


################################# Main Functions #################################
len_0 <- function(x)
  length(x[x > 0])

cnt <- function(x)
  length(x[!is.na(x) & x > 0])

cnt_int <- function(x) 
  sum(!is.na(x))

unq <- function(x)
  length(unique(x[!is.na(x)]))

sson1 <- function(x) {
  x <- as.numeric(x)
  if (x %in% c(3, 4, 5))
    "MarAprMay"
  else if (x %in% c(6, 7, 8))
    "JunJulAug"
  else if (x %in% c(9, 10, 11))
    "SepOctNov"
  else if (x %in% c(1, 2, 12))
    "DecJanFeb"
}

sson2 <- function(m, y, s) {
  m <- as.numeric(m)
  y <- as.numeric(y)
  if (s != "DecJanFeb")
    paste0(s, "_", y)
  else if (m == 12)
    paste0(s, "_", y, "/", y + 1)
  else
    paste0(s, "_", y - 1, "/", y)
}

ws <- function(x) {
  x <- x[x > 0]
  if (length(x) > 0) {
    l <- 1 / (length(x) * x)
    l1 <- NA
    for (j in 1:length(x)) {
      l1 <- c(l1, rep(l[j], times = x[j]))
    }
    l1 <- l1[-1]
    return(l1)
  } else {
    return(NA)
  }
}

w_median <- function(x, w) {
  ind <- which(!is.na(x) & !is.na(w))
  if (length(ind) > 0)
    weighted.median(x[ind], w[ind])
  else
    NA
}

w_median_int <- function(x, w) {  
  x <- unlist(x)
  w <- unlist(w)
  x <- x[!is.na(x)]
  if (length(x) > 0) {
    return(weighted.median(x, w))
  } else {
    return(NA)
  }
}

w_sd <- function(x, w) {
  ind <- which(!is.na(x) & !is.na(w))
  if (length(ind) > 1)
    weighted.sd(x[ind], w[ind])
  else
    NA
}

w_sd_int <- function(x, w) {
  x <- unlist(x)
  w <- unlist(w)
  x <- x[!is.na(x)]
  if (length(x) > 1) {
    return(weighted.sd(x, w))
  } else {
    return(NA)
  }
}

seas <- function(m, y) {
  m <- as.numeric(m)
  y <- as.numeric(y)
  if (m %in% c(3, 4, 5)) {
    return(paste0("MarAprMay_", y))
  } 
  else if (m %in% c(6, 7, 8)) {
    return(paste0("JunJulAug_", y))
  }
  else if (m %in% c(9, 10, 11)) {
    return(paste0("SepOctNov_", y))
  }
  else if (m %in% c(1, 2)) {
    return(paste0("DecJanFeb_", y - 1, "/", y))
  }
  else if (m == 12) {
    return(paste0("DecJanFeb_", y, "/", y + 1))
  }
}

seas1 <- function(m) {  # Funtion for lev3 integrates, similar to sson1
  m <- as.numeric(m)
  if (m %in% c(3, 4, 5)) {
    return("MarAprMay")
  }
  else if (m %in% c(6, 7, 8)) {
    return("JunJulAug")
  }
  else if (m %in% c(9, 10, 11)) {
    return("SepOctNov")
  }
  else if (m %in% c(1, 2, 12)) {
    return("DecJanFeb")
  }
}

f1 <- function(x, i) {
  ifelse(!is.na(x), i, NA)
}

# Function to determine if a year is a leap year
is_leap_year <- function(year) {
  (year %% 4 == 0 & year %% 100 != 0) | (year %% 400 == 0)
}

# Function to find the last non-null value in a row
last_non_na <- function(row) {
  non_na_indices <- which(!is.na(row) & row != "")
  if (length(non_na_indices) == 0) {
    return(NA)
  } else {
    return(row[max(non_na_indices)])
  }
}

# Function to calculate the number of days in February
days_in_february <- function(year) {
  if ((year %% 4 == 0 && year %% 100 != 0) || (year %% 400 == 0)) {
    return(29)
  } else {
    return(28)
  }
}

ang <- function(x1, x2) {
  log(x1 / x2) / log(532 / 355)
}

ang_err <- function(x1, x2, e1, e2) {
  ((e1 / x1) + (e2 / x2)) / log(532 / 355)
}

# Functions to filter and merge data
filter_and_merge <- function(db, station, app) {
  db_filtered <- db[db$Station == station, ]
  db_merged <- merge(db_filtered, app, by = "Month", all = TRUE)
  return(db_merged)
}

filter_and_merge_seas <- function(db, station, app) {
  db_filtered <- db[db$Station == station, ]
  db_merged <- merge(db_filtered, app, by = "Season", all = TRUE)
  return(db_merged)
}

filter_and_merge_seas_y <- function(db, station, year, app) {
  filtered_db <- db[db$Station == station & substr(db$Season_Year, nchar(db$Season_Year) - 3, nchar(db$Season_Year)) == as.character(year), ]
  merged_db <- merge(filtered_db, app, by = "Season_Year", all = TRUE)
  return(merged_db)
}

filter_and_merge_y <- function(db, station, year) {
  app <- data.frame(Year = year, Null = NA)
  filtered_db <- db[db$Station == station & db$Year == year, ]
  merged_db <- merge(filtered_db, app, by = "Year", all = TRUE)
  return(merged_db)
}

# Function to get NetCDF variables with error handling
get_nc_var <- function(nc_file, var_name) {
  var_data <- try(ncvar_get(nc_file, var_name), silent = TRUE)
  
  if (!inherits(var_data, "try-error")) {
    if (length(class(var_data)) > 1) {
      var_data <- rowMeans(var_data, na.rm = TRUE)
    } 
  } else {
    var_data <- rep(NA, times = length(ncvar_get(nc_file, "altitude")))
  }
  return(var_data)
}

# Function to calculate weighted averages and errors
calculate_weighted_mean <- function(values, errors) {
  weights <- abs(values / errors)
  weights <- weights / sum(weights, na.rm = TRUE)
  weighted_mean <- weighted.mean(values, weights)
  mean_error <- mean(errors, na.rm = TRUE)
  return(list(weighted_mean = weighted_mean, mean_error = mean_error))
}

# Support function for safe histograms
safe_hist_counts <- function(x, breaks) {
  x <- as.numeric(x)
  x <- x[!is.na(x)]
  if (length(x) > 0 && is.numeric(x)) {
    return(hist(x, breaks = breaks, plot = FALSE)$counts)
  } else {
    return(rep(NA, length(breaks) - 1))  # or rep(0, ...) 
  }
}

# Support function for safe histograms for lidar ratio
safe_hist_lr_counts <- function(x, breaks) {
  x <- as.numeric(x)
  x <- x[!is.na(x)]
  x <- x[x >= 5 & x <= 160]  # only values between the bin range
  if (length(x) > 0 && is.numeric(x)) {
    return(hist(x, breaks = breaks, plot = FALSE)$counts)
  } else {
    return(rep(NA, length(breaks) - 1))  # or rep(0, ...) 
  }
}

# Support function for safe histograms for particle depolarizzation
safe_hist_pd_counts <- function(x, breaks) {
  x <- as.numeric(x)
  x <- x[!is.na(x)]
  x <- x[x >= 0 & x <= 0.55]  # only values between the bin range
  if (length(x) > 0 && is.numeric(x)) {
    return(hist(x, breaks = breaks, plot = FALSE)$counts)
  } else {
    return(rep(NA, length(breaks) - 1))  # or rep(0, ...) 
  }
}

# Function for correct variable assignment
assign_with_fallback <- function(target_df, source_df, weights_df, j, fun) {
  
  flag <- try(target_df[, j] <- mapply(fun, source_df[, j], weights_df[, j]), silent = TRUE)
  if (class(flag) == "try-error") {
    target_df[, j] <- fun(source_df[, j], weights_df[, j])
  }
  
  return(target_df)
}

# This script is developed to execute all project scripts in an orderly manner

result <- tryCatch({
  # Definition of release variable setting the start and end years
  release <- c(2000, 2021)                          
  release_suffix <- sprintf("_%02d%02d_", release[1] %% 100, release[2] %% 100)
  
  # Generation of a table with the list of level 2 files and their characteristics.
  source("./Lev2database.R")
  source("./Lev2database_layers.R")
  source("./system.R")
  
  # Construction of tables with profile data, divided by wavelength.
  source("./lev3pro_355.R")
  source("./lev3pro_532.R")
  source("./lev3pro_1064.R")
  
  # Generation of netCDF files containing profile data.
  flag_355 <- if ((exists("db_prof_nm_355") && nrow(db_prof_nm_355) > 0) && (exists("db_prof_ns_355") && nrow(db_prof_ns_355) > 0) && (exists("db_prof_season_355") && nrow(db_prof_season_355) > 0) && (exists("db_prof_year_355") && nrow(db_prof_year_355) > 0)) { TRUE } else { FALSE }
  flag_532 <- if ((exists("db_prof_nm_532") && nrow(db_prof_nm_532) > 0) && (exists("db_prof_ns_532") && nrow(db_prof_ns_532) > 0) && (exists("db_prof_season_532") && nrow(db_prof_season_532) > 0) && (exists("db_prof_year_532") && nrow(db_prof_year_532) > 0)) { TRUE } else { FALSE }
  flag_1064 <- if ((exists("db_prof_nm_1064") && nrow(db_prof_nm_1064) > 0) && (exists("db_prof_ns_1064") && nrow(db_prof_ns_1064) > 0) && (exists("db_prof_season_1064") && nrow(db_prof_season_1064) > 0) && (exists("db_prof_year_1064") && nrow(db_prof_year_1064) > 0)) { TRUE } else { FALSE }
  
  # Read CSV file
  stations <- read.csv("station.csv", stringsAsFactors = FALSE)
  
  # Get a list of folders in the ./New/ directory
  station_dirs <- list.dirs("./New", full.names = FALSE, recursive = FALSE)
  
  # Check which folders contain at least one file
  non_empty_dirs <- station_dirs[sapply(station_dirs, function(dir) {
    length(list.files(file.path("./New", dir))) > 0
  })]
  
  # Filter the dataframe to keep only stations with non-empty folders
  stations_filtered <- stations[stations$Code %in% non_empty_dirs, ]
  
  # Access station information
  loc <- stations_filtered$Location
  loc2 <- unique(lev2db$Station)
  latit <- stations_filtered$Latitude
  longit <- stations_filtered$Longitude
  station_alt <- stations_filtered$Altitude
  institution <- stations_filtered$Institution
  acronym <- stations_filtered$Acronym
  PI <- stations_filtered$PI
  PI_contact <- stations_filtered$PI_Contact
  
  syst <- as.matrix(syst_df)
  
  # Apply the function to each row of the matrix
  last_values <- apply(syst, 1, last_non_na)
  
  sist <- rep(NA, times = length(last_values))
  
  for (i in 1:length(last_values)) {
    if (row.names(syst_df)[i] == "azt") {
      sist[i] <- "NTUA Raman Lidar System ; EOLE"
    } else if (row.names(syst_df)[i] == "ipr") {
      sist[i] <- "JRC Ispra ; ADAM noew"
    } else if (row.names(syst_df)[i] == "lle") {
      sist[i] <- paste0(last_values[i], " ; Arielle")
    } else if (row.names(syst_df)[i] == "mel") {
      sist[i] <- paste0(last_values[i], " ; MSTL-2")
    } else if (row.names(syst_df)[i] == "nap") {
      sist[i] <- paste0(last_values[i], " ; ", syst[i, 18])
    } else if (row.names(syst_df)[i] == "sal") {
      sist[i] <- "PEARL ; MUSA"
    } else if (row.names(syst_df)[i] == "the") {
      sist[i] <- paste0(last_values[i], " ; ", syst[i, 19])
    } else {
      # Assign standard values
      sist[i] <- last_values[i]
    }
  }
  
  source("./lev3pro_nm_files.R")
  source("./lev3pro_ns_files.R")
  source("./lev3pro_s_files.R")  # update leap year check
  source("./lev3pro_y_files.R")
  
  # Construction of tables with integrated values, divided by type of source Level 2 files (e or b).
  source("./lev3int_b_e_pbl.R")
  
  # Reorganization of the integrated values, dividing them by wavelength.
  flag_355 <- if ((exists("db_int_e") && nrow(db_int_e) > 0 && '0355' %in% db_int_e$Wavelength) && (exists("db_int_b") && nrow(db_int_b) > 0 && '0355' %in% db_int_b$Wavelength)) { TRUE } else { FALSE }
  flag_532 <- if ((exists("db_int_e") && nrow(db_int_e) > 0 && '0532' %in% db_int_e$Wavelength) && (exists("db_int_b") && nrow(db_int_b) > 0 && '0532' %in% db_int_b$Wavelength)) { TRUE } else { FALSE } 
  flag_1064 <- if (exists("db_int_b") && nrow(db_int_b) > 0 && '1064' %in% db_int_b$Wavelength) { TRUE } else { FALSE }
  
  source("./lev3int_355.R")
  source("./lev3int_532.R")
  source("./lev3int_1064.R")
  
  # Generation of netCDF files containing profile data.
  flag_355 <- if ((exists("db_b_nm_355") && nrow(db_b_nm_355) > 0) && (exists("db_e_nm_355") && nrow(db_e_nm_355) > 0) && (exists("db_pd_nm_355") && nrow(db_pd_nm_355) > 0) && (exists("db_b_ns_355") && nrow(db_b_ns_355) > 0) && (exists("db_e_ns_355") && nrow(db_e_ns_355) > 0) && (exists("db_pd_ns_355") && nrow(db_pd_ns_355) > 0) && (exists("db_b_s_355") && nrow(db_b_s_355) > 0) && (exists("db_e_s_355") && nrow(db_e_s_355) > 0) && (exists("db_pd_s_355") && nrow(db_pd_s_355) > 0) && (exists("db_b_y_355") && nrow(db_b_y_355) > 0) && (exists("db_e_y_355") && nrow(db_e_y_355) > 0) && (exists("db_pd_y_355") && nrow(db_pd_y_355) > 0)) { TRUE } else { FALSE }
  flag_532 <- if ((exists("db_b_nm_532") && nrow(db_b_nm_532) > 0) && (exists("db_e_nm_532") && nrow(db_e_nm_532) > 0) && (exists("db_pd_nm_532") && nrow(db_pd_nm_532) > 0) && (exists("db_b_ns_532") && nrow(db_b_ns_532) > 0) && (exists("db_e_ns_532") && nrow(db_e_ns_532) > 0) && (exists("db_pd_ns_532") && nrow(db_pd_ns_532) > 0) && (exists("db_b_s_532") && nrow(db_b_s_532) > 0) && (exists("db_e_s_532") && nrow(db_e_s_532) > 0) && (exists("db_pd_s_532") && nrow(db_pd_s_532) > 0) && (exists("db_b_y_532") && nrow(db_b_y_532) > 0) && (exists("db_e_y_532") && nrow(db_e_y_532) > 0) && (exists("db_pd_y_532") && nrow(db_pd_y_532) > 0)) { TRUE } else { FALSE }
  flag_1064 <- if ((exists("db_b_nm_1064") && nrow(db_b_nm_1064) > 0) && (exists("db_pd_nm_1064") && nrow(db_pd_nm_1064) > 0) && (exists("db_b_ns_1064") && nrow(db_b_ns_1064) > 0) && (exists("db_pd_ns_1064") && nrow(db_pd_ns_1064) > 0) && (exists("db_b_s_1064") && nrow(db_b_s_1064) > 0) && (exists("db_pd_s_1064") && nrow(db_pd_s_1064) > 0) && (exists("db_b_y_1064") && nrow(db_b_y_1064) > 0) && (exists("db_pd_y_1064") && nrow(db_pd_y_1064) > 0)) { TRUE } else { FALSE }
  flag_ang_nm <- if (exists("db_ang_nm_tot") && nrow(db_ang_nm_tot) > 0) { TRUE } else { FALSE } 
  flag_pbl_nm <- if (exists("db_pbl_nm_tot") && nrow(db_pbl_nm_tot) > 0) { TRUE } else { FALSE } 
  flag_ang_ns <- if (exists("db_ang_ns_tot") && nrow(db_ang_ns_tot) > 0) { TRUE } else { FALSE } 
  flag_pbl_ns <- if (exists("db_pbl_ns_tot") && nrow(db_pbl_ns_tot) > 0) { TRUE } else { FALSE }
  flag_ang_s <- if (exists("db_ang_s_tot") && nrow(db_ang_s_tot) > 0) { TRUE } else { FALSE } 
  flag_pbl_s <- if (exists("db_pbl_s_tot") && nrow(db_pbl_s_tot) > 0) { TRUE } else { FALSE }
  flag_ang_y <- if (exists("db_ang_y_tot") && nrow(db_ang_y_tot) > 0) { TRUE } else { FALSE } 
  flag_pbl_y <- if (exists("db_pbl_y_tot") && nrow(db_pbl_y_tot) > 0) { TRUE } else { FALSE }
    
  source("./lev3int_nm_files.R")
  source("./lev3int_ns_files.R")
  source("./lev3int_s_files.R")   # update leap year check
  source("./lev3int_y_files.R")

  # Generation of the table with the layer values.
  source("./lev3int_layers.R")
  
  # Generation of netCDF files containing the layer values.
  db_int_layers$Extinction <- db_int_layers$Extinction * 10 ^ 3
  db_int_layers$Error_Extinction <- db_int_layers$Error_Extinction * 10 ^ 3
  db_int_layers$IntBs <- db_int_layers$IntBs * 10 ^ 3
  db_int_layers$Error_IntBs <- db_int_layers$Error_IntBs * 10 ^ 3
  db_int_layers$Backscatter <- db_int_layers$Backscatter * 10 ^ 6
  db_int_layers$Error_Backscatter <- db_int_layers$Error_Backscatter * 10 ^ 6
  
  altitude_breaks <- seq(from = 0, to = 20000, by = 1000)
  altitude_intervals_data <- altitude_breaks[1:20]
  
  # # Old version with dynamic bin calculation
  # vvv <- db_int_layers$LidarRatio
  # vvv <- vvv[!is.na(vvv)]
  # if (length(vvv) > 0) {
  #   nnn <- as.numeric(quantile(vvv, probs = c(0.1, 0.9)))
  #   lr_breaks <- c(min(vvv), seq(from = nnn[1], to = nnn[2], by = ((nnn[2] - nnn[1]) / 18)), max(vvv))
  #   lr_intervals_data <- lr_breaks[1:20]
  # } else {
  #   lr_breaks <- seq(0, 0, length.out = 21)
  #   lr_intervals_data <- lr_breaks[1:20]
  # }
  # 
  # vvv <- db_int_layers$ParticleDep
  # vvv <- vvv[!is.na(vvv)]
  # if (length(vvv) > 0) {
  #   nnn <- as.numeric(quantile(vvv, probs = c(0.1, 0.9)))
  #   pd_breaks <- c(min(vvv), seq(from = nnn[1], to = nnn[2], by = ((nnn[2] - nnn[1]) / 18)), max(vvv))
  #   pd_intervals_data <- pd_breaks[1:20]
  # } else {
  #   pd_breaks <- seq(0, 0, length.out = 21)  
  #   pd_intervals_data <- pd_breaks[1:20]
  # }
  
  vvv <- db_int_layers$LidarRatio
  vvv <- vvv[!is.na(vvv)]
  lr_breaks <- seq(5, 160, by = 5)
  if (length(vvv) > 0) {
    lr_intervals_data <- lr_breaks
  } else {
    lr_breaks <- seq(0, 0, length.out = length(lr_breaks))
    lr_intervals_data <- rep(0, length(lr_breaks))
  }
  
  vvv <- db_int_layers$ParticleDep
  vvv <- vvv[!is.na(vvv)]
  pd_breaks <- seq(0, 0.55, by = 0.025)
  if (length(vvv) > 0) {
    pd_intervals_data <- pd_breaks
  } else {
    pd_breaks <- seq(0, 0, length.out = length(pd_breaks)) 
    pd_intervals_data <- rep(0, length(pd_breaks))
  }
  
  vvv <- db_int_layers$AOD
  vvv <- vvv[!is.na(vvv)]
  if (length(vvv) > 0) {
    nnn <- as.numeric(quantile(vvv, probs = c(0.1, 0.9)))
    aod_breaks <- c(min(vvv), seq(from = nnn[1], to = nnn[2], by = ((nnn[2] - nnn[1]) / 18)), max(vvv))
    aod_intervals_data <- aod_breaks[1:20]
  } else {
    aod_breaks <- seq(0, 0, length.out = 21)
    aod_intervals_data <- aod_breaks[1:20]
  }
  
  vvv <- db_int_layers$Extinction
  vvv <- vvv[!is.na(vvv)]
  if (length(vvv) > 0) {
    nnn <- as.numeric(quantile(vvv, probs = c(0.1, 0.9)))
    extinction_breaks <- c(min(vvv), seq(from = nnn[1], to = nnn[2], by = ((nnn[2] - nnn[1]) / 18)), max(vvv))
    extinction_intervals_data <- extinction_breaks[1:20]
  } else {
    extinction_breaks <- seq(0, 0, length.out = 21)
    extinction_intervals_data <- extinction_breaks[1:20]
  }
  
  vvv <- db_int_layers$IntBs
  vvv <- vvv[!is.na(vvv)]
  if (length(vvv) > 0) {
    nnn <- as.numeric(quantile(vvv, probs = c(0.1, 0.9)))
    ib_breaks <- c(min(vvv), seq(from = nnn[1], to = nnn[2], by = ((nnn[2] - nnn[1]) / 18)), max(vvv))
    ib_intervals_data <- ib_breaks[1:20]
  } else {
    ib_breaks <- seq(0, 0, length.out = 21)
    ib_intervals_data <- ib_breaks[1:20]
  }
  
  vvv <- db_int_layers$Backscatter
  vvv <- vvv[!is.na(vvv)]
  if (length(vvv) > 0) {
    nnn <- as.numeric(quantile(vvv, probs = c(0.1, 0.9)))
    backscatter_breaks <- c(min(vvv), seq(from = nnn[1], to = nnn[2], by = ((nnn[2] - nnn[1]) / 18)), max(vvv))
    backscatter_intervals_data <- backscatter_breaks[1:20]
  } else {
    backscatter_breaks <- seq(0, 0, length.out = 21)
    backscatter_intervals_data <- backscatter_breaks[1:20]
  }
  
  source("./lev3layers_nm_files.R")  
  source("./lev3layers_ns_files.R")  
  source("./lev3layers_s_files.R")  
  source("./lev3layers_y_files.R")  
  log(-1)
}, error = function(e) {
  cat("Errore catturato: ", e$message, "\n")
  NA
}, finally = {
  closeAllConnections()
  cat("Esecuzione del blocco finally\n")
})

# , warning = function(w) {
#  cat("Avvertimento: ", w$message, "\n")
#  invokeRestart("muffleWarning")
# }
