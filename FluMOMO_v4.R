# Install required packages (if not already installed)
# install.packages("readstata13")


# Clear enviromnet i.e. delete all data and variabels
rm(list = ls())
# Clear console
cat("\014")

# Set default work directory
setwd("C:/users/public")
# Use RAM disk
if (dir.exists("W:/")) setwd("W:/")

# Country
country <- "Denmark"

# Work directory - directory where the Stata programs are placed
#wdir <- "."
wdir="H:/SFSD/INFEPI/Projekter/AKTIVE/MOMO/FluMOMO version 4/FluMOMO v4 R"

# Period: Start and end (both included)
start_year <- 2012
start_week <- 27
end_year <- 2018
end_week <- 5

# Deaths data from A-MOMO
# 1 = A-MOMO complete file, renamed to A-MOMO data.txt
# 0 = you provide a Stata data file: deaths.dta, containing the variable: agegrp, year, week, deaths
A_MOMO <- 1

# Population data (TRUE/FALSE)
population <- TRUE

# Restrict IA to only positive (TRUE/FALSE)
IArest <- FALSE

# Number of IA lags
# 0 = no lag, 1 = one week lag, ...(max=9)
IAlags <- 2

# Number of OF lags
# 0 = no lag, 1 = one week lag, ...(max=9)
ETlags <- 2


### Direcory setup #######################################################################

### Create general output dir
if (!dir.exists(paste0(wdir,"/FluMOMO_",end_year,"w",formatC(end_week, width=2, flag="0"))))
  { dir.create(paste0(wdir,"/FluMOMO_",end_year,"w",formatC(end_week, width=2, flag="0"))) }

### Copy data directory - directory where input data are stored
indir <- paste0(wdir,"/FluMOMO_",end_year,"w",formatC(end_week, width=2, flag="0"),"/data")
if (!dir.exists(indir)) { dir.create(indir) }
file.copy(from = list.files(paste0(wdir,"/data"), all.files = TRUE, full.names = TRUE, no.. = TRUE),
          to = indir, overwrite = TRUE, recursive = FALSE, copy.mode = TRUE)
#file.remove(list.files(paste0(wdir,"/data"), all.files = TRUE, full.names = TRUE, no.. = TRUE))

### Create output directory - directory where output are created
outdir <- paste0(wdir,"/FluMOMO_",end_year,"w",formatC(end_week, width=2, flag="0"),"/output")
if (!dir.exists(outdir)) { dir.create(outdir) }
#########################################################################################

source(paste0(wdir,"/Estimation_v41.R"), echo = FALSE)
# Output: csv-files
source(paste0(wdir,"/Output_csv_v4.R"))
# Output: graphs IA and ET
source(paste0(wdir,"/Output_IA_ET_v4.R"))
# Output: graphs over calendar time
source(paste0(wdir,"/Output_calendar_v4.R"))
# Output: graphs cumulated IA
source(paste0(wdir,"/Output_cumulated_v4.R"))
