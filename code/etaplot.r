## READ DATA

library(tidyverse)

files <- list.files("../data/k_runs/")
pump <- str_detect(files, "pump")
trap <- str_detect(files, "trap")

pumpdata <- lapply(paste0("../data/k_runs/",files[pump]), read_csv)
trapdata <- lapply(paste0("../data/k_runs/",files[trap]), read_csv)

names(pumpdata) <- str_sub(files[pump],,6)
names(trapdata) <- str_sub(files[trap],,6)

metadata <- read_csv("../data/flumeexperimentsettings.csv")

## K_TOT MODELS

colnames(pumpdata) <- lapply(pumpdata, function(x) colnames(x) <- tolower(colnames(x)))

lapply(pumpdata, function(x) rename(.data = x, t = "Time (s)", mvc = "mass volume concentration (ppm)"))
