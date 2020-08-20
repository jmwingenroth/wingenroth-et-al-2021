## READ DATA

library(tidyverse)
library(lme4)

files <- list.files("../data/k_runs/")
pump <- str_detect(files, "pump")
trap <- str_detect(files, "trap")

pumpdata <- lapply(paste0("../data/k_runs/",files[pump]), read_csv)
trapdata <- lapply(paste0("../data/k_runs/",files[trap]), read_csv)

metadata <- read_csv("../data/flumeexperimentsettings.csv")

## TIDY DATA

names(pumpdata) <- str_sub(files[pump],,6)
names(trapdata) <- str_sub(files[trap],,6)

pumpdata <- bind_rows(pumpdata, .id = "run")

pumpdata