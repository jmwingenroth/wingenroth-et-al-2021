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

pumpdata <- pumpdata %>%
  rename(mvc = `mass volume concentration (ppm)`, t = `Time (s)`)



## K_TOT MODELS

k_fits <- lmList(log(mvc) ~ t | run, data = pumpdata)

lm_results <- as.tibble(summary(k_fits)$coefficients[,,'t']) %>%
  bind_cols(id = names(k_fits)) %>%
  left_join(metadata, by = c("id" = "Date")) %>%
  select(id, 
         k_t = Estimate, 
         dk_t = `Std. Error`, 
         density = Density, 
         freq = Velocity, 
         biofilm = Biofouled)



## K_C AND ETA CALCULATIONS

df <- tibble(id = names(trapdata),
       m_s = unlist(lapply(trapdata, function(x) mean(as.numeric(x$Sediment),na.rm = T)))) %>%
  left_join(lm_results, by = "id")



## PLOTS