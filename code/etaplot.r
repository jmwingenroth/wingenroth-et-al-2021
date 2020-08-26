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

trapdata <- lapply(trapdata, function(x) select(.data = x, trap = ID, m_trap = Sediment)) %>%
  bind_rows(.id = "run") %>%
  filter(!is.na(trap), !is.na(m_trap))



## K_TOT MODELS

k_fits <- lmList(log(mvc) ~ t | run, data = pumpdata)

lm_results <- as.tibble(summary(k_fits)$coefficients[,,'t']) %>%
  bind_cols(id = names(k_fits)) %>%
  left_join(metadata, by = c("id" = "Date")) %>%
  transmute(id, 
         k_t = -Estimate, 
         dk_t = `Std. Error`, 
         density = Density, 
         freq = Velocity, 
         biofilm = Biofouled)



## K_C AND ETA CALCULATIONS

trapmeans <- trapdata %>%
  group_by(run) %>%
  filter(m_trap > 1) %>%
  summarise(m_trap = mean(m_trap))

df <- trapmeans %>%
  left_join(lm_results, by = c("run" = "id"))

# backgrounds: old 30 Hz one was soon after I started in the lab and more prone to error

df %>%
  filter(density == 0)

# 041719 was also an outlier, 232 density treatments were in error, and biofouled runs not of interest

df <- df %>%
  filter(!run %in% c("101918", "041719"), density != 232, biofilm==0)
  
df <- df %>%
  mutate(m_s = m_trap/1000*1.95*.6/(pi*.0127^2),
         k_s = m_s*k_t/200/(1-exp(-k_t*6000)),
         k_bg = ifelse(density == 0, k_t-k_s, NA)) %>%
  group_by(freq) %>%
  mutate(k_bg = mean(k_bg, na.rm = TRUE),
         k_c = k_t - k_s - k_bg)

df %>%
  ggplot(aes(x = freq, y = k_t, color = factor(density))) +
  geom_line()

## PLOTS