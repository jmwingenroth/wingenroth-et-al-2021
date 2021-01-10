library(tidyverse)

pumpfiles <- list.files("../data/peristaltic pumps/")
trapfiles <- list.files("../data/sediment traps/")

pumpdate <- str_sub(pumpfiles,0,6)
trapdate <- str_sub(trapfiles,0,6)



#PERISTALTIC PUMPS

pumpdata <- lapply(pumpfiles, function(x) read_csv(paste0("../data/peristaltic pumps/",x)))

names(pumpdata) <- pumpdate

x <- pumpdata

tidypump <- lapply(seq_along(x), function(i) {
  select(x[[i]], 
         loc = Location, 
         ht = Height,
         t = `time series`,
         mvc = contains("(ppm)")) %>%
    filter(t < 21) %>% #filter a few timepoints outside the normal window 
    mutate(t = (t-min(t)+1)*300, #convert from timestep to seconds
           mvc = as.numeric(mvc), 
           date = as.numeric(names(x)[[i]])) %>%
    filter(mvc<80, mvc>1) #outliers were removed based on the residual graph
}
)

pump <- bind_rows(tidypump)

pump <- pump %>%
  
  #first 3 runs had starting sediment mass of 100g rather than 200g
  
  filter(date > 181005) %>%
  
  #we can't use runs without sediment mass
  
  filter(as.character(date) %in% trapdate)


metadata <- read_csv("../data/run_metadata.csv")

pump <- left_join(pump, metadata, by = "date")

pumpfinal <- pump %>%
  filter(date != 190417,
         dowel_density != "0232")



# SEDIMENT TRAPS

trapdata <- lapply(trapfiles, function(x) read_csv(paste0("../data/sediment traps/",x)))

names(trapdata) <- trapdate

x <- trapdata

tidytrap <- lapply(seq_along(x), function(i) {
  select(x[[i]], station = 1, pre = contains("pre"), post = contains("post"), sed = contains("sed")) %>%
    mutate(date = names(x)[i]) %>%
    mutate_at(vars(pre,post,sed,date), as.numeric)
}
)

trap <- bind_rows(tidytrap) %>%
  filter(!is.na(sed))


trapfinal <- left_join(trap, metadata, by = "date") %>%
  filter(date %in% pumpfinal$date)



#ERROR PROPOGATION FUNCTION

mutate_with_error = function(.data, f) {
  exprs = list(
    # expression to compute new variable values
    deparse(f[[3]]),
    
    # expression to compute new variable errors
    sapply(all.vars(f[[3]]), function(v) {
      dfdp = deparse(D(f[[3]], v))
      sprintf('(d%s*(%s))^2', v, dfdp)
    }) %>%
      paste(collapse='+') %>%
      sprintf('sqrt(%s)', .)
  )
  names(exprs) = c(
    deparse(f[[2]]),
    sprintf('d%s', deparse(f[[2]]))
  )
  
  .data %>%
    # the standard evaluation alternative of mutate()
    mutate_(.dots=exprs)
}

# K_TOT

library(lme4)

fits <- lmList(data = pumpfinal, log(mvc)~t | date)

summary(fits)

cbind(pumpfinal, resid = residuals(fits)) %>%
  ggplot(aes(x = t, y = resid)) +
  geom_point() +
  geom_smooth()

fitdata <- coefficients(fits) %>%
  mutate(k_t = -t) %>%
  cbind(dk_t = summary(fits)[[4]][,,'t'][,2], date = as.numeric(row.names(coefficients(fits)))) %>%
  left_join(metadata, by = "date") %>%
  mutate(ad = as.numeric(dowel_density)*.003175^2*2/1.95, Re = pump_freq/500*.003175/9.509e-7) #DOWEL DENSITY CORRECTION INCLUDED HERE

fitdata %>%
  filter(pump_freq == 30, dowel_density!="0000") %>%
  ggplot(aes(x = growth_days, y = k_t, ymin = k_t - 1.96*dk_t, ymax = k_t + 1.96*dk_t, color = as.character(round(ad,5)))) +
  geom_line(lty = "dotted", size = 2, alpha = .3) + 
  geom_point() +
  geom_errorbar() +  
  scale_color_manual(values = c("blue","green","red")) +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  labs(y = expression(paste(italic(k[t])," (",s^{-1},")")), x = "Biofilm growth (days)", color = expression(italic(ad))) +
  theme_minimal()

# K_S

final <- trapfinal %>%
  group_by(date) %>%
  summarise(m_trap = mean(sed)/1000,
            dm_trap = sd(sed)/1000) %>% #average sediment in trap converted from mg to g
  mutate_with_error(m_s ~ m_trap*1.95*.6/(3.14159*.0127^2)) %>% #*A_TS/a_trap
  left_join(fitdata, by = "date") %>%
  mutate_with_error(c_s ~ m_s/200/(1-exp(-k_t*6000))) %>% #intermediate step to squash error 
  mutate_with_error(k_s ~ c_s*k_t)

final %>%
  filter(pump_freq == 30, dowel_density!="0000") %>%
  ggplot(aes(x = growth_days, y = k_s, ymin = k_s - 1.96*dk_s, ymax = k_s + 1.96*dk_s, color = as.character(round(ad,5)))) +
  geom_line(lty = "dotted", size = 2, alpha = .3) + 
  geom_point() +
  geom_errorbar() +  
  scale_color_manual(values = c("blue","green","red")) +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  labs(y = expression(paste(italic(k[s])," (",s^{-1},")")), x = "Biofilm growth (days)", color = expression(italic(ad))) +
  theme_minimal()

# K_BG

bgvals <- final %>%
  filter(dowel_density=="0000", date != 190802) %>% #took out control to avoid join error
  mutate_with_error(k_bg ~ k_t - k_s) %>%
  select(pump_freq, k_bg, dk_bg) %>%
  arrange(pump_freq)

bgvals

final <- left_join(final, bgvals, by = "pump_freq")

# K_C

final <- final %>%
  mutate_with_error(k_c ~ k_t - k_s - k_bg)

final %>%
  filter(pump_freq == 30, dowel_density!="0000") %>%
  ggplot(aes(x = growth_days, y = k_c, ymin = k_c - 1.96*dk_c, ymax = k_c + 1.96*dk_c, color = as.character(round(ad,4)))) +
  geom_line(lty = "dotted", size = 2, alpha = .3) + 
  geom_point() +
  geom_errorbar() +  
  scale_color_manual(values = c("blue","green","red")) +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  labs(y = expression(paste(italic(k[c])," (",s^{-1},")")), x = "Biofilm growth (days)", color = expression(italic(ad))) +
  theme_minimal()

# ETA


#temp in flume measured at 22.2C with a calibrated thermometer
#plugged into https://www.engineeringtoolbox.com/water-dynamic-kinematic-viscosity-d_596.html
visc = 9.509e-7 #kinematic viscosity, m2/s

d = .003175 #dowel diameter

final <- final %>%
  mutate(frontal_area = as.numeric(dowel_density)*.003175, 
         u = as.numeric(pump_freq)/500, #velocity 
         Re = u*.003175/visc, #Reynolds #
         k_t = -t, 
         k_c, k_s, k_bg) %>%
  mutate(eta = k_c/u/frontal_area,
         deta = dk_c/u/frontal_area) %>%
  mutate(eta = eta * 2.43/(1.95*.4*.6),
         deta = deta * 2.43/(1.95*.4*.6)) # don't forget to correct for time out of test section!

biofilmplot <- final %>%
  filter(pump_freq == 30, dowel_density!="0000",growth_days!=30) %>% #the 30-day growth period was not controlled correctly, invalidating that run's data
  ggplot(aes(x = growth_days, y = eta, ymin = eta - 1.96*deta, ymax = eta + 1.96*deta, color = factor(ad, labels = c("0.22%","0.64%","1.17%")))) +
  geom_line() + 
  geom_point() +
  geom_ribbon(aes(fill = factor(ad)), alpha = .2, color = NA, show.legend = FALSE) +
  scale_y_continuous(labels = scales::percent) +
  labs(y = "Effective Capture Efficiency", x = "Biofilm Growth (Days)", color = "Collector\nSolid\nVolume\nFraction") +
  theme_bw() +
  theme(legend.background = element_rect(colour = "black", size = .2))

ggsave("../pics/biofilm.png",biofilmplot,dpi = 600, width = 7, height = 5)

