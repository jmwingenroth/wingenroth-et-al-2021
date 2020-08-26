library(tidyverse)

data <- read_csv("../data/candace results.csv")

data <- data %>%
  mutate(velocity = Velocity/500,
         eta = k_c/velocity/.003175/Density*2.43/(1.95*.4*.6))

etafig <- data %>%
  filter(Biofouled == 0, `Frontal Area`>0) %>%
  ggplot(aes(x = Re_c, color = factor(sprintf("%.2f", `Frontal Area`)), y = eta)) +
  geom_line(size = 1.2) +
  geom_point(size = 3.3, color = "black") +
  geom_point(size = 1.5) +
  theme_bw(base_size = 16) +
  labs(y = "Effective capture efficiency (%)", 
       x = "Collector Reynolds number", 
       color = expression(Frontal~Area~(m^-1))) +
  scale_x_continuous(breaks = c(68.1,133.6,199.0)) +
  theme(legend.position = c(.753,.78), 
        legend.background = element_rect(fill = "white", colour = "black"), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "black", size = .3),
        panel.background = element_rect(colour = "black", size = 2))

etafig

ggsave("../pics/etafig.png", width = 7, height = 5, units = "in", dpi = 300)
  
