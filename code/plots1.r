library(tidyverse)

data <- read_csv("../data/candace results.csv")

data %>%
  mutate(velocity = Velocity/500,
         eta = k_c/velocity/.003175/Density)

etafig <- data %>%
  filter(Biofouled == 0, `Frontal Area`>0) %>%
  ggplot(aes(x = Re_c, color = factor(`Frontal Area`), y = `ECE %`/100*2.43/1.95/.4/.6)) +
  geom_line(size = 1.2) +
  geom_point(size = 3.3, color = "black") +
  geom_point(size = 1.5) +
  theme_bw(base_size = 16) +
  scale_y_continuous(limits = c(0,.2)) +
  scale_x_continuous(limits = c(50,200)) +
  labs(y = expression(eta*minute), 
       x = expression(Re[c]), 
       color = expression(Frontal~Area~(m^2/m^3))) +
  theme(legend.position = c(.753,.78), 
        legend.background = element_rect(fill = "white", colour = "black"), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "black", size = .3),
        panel.background = element_rect(colour = "black", size = 2))

etafig

ggsave("../pics/etafig.png", width = 7, height = 5, units = "in", dpi = 300)
  
