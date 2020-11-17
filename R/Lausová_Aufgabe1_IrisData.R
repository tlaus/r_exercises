## Aufgabe 1 - Iris Data
## Tereza Lausová
## 14.9.2020

## Load libraries
library(tidyverse)
library(here)
library(patchwork)

## Load and explore data
data(iris)
# plot(iris)

# Sepal plot
sepal_c <- ggplot(iris) +
  theme_minimal() +
  geom_rect(aes(xmax = Sepal.Width+Sepal.Width/20, ymin = Sepal.Length-Sepal.Length/20, xmin=Sepal.Width-Sepal.Width/20, ymax=Sepal.Length+Sepal.Length/20, fill = Species), color = "black", alpha = 0.5) +
  xlim(0,9) +
  ylim(0,9) +
  theme(legend.position = c(0.75, 0.5), 
        legend.background = element_rect(colour = NA, fill = "white")) +
  coord_fixed() +
  labs(title = "Sepal size in the iris dataset in proportion")

sepal_b <- ggplot(iris) +
  theme_minimal() +
  geom_density(aes(x = Sepal.Width, fill = Species), alpha = 0.5) +
  xlim(0,9) +
  guides(fill = "none")

sepal_l <- ggplot(iris) +
  theme_minimal() +
  geom_density(aes(y = Sepal.Length, fill = Species), alpha = 0.5) +
  ylim(0,9) +
  guides(fill = "none")


sepal_plot <- sepal_l + sepal_c + plot_spacer() + sepal_b + plot_layout(widths = c(1,6), heights = c(6,1))

pdf(file = here("plots/Lausová_Aufgabe1_IrisData1.pdf"))
sepal_plot
dev.off()

# Petal plot
petal_c <- ggplot(iris) +
  theme_minimal() +
  # geom_ellipse(aes(x0 = Petal.Width, y0 = Petal.Length, a = Petal.Width/30, b = Petal.Length/30, angle = 0, fill = Species), alpha = 0.5) +
  # geom_linerange(aes(x = Petal.Width, size = Petal.Width, ymin = Petal.Length-(Petal.Length/20), ymax = Petal.Length+(Petal.Length/20), color = Species)) +
  geom_rect(aes(xmax = Petal.Width+Petal.Width/20, ymin = Petal.Length-Petal.Length/20, xmin=Petal.Width-Petal.Width/20, ymax=Petal.Length+Petal.Length/20, fill = Species), color = "black", alpha = 0.5) +
  # coord_fixed() +
  theme(legend.position = c(0.75, 0.5), 
        legend.background = element_rect(colour = NA, fill = "white")) +
  xlim(0,8) +
  ylim(0,8) +
  labs(title = "Petal size in the iris dataset in proportion")

petal_b <- ggplot(iris) +
  theme_minimal() +
  geom_density(aes(x = Petal.Width, fill = Species), alpha = 0.5) +
  xlim(0,8) +
  # coord_fixed(1/8) +
  guides(fill = "none")

petal_l <- ggplot(iris) +
  theme_minimal() +
  geom_density(aes(y = Petal.Length, fill = Species), alpha = 0.5) +
  ylim(0,8) +
  # coord_fixed(8/1) +
  guides(fill = "none")


petal_plot <- petal_l + petal_c + plot_spacer() + petal_b + plot_layout(widths = c(1,6), heights = c(6,1))


pdf(file = here("plots/Lausová_Aufgabe1_IrisData2.pdf"))
petal_plot
dev.off()
