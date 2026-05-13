library(tidyverse)
library(lme4)
library(car)
library(emmeans)
library(ggprism)
library(forcats)
library(dplyr)
library(ggsignif)

actin_live <- read.csv("actin tp level full.csv")
scene1 <- filter(actin_live, scene %in% c("scene 1"))
scene2 <- filter(actin_live, scene %in% c("scene 2"))
scene3 <- filter(actin_live, scene %in% c("scene 3"))
scene4 <- filter(actin_live, scene %in% c("scene 4"))
scene5 <- filter(actin_live, scene %in% c("scene 5"))
live18yo <- filter(data_trimmed, age %in% "18yo")
live29yo <- filter(data_trimmed, age %in% "29yo")
live32yo <- filter(data_trimmed, age %in% "32yo")
live59yo <- filter(data_trimmed, age %in% "59yo")

data_trimmed <- actin_live %>%
  filter(time_min <= 70)

ggplot(mutate(, scene = fct_relevel(age, "18yo", "29yo", "32yo", "59yo")), aes(time_min, sum_length_2d_pct_change, color = age)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_point(size = 0.5, position = position_jitterdodge()) +
  labs(x = "Time (min)", y = "Total 2D filament length % change", color = "Age") +
  theme_classic() +
  theme(legend.position = "top") +
  scale_color_brewer(palette = "Set2", labels = c("18 y/o", "29 y/o", "32 y/o", "59 y/o")) +
  guides(y = guide_axis(minor.ticks = TRUE), x = guide_axis(minor.ticks = TRUE)) +
  scale_y_continuous(minor_breaks = scales::breaks_width(10)) +
  scale_x_continuous(minor_breaks = scales::breaks_width(10)) 
ggsave("sum2dpctchangeScene1.tiff", device = "tiff", dpi = 300, width = 177, height = 150, units = "mm")

ggplot(mutate(, scene = fct_relevel(age, "18yo", "29yo", "32yo", "59yo")), aes(time_min, sum_length_2d_pct_change, color = age)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_point(size = 0.5, position = position_jitterdodge()) +
  labs(x = "Time (min)", y = "Total 2D filament length % change", color = "Age") +
  theme_classic() +
  theme(legend.position = "top") +
  scale_color_brewer(palette = "Set2", labels = c("18 y/o", "29 y/o", "32 y/o", "59 y/o")) +
  guides(y = guide_axis(minor.ticks = TRUE), x = guide_axis(minor.ticks = TRUE)) +
  scale_y_continuous(minor_breaks = scales::breaks_width(10)) +
  scale_x_continuous(minor_breaks = scales::breaks_width(10)) 
ggsave("sum2dpctchangeScene1.tiff", device = "tiff", dpi = 300, width = 177, height = 150, units = "mm")

ggplot(mutate(scene1, scene = fct_relevel(age, "18yo", "29yo", "32yo", "59yo")), aes(time_min, sum_length_2d_pct_change, color = age)) + 
  geom_line() +
  geom_point(size = 0.5, position = position_jitterdodge()) +
  labs(x = "Time (min)", y = "Total 2D filament length % change", color = "Age") +
  theme_classic() +
  theme(legend.position = "top") +
  scale_color_brewer(palette = "Set2", labels = c("18 y/o", "29 y/o", "32 y/o", "59 y/o")) +
  guides(y = guide_axis(minor.ticks = TRUE), x = guide_axis(minor.ticks = TRUE)) +
  scale_y_continuous(minor_breaks = scales::breaks_width(10)) +
  scale_x_continuous(minor_breaks = scales::breaks_width(10)) 
ggsave("sum2dpctchangeScene1.tiff", device = "tiff", dpi = 300, width = 177, height = 150, units = "mm")

ggplot(mutate(scene2, scene = fct_relevel(age, "18yo", "29yo", "32yo", "59yo")), aes(time_min, sum_length_2d_pct_change, color = age)) + 
  geom_line() +
  geom_point(size = 0.5, position = position_jitterdodge()) +
  labs(x = "Time (min)", y = "Total 2D filament length % change", color = "Age") +
  theme_classic() +
  theme(legend.position = "top") +
  scale_color_brewer(palette = "Set2", labels = c("18 y/o", "29 y/o", "32 y/o", "59 y/o")) +
  guides(y = guide_axis(minor.ticks = TRUE), x = guide_axis(minor.ticks = TRUE)) +
  scale_y_continuous(minor_breaks = scales::breaks_width(10)) +
  scale_x_continuous(minor_breaks = scales::breaks_width(10)) 
ggsave("sum2dpctchangeScene2.tiff", device = "tiff", dpi = 300, width = 177, height = 150, units = "mm")

ggplot(mutate(scene3, scene = fct_relevel(age, "18yo", "29yo", "32yo", "59yo")), aes(time_min, sum_length_2d_pct_change, color = age)) + 
  geom_line() +
  geom_point(size = 0.5, position = position_jitterdodge()) +
  labs(x = "Time (min)", y = "Total 2D filament length % change", color = "Age") +
  theme_classic() +
  theme(legend.position = "top") +
  scale_color_brewer(palette = "Set2", labels = c("18 y/o", "29 y/o", "32 y/o", "59 y/o")) +
  guides(y = guide_axis(minor.ticks = TRUE), x = guide_axis(minor.ticks = TRUE)) +
  scale_y_continuous(minor_breaks = scales::breaks_width(10)) +
  scale_x_continuous(minor_breaks = scales::breaks_width(10)) 
ggsave("sum2dpctchangeScene3.tiff", device = "tiff", dpi = 300, width = 177, height = 150, units = "mm")

ggplot(mutate(scene4, scene = fct_relevel(age, "18yo", "29yo", "59yo")), aes(time_min, sum_length_2d_pct_change, color = age)) + 
  geom_line() +
  geom_point(size = 0.5, position = position_jitterdodge()) +
  labs(x = "Time (min)", y = "Total 2D filament length % change", color = "Age") +
  theme_classic() +
  theme(legend.position = "top") +
  scale_color_brewer(palette = "Set2", labels = c("18 y/o", "29 y/o", "59 y/o")) +
  guides(y = guide_axis(minor.ticks = TRUE), x = guide_axis(minor.ticks = TRUE)) +
  scale_y_continuous(minor_breaks = scales::breaks_width(10)) +
  scale_x_continuous(minor_breaks = scales::breaks_width(10)) 
ggsave("sum2dpctchangeScene4.tiff", device = "tiff", dpi = 300, width = 177, height = 150, units = "mm")

ggplot(mutate(scene5, scene = fct_relevel(age, "18yo", "29yo", "32yo", "59yo")), aes(time_min, sum_length_2d_pct_change, color = age)) + 
  geom_line() +
  geom_point(size = 0.5, position = position_jitterdodge()) +
  labs(x = "Time (min)", y = "Total 2D filament length % change", color = "Age") +
  theme_classic() +
  theme(legend.position = "top") +
  scale_color_brewer(palette = "Set2", labels = c("18 y/o", "29 y/o", "32 y/o", "59 y/o")) +
  guides(y = guide_axis(minor.ticks = TRUE), x = guide_axis(minor.ticks = TRUE)) +
  scale_y_continuous(minor_breaks = scales::breaks_width(10)) +
  scale_x_continuous(minor_breaks = scales::breaks_width(10)) 
ggsave("sum2dpctchangeScene5.tiff", device = "tiff", dpi = 300, width = 177, height = 150, units = "mm")


summary_all_2dtotal <- data_trimmed %>%
  mutate(age = fct_relevel(age, "18yo", "29yo", "32yo", "59yo"),
         scene = fct_relevel(scene, "scene 1", "scene 2", "scene 3", "scene 4", "scene 5")) %>%
  group_by(age, scene) %>%
  summarise(
    mean_dens = mean(abs(sum_length_2d_pct_change), na.rm = TRUE),
    sd_dens = sd(abs(sum_length_2d_pct_change), na.rm = TRUE),
    cv_dens = sd_dens / mean_dens
  )

model <- lm(cv_dens ~ age, summary_all_2dtotal)
summary(model)
Anova(model)
em1 <- emmeans(model, ~ age)
pairs(em1)

simulationOutput <- simulateResiduals(fittedModel = model, plot = T)

ggplot(mutate(live18yo, scene = fct_relevel(scene, "scene 1", "scene 2", "scene 3", "scene 4", "scene 5")), aes(time_min, sum_length_2d_pct_change, color = scene)) + 
  geom_line() +
  geom_point(size = 0.5, position = position_jitterdodge()) +
  labs(x = "Time (min)", y = "Total 2D filament length (μm)", color = "Cell") +
  theme_classic() +
  theme(legend.position = "top") +
  scale_color_prism(palette = "floral", , labels = c("Cell 1", "Cell 2", "Cell 3", "Cell 4", "Cell 5")) +
  guides(y = guide_axis(minor.ticks = TRUE), x = guide_axis(minor.ticks = TRUE)) +
  scale_y_continuous(minor_breaks = scales::breaks_width(10)) +
  scale_x_continuous(minor_breaks = scales::breaks_width(10)) 
ggsave("avgfilwidth29.tiff", device = "tiff", dpi = 300, width = 177, height = 150, units = "mm")


ggplot(mutate(live29yo, scene = fct_relevel(scene, "scene 1", "scene 2", "scene 3", "scene 4", "scene 5")), aes(time_min, sum_length_2d, color = scene)) + 
  geom_line() +
  geom_point(size = 0.5, position = position_jitterdodge()) +
  labs(x = "Time (min)", y = "Total 2D filament length (μm)", color = "Cell") +
  theme_classic() +
  theme(legend.position = "top") +
  scale_color_prism(palette = "floral", , labels = c("Cell 1", "Cell 2", "Cell 3", "Cell 4", "Cell 5")) +
  guides(y = guide_axis(minor.ticks = TRUE), x = guide_axis(minor.ticks = TRUE)) +
  scale_y_continuous(minor_breaks = scales::breaks_width(1000)) +
  scale_x_continuous(minor_breaks = scales::breaks_width(10)) 
ggsave("avgfilwidth29.tiff", device = "tiff", dpi = 300, width = 177, height = 150, units = "mm")


# CV plots

ggplot(mutate(summary_all_2dtotal, age = fct_relevel(age, "18yo", "29yo", "32yo", "59yo")), aes(age, cv_dens, color = age, fill = age)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_point(size = 1.5) +
  labs(x = "Age", y = "Coefficient of variance for total 2D filament length", color = "Age") +
  theme_classic() +
  scale_color_prism(palette = "floral", labels = c("18 y/o", "29 y/o", "32 y/o", "59 y/o")) +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('18 y/o', '29 y/o', '32 y/o', '59 y/o')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(0.1)) 
ggsave("/Volumes/T7/Live data quantitative analysis/final graphs/CV_density_actin_trim.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")
