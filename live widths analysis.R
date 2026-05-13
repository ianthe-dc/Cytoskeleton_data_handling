library(tidyverse)
library(lme4)
library(car)
library(emmeans)
library(ggprism)
library(forcats)
library(dplyr)
library(ggsignif)


actin_live_widths <- read.csv("actin_fil_widths_tidy.csv")
scene1_widths <- filter(actin_live_widths, scene %in% c("scene 1"))
cell_actin_live_widths <- read.csv("actin widths tp level.csv")
scene1_avgwidths <- filter(cell_actin_live_widths, scene %in% c("scene 1"))
scene2_avgwidths <- filter(cell_actin_live_widths, scene %in% c("scene 2"))
scene3_avgwidths <- filter(cell_actin_live_widths, scene %in% c("scene 3"))
scene4_avgwidths <- filter(cell_actin_live_widths, scene %in% c("scene 4"))
scene5_avgwidths <- filter(cell_actin_live_widths, scene %in% c("scene 5"))

yo29_avg_widths <- filter(cell_actin_live_widths, age %in% c("29yo"))

ggplot(mutate(scene1_avgwidths, scene = fct_relevel(age, "18yo", "29yo", "32yo", "59yo")), aes(time_min, mean_filament_width_pct_change, color = age)) + 
  geom_line() +
  geom_point(size = 0.5, position = position_jitterdodge()) +
  labs(x = "Time (min)", y = "Change in filament width (%)", color = "Cell") +
  theme_classic() +
  theme(legend.position = "top") +
  scale_color_brewer(palette = "Set2", , labels = c("18 y/o", "29 y/o", "32 y/o", "59 y/o")) +
  guides(y = guide_axis(minor.ticks = TRUE), x = guide_axis(minor.ticks = TRUE)) +
  scale_y_continuous(minor_breaks = scales::breaks_width(1)) +
  scale_x_continuous(minor_breaks = scales::breaks_width(10)) 
ggsave("avgfilwidthPctchgSc1.tiff", device = "tiff", dpi = 300, width = 177, height = 150, units = "mm")

ggplot(mutate(scene2_avgwidths, scene = fct_relevel(age, "18yo", "29yo", "32yo", "59yo")), aes(time_min, mean_filament_width_pct_change, color = age)) + 
  geom_line() +
  geom_point(size = 0.5, position = position_jitterdodge()) +
  labs(x = "Time (min)", y = "Change in filament width (%)", color = "Cell") +
  theme_classic() +
  theme(legend.position = "top") +
  scale_color_brewer(palette = "Set2", , labels = c("18 y/o", "29 y/o", "32 y/o", "59 y/o")) +
  guides(y = guide_axis(minor.ticks = TRUE), x = guide_axis(minor.ticks = TRUE)) +
  scale_y_continuous(minor_breaks = scales::breaks_width(1)) +
  scale_x_continuous(minor_breaks = scales::breaks_width(10)) 
ggsave("avgfilwidthPctchgSc2.tiff", device = "tiff", dpi = 300, width = 177, height = 150, units = "mm")

ggplot(mutate(scene3_avgwidths, scene = fct_relevel(age, "18yo", "29yo", "32yo", "59yo")), aes(time_min, mean_filament_width_pct_change, color = age)) + 
  geom_line() +
  geom_point(size = 0.5, position = position_jitterdodge()) +
  labs(x = "Time (min)", y = "Change in filament width (%)", color = "Cell") +
  theme_classic() +
  theme(legend.position = "top") +
  scale_color_brewer(palette = "Set2", , labels = c("18 y/o", "29 y/o", "32 y/o", "59 y/o")) +
  guides(y = guide_axis(minor.ticks = TRUE), x = guide_axis(minor.ticks = TRUE)) +
  scale_y_continuous(minor_breaks = scales::breaks_width(1)) +
  scale_x_continuous(minor_breaks = scales::breaks_width(10)) 
ggsave("avgfilwidthPctchgSc3.tiff", device = "tiff", dpi = 300, width = 177, height = 150, units = "mm")

ggplot(mutate(scene4_avgwidths, scene = fct_relevel(age, "18yo", "29yo", "59yo")), aes(time_min, mean_filament_width_pct_change, color = age)) + 
  geom_line() +
  geom_point(size = 0.5, position = position_jitterdodge()) +
  labs(x = "Time (min)", y = "Change in filament width (%)", color = "Cell") +
  theme_classic() +
  theme(legend.position = "top") +
  scale_color_brewer(palette = "Set2", , labels = c("18 y/o", "29 y/o", "59 y/o")) +
  guides(y = guide_axis(minor.ticks = TRUE), x = guide_axis(minor.ticks = TRUE)) +
  scale_y_continuous(minor_breaks = scales::breaks_width(1)) +
  scale_x_continuous(minor_breaks = scales::breaks_width(10)) 
ggsave("avgfilwidthPctchgSc4.tiff", device = "tiff", dpi = 300, width = 177, height = 150, units = "mm")

ggplot(mutate(scene5_avgwidths, scene = fct_relevel(age, "18yo", "29yo", "32yo", "59yo")), aes(time_min, mean_filament_width_pct_change, color = age)) + 
  geom_line() +
  geom_point(size = 0.5, position = position_jitterdodge()) +
  labs(x = "Time (min)", y = "Change in filament width (%)", color = "Cell") +
  theme_classic() +
  theme(legend.position = "top") +
  scale_color_brewer(palette = "Set2", , labels = c("18 y/o", "29 y/o", "32 y/o", "59 y/o")) +
  guides(y = guide_axis(minor.ticks = TRUE), x = guide_axis(minor.ticks = TRUE)) +
  scale_y_continuous(minor_breaks = scales::breaks_width(1)) +
  scale_x_continuous(minor_breaks = scales::breaks_width(10)) 
ggsave("avgfilwidthPctchgSc5.tiff", device = "tiff", dpi = 300, width = 177, height = 150, units = "mm")
