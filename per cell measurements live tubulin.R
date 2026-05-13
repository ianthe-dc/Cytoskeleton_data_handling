library(tidyverse)
library(lme4)
library(car)
library(emmeans)
library(ggprism)
library(forcats)
library(dplyr)
library(ggsignif)
library(DHARMa)
library(flextable)

tubulin_live <- read.csv("tubulin live cell meas.csv")

cols <- c(
  "Cell.Volume..um.3.",
  "Skeleton.Density",
  "Branch.Ratio",
  "Branch.Point.Density...um.3."
)

tubulin_live <- tubulin_live %>%
  group_by(age, scene) %>%
  mutate(across(all_of(cols),
                ~ (. - .[time_min == 0][1]) / .[time_min == 0][1] * 100,
                .names = "{.col}_pct_change"
  )) %>%
  ungroup()

write.csv(tubulin_live, "/Users/vbss76/Desktop/CytoAge/Thesis/image analysis/live analysis/tubulin live cell meas pct changes.csv", row.names = FALSE)

tubulin_live_pct <- read.csv("tubulin live cell meas pct changes.csv")


scene1 <- filter(tubulin_live_pct, scene %in% c("scene 1"))
scene2 <- filter(tubulin_live_pct, scene %in% c("scene 2"))
scene3 <- filter(tubulin_live_pct, scene %in% c("scene 3"))
scene4 <- filter(tubulin_live_pct, scene %in% c("scene 4"))
scene5 <- filter(tubulin_live_pct, scene %in% c("scene 5"))
live18yo <- filter(data_trimmed, age %in% "18yo")
live29yo <- filter(data_trimmed, age %in% "29yo")
live32yo <- filter(data_trimmed, age %in% "32yo")
live59yo <- filter(data_trimmed, age %in% "59yo")


max_common_time <- tubulin_live_pct %>%
  group_by(scene, age) %>%
  summarise(max_time = max(time_min)) %>%
  summarise(min(max_time)) %>%
  pull()

data_trimmed <- tubulin_live_pct %>%
  filter(time_min <= 70)


# summary stats for network density

summary_all_density <- data_trimmed %>%
  mutate(age = fct_relevel(age, "18yo", "29yo", "32yo", "59yo"),
         scene = fct_relevel(scene, "scene 1", "scene 2", "scene 3", "scene 4", "scene 5", "scene 6")) %>%
  group_by(age, scene) %>%
  summarise(
    mean_dens = mean(abs(Skeleton.Density_pct_change), na.rm = TRUE),
    sd_dens = sd(Skeleton.Density, na.rm = TRUE),
    cv_dens = sd_dens/mean_dens
  )

summary_all_density

model <- lm(cv_dens ~ age, summary_all_density)
summary(model)
Anova(model)
em1 <- emmeans(model, ~ age)
pairs(em1)

simulationOutput <- simulateResiduals(fittedModel = model, plot = T)

summary_all_density_fmt <- summary_all_density %>%
  mutate(
    mean_dens = signif(mean_dens, 3),
    sd_dens   = signif(sd_dens, 3),
    cv_dens   = signif(cv_dens, 3)
  )

summary_all_density_fmt <- summary_all_density_fmt %>%
  mutate(
    mean_dens = format(mean_dens, trim = TRUE, scientific = FALSE),
    sd_dens   = format(sd_dens, trim = TRUE, scientific = FALSE),
    cv_dens   = format(cv_dens, trim = TRUE, scientific = FALSE)
  )

# table of means for network density

denst <- flextable(summary_all_density_fmt) %>%
  set_header_labels(
    age = "Age",
    scene = "Cell",
    mean_dens = "Mean",
    sd_dens = "SD",
    cv_dens = "CV"
  )

denst


save_as_image(denst, path = "/Volumes/T7/Live data quantitative analysis/final graphs/tubulin/summary_all_density.png", dpi = 300, width = 100, height = 190, units = "mm")


summary_all_stats_dens <- summary_all_density %>%
  mutate(age = fct_relevel(age, "18yo", "29yo", "32yo", "59yo")) %>%
  group_by(age) %>%
  summarise(
    mean_all_dens = mean(mean_dens, na.rm = TRUE),
    sd_all_dens = sd(mean_dens, na.rm = TRUE)
  )

# summary stats for branch point density

summary_all_bpd <- data_trimmed %>%
  mutate(age = fct_relevel(age, "18yo", "29yo", "32yo", "59yo"),
         scene = fct_relevel(scene, "scene 1", "scene 2", "scene 3", "scene 4", "scene 5")) %>%
  group_by(age, scene) %>%
  summarise(
    mean_bpd = mean(abs(Cell.Volume..um.3._pct_change), na.rm = TRUE),
    sd_bpd = sd(abs(Cell.Volume..um.3._pct_change), na.rm = TRUE),
    cv_bpd = sd_bpd / mean_bpd
  )

model2 <- lm(cv_bpd ~ age, summary_all_bpd)
summary(model2)
Anova(model2)
em2 <- emmeans(model2, ~ age)
pairs(em2)

simulationOutput <- simulateResiduals(fittedModel = model2, plot = T)


summary_all_stats_bpd <- summary_all_bpd %>%
  mutate(age = fct_relevel(age, "18yo", "29yo", "32yo", "59yo")) %>%
  group_by(age) %>%
  summarise(
    mean_all_bpd = mean(mean_bpd, na.rm = TRUE),
    sd_all_bpd = sd(mean_bpd, na.rm = TRUE)
  )


# Overall barplots

ggplot(mutate(summary_all_stats_dens, age = fct_relevel(age, "18yo", "29yo", "32yo", "59yo")), aes(age, mean_all_dens, fill = age, colour = age)) + 
  geom_bar(stat = "identity", width = 0.8) +
  geom_errorbar(aes(ymin = mean_all_dens - sd_all_dens, ymax = mean_all_dens + sd_all_dens), width = 0.2, size = 0.5) +
  labs(x = "Age", y = "Average density change (%)") +
  theme_bw() +
  theme(legend.position="none") +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(y = guide_axis(minor.ticks = TRUE), x = guide_axis(minor.ticks = TRUE)) 
ggsave("/Volumes/T7/Live data quantitative analysis/final graphs/tubulin/densitypctchangeAll.tiff", device = "tiff", dpi = 300, width = 125, height = 100, units = "mm")

ggplot(mutate(summary_all_stats_bpd, age = fct_relevel(age, "18yo", "29yo", "32yo", "59yo")), aes(age, mean_all_bpd, fill = age, colour = age)) + 
  geom_bar(stat = "identity", width = 0.8) +
  geom_errorbar(aes(ymin = mean_all_bpd - sd_all_bpd, ymax = mean_all_bpd + sd_all_bpd), width = 0.2, size = 0.5) +
  labs(x = "Age", y = "Average branch point density change (%)") +
  theme_bw() +
  theme(legend.position="none") +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(y = guide_axis(minor.ticks = TRUE), x = guide_axis(minor.ticks = TRUE)) 
ggsave("/Volumes/T7/Live data quantitative analysis/final graphs/tubulin/bpdensitypctchangeAll.tiff", device = "tiff", dpi = 300, width = 125, height = 100, units = "mm")



# density plots

ggplot(mutate(scene1, scene = fct_relevel(age, "18yo", "29yo", "32yo", "59yo")), aes(time_min, Skeleton.Density_pct_change, color = age)) + 
  geom_line() +
  geom_point(size = 0.5) +
  labs(x = "Time (min)", y = "Network density", color = "Age") +
  theme_classic() +
  theme(legend.position = "top") +
  scale_color_brewer(palette = "Set2", labels = c("18 y/o", "29 y/o", "32 y/o", "59 y/o")) +
  guides(y = guide_axis(minor.ticks = TRUE), x = guide_axis(minor.ticks = TRUE)) +
  scale_y_continuous(minor_breaks = scales::breaks_width(10)) +
  scale_x_continuous(minor_breaks = scales::breaks_width(10)) 
ggsave("/Volumes/T7/Live data quantitative analysis/final graphs/tubulin/densitypctchangeScene1.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")

ggplot(mutate(scene2, scene = fct_relevel(age, "18yo", "29yo", "32yo", "59yo")), aes(time_min, Skeleton.Density_pct_change, color = age)) + 
  geom_line() +
  geom_point(size = 0.5) +
  labs(x = "Time (min)", y = "Network density", color = "Age") +
  theme_classic() +
  theme(legend.position = "top") +
  scale_color_brewer(palette = "Set2", labels = c("18 y/o", "29 y/o", "32 y/o", "59 y/o")) +
  guides(y = guide_axis(minor.ticks = TRUE), x = guide_axis(minor.ticks = TRUE)) +
  scale_y_continuous(minor_breaks = scales::breaks_width(10)) +
  scale_x_continuous(minor_breaks = scales::breaks_width(10)) 
ggsave("/Volumes/T7/Live data quantitative analysis/final graphs/tubulin/densitypctchangeScene2.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")

ggplot(mutate(scene3, scene = fct_relevel(age, "18yo", "29yo", "32yo", "59yo")), aes(time_min, Skeleton.Density_pct_change, color = age)) + 
  geom_line() +
  geom_point(size = 0.5) +
  labs(x = "Time (min)", y = "Network density", color = "Age") +
  theme_classic() +
  theme(legend.position = "top") +
  scale_color_brewer(palette = "Set2", labels = c("18 y/o", "29 y/o", "32 y/o", "59 y/o")) +
  guides(y = guide_axis(minor.ticks = TRUE), x = guide_axis(minor.ticks = TRUE)) +
  scale_y_continuous(minor_breaks = scales::breaks_width(10)) +
  scale_x_continuous(minor_breaks = scales::breaks_width(10)) 
ggsave("/Volumes/T7/Live data quantitative analysis/final graphs/tubulin/densitypctchangeScene3.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")

ggplot(mutate(scene4, scene = fct_relevel(age, "18yo", "29yo", "32yo", "59yo")), aes(time_min, Skeleton.Density_pct_change, color = age)) + 
  geom_line() +
  geom_point(size = 0.5) +
  labs(x = "Time (min)", y = "Network density", color = "Age") +
  theme_classic() +
  theme(legend.position = "top") +
  scale_color_brewer(palette = "Set2", labels = c("18 y/o", "29 y/o", "32 y/o", "59 y/o")) +
  guides(y = guide_axis(minor.ticks = TRUE), x = guide_axis(minor.ticks = TRUE)) +
  scale_y_continuous(minor_breaks = scales::breaks_width(10)) +
  scale_x_continuous(minor_breaks = scales::breaks_width(10)) 
ggsave("/Volumes/T7/Live data quantitative analysis/final graphs/tubulin/densitypctchangeScene4.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")

ggplot(mutate(scene5, scene = fct_relevel(age, "18yo", "29yo", "32yo", "59yo")), aes(time_min, Skeleton.Density_pct_change, color = age)) + 
  geom_line() +
  geom_point(size = 0.5) +
  labs(x = "Time (min)", y = "Network density", color = "Age") +
  theme_classic() +
  theme(legend.position = "top") +
  scale_color_brewer(palette = "Set2", labels = c("18 y/o", "29 y/o", "32 y/o", "59 y/o")) +
  guides(y = guide_axis(minor.ticks = TRUE), x = guide_axis(minor.ticks = TRUE)) +
  scale_y_continuous(minor_breaks = scales::breaks_width(10)) +
  scale_x_continuous(minor_breaks = scales::breaks_width(10)) 
ggsave("/Volumes/T7/Live data quantitative analysis/final graphs/tubulin/densitypctchangeScene5.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")


# density plots by age


ggplot(mutate(live18yo, scene = fct_relevel(scene, "scene 1", "scene 2", "scene 3", "scene 4", "scene 5")), aes(time_min, Skeleton.Density_pct_change, color = scene)) + 
  geom_line() +
  geom_point(size = 0.5) +
  labs(x = "Time (min)", y = "Network density (% change)", color = "Cell") +
  theme_classic() +
  theme(legend.position = "top") +
  scale_color_brewer(palette = "Set2", labels = c("Cell 1", "Cell 2", "Cell 3", "Cell 4", "Cell 5")) +
  guides(y = guide_axis(minor.ticks = TRUE), x = guide_axis(minor.ticks = TRUE)) +
  scale_y_continuous(minor_breaks = scales::breaks_width(10), breaks = breaks_width(100), limits = c(-50, 340)) +
  scale_x_continuous(minor_breaks = scales::breaks_width(10)) 
ggsave("/Volumes/T7/Live data quantitative analysis/final graphs/tubulin/densitypctchange18yo_trim.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")

ggplot(mutate(live29yo, scene = fct_relevel(scene, "scene 1", "scene 2", "scene 3", "scene 4", "scene 5")), aes(time_min, Skeleton.Density_pct_change, color = scene)) + 
  geom_line() +
  geom_point(size = 0.5) +
  labs(x = "Time (min)", y = "Network density (% change)", color = "Cell") +
  theme_classic() +
  theme(legend.position = "top") +
  scale_color_brewer(palette = "Set2", labels = c("Cell 1", "Cell 2", "Cell 3", "Cell 4", "Cell 5")) +
  guides(y = guide_axis(minor.ticks = TRUE), x = guide_axis(minor.ticks = TRUE)) +
  scale_y_continuous(minor_breaks = scales::breaks_width(10), breaks = breaks_width(100), limits = c(-50, 340)) +
  scale_x_continuous(minor_breaks = scales::breaks_width(10)) 
ggsave("/Volumes/T7/Live data quantitative analysis/final graphs/tubulin/densitypctchange29yo_trim.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")

ggplot(mutate(live32yo, scene = fct_relevel(scene, "scene 1", "scene 2", "scene 3", "scene 4", "scene 5", "scene 6")), aes(time_min, Skeleton.Density_pct_change, color = scene)) + 
  geom_line() +
  geom_point(size = 0.5) +
  labs(x = "Time (min)", y = "Network density (% change)", color = "Cell") +
  theme_classic() +
  theme(legend.position = "top") +
  scale_color_brewer(palette = "Set2", labels = c("Cell 1", "Cell 2", "Cell 3", "Cell 4", "Cell 5", "Cell 6")) +
  guides(y = guide_axis(minor.ticks = TRUE), x = guide_axis(minor.ticks = TRUE)) +
  scale_y_continuous(minor_breaks = scales::breaks_width(10), breaks = breaks_width(100), limits = c(-50, 340)) +
  scale_x_continuous(minor_breaks = scales::breaks_width(10)) 
ggsave("/Volumes/T7/Live data quantitative analysis/final graphs/tubulin/densitypctchange32yo_trim.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")

ggplot(mutate(live59yo, scene = fct_relevel(scene, "scene 1", "scene 2", "scene 3", "scene 4", "scene 5")), aes(time_min, Skeleton.Density_pct_change, color = scene)) + 
  geom_line() +
  geom_point(size = 0.5) +
  labs(x = "Time (min)", y = "Network density (% change)", color = "Cell") +
  theme_classic() +
  theme(legend.position = "top") +
  scale_color_brewer(palette = "Set2", labels = c("Cell 1", "Cell 2", "Cell 3", "Cell 4", "Cell 5")) +
  guides(y = guide_axis(minor.ticks = TRUE), x = guide_axis(minor.ticks = TRUE)) +
  scale_y_continuous(minor_breaks = scales::breaks_width(10), breaks = breaks_width(100), limits = c(-50, 340)) +
  scale_x_continuous(minor_breaks = scales::breaks_width(10)) 
ggsave("/Volumes/T7/Live data quantitative analysis/final graphs/tubulin/densitypctchange59yo_trim.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")


# CV plots

ggplot(mutate(summary_all_density, age = fct_relevel(age, "18yo", "29yo", "32yo", "59yo")), aes(age, cv_dens, color = age, fill = age)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_point(size = 1.5) +
  labs(x = "Age", y = "Coefficient of variance for density", color = "Age") +
  theme_classic() +
  scale_color_prism(palette = "floral", labels = c("18 y/o", "29 y/o", "32 y/o", "59 y/o")) +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('18 y/o', '29 y/o', '32 y/o', '59 y/o')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(0.05)) 
ggsave("/Volumes/T7/Live data quantitative analysis/final graphs/CV_density_tubulin_OG.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")

ggplot(mutate(summary_all_bpd, age = fct_relevel(age, "18yo", "29yo", "32yo", "59yo")), aes(age, cv_bpd, color = age, fill = age)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_point(size = 1.5) +
  labs(x = "Age", y = "Coefficient of variance for branch density", color = "Age") +
  theme_classic() +
  scale_color_prism(palette = "floral", labels = c("18 y/o", "29 y/o", "32 y/o", "59 y/o")) +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('18 y/o', '29 y/o', '32 y/o', '59 y/o')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(0.1)) 
ggsave("/Volumes/T7/Live data quantitative analysis/final graphs/CV_branchdens_tubulin_trim.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")
