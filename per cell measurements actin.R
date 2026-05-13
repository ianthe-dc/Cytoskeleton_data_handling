library(tidyverse)
library(lme4)
library(car)
library(emmeans)
library(ggprism)
library(forcats)
library(dplyr)
library(ggsignif)
library(DHARMa)

actin_df <- read.csv("actin cell meas.csv")
actin_BJ <- filter(actin_df, condition_clean %in% c("bj p8.5 act", "bj p9.5 etop vim", "bj p11.5 act", "bj p12.5 -uv act", "bj p12.5 +uv act", "bj p26.5 act", "bj p31.5 act", "bj p36.5 act", "bj p38.5 act", "bj p40.5 act"))
actin_primary <- filter(actin_df, condition_clean %in% c("primary fibro 18yo act", "primary fibro 29yo act", "primary fibro 32yo act", "primary fibro 59yo act"))
actin_UV <- filter(actin_BJ, condition_clean %in% c("bj p12.5 -uv act", "bj p12.5 +uv act"))
actin_etop <- filter(actin_BJ, condition_clean %in% c("bj p8.5 act", "bj p9.5 etop vim"))
actin_pass <- filter(actin_df, condition_clean %in% c("bj p11.5 act", "bj p26.5 act", "bj p31.5 act", "bj p36.5 act", "bj p38.5 act", "bj p40.5 act"))


#stats for passage only
model_pass <- lmer(Branch.Point.Density...um.3. ~ condition_clean + (1|condition_clean:biorep:techrep), actin_pass)
summary(model_pass)
Anova(model_pass)
em_pass <- emmeans(model_pass, ~condition_clean)
pairs(em_pass)

simulationOutput <- simulateResiduals(fittedModel = model_pass, plot = T)

kruskal <- kruskal.test(Branch.Point.Density...um.3. ~ condition_clean, data = actin_pass)
kruskal
pairwise.wilcox.test(actin_pass$Branch.Point.Density...um.3., actin_pass$condition_clean, p.adjust.method = "BH")


#stats for UV only
model_UV <- lmer(Branch.Point.Density...um.3. ~ condition_clean + (1|condition_clean:biorep:techrep), actin_UV)
summary(model_UV)
Anova(model_UV)
em_UV <- emmeans(model_UV, ~condition_clean)
pairs(em_UV)

simulationOutput <- simulateResiduals(fittedModel = model_UV, plot = T)


#stats for etop only
model_etop <- lmer(Branch.Point.Density...um.3. ~ condition_clean + (1|condition_clean:biorep:techrep), actin_etop)
summary(model_etop)
Anova(model_etop)
em_etop <- emmeans(model_etop, ~condition_clean)
pairs(em_etop)

simulationOutput <- simulateResiduals(fittedModel = model_etop, plot = T)

kruskal <- kruskal.test(Branch.Point.Density...um.3. ~ condition_clean, data = actin_etop)
kruskal

#stats for actin primary
model13 <- lm(Branch.Point.Density...um.3. ~ condition_clean + techrep, actin_primary)
summary(model13)
Anova(model13)
em13 <- emmeans(model13, ~condition_clean)
pairs(em13)

simulationOutput <- simulateResiduals(fittedModel = model13, plot = T)

kruskal <- kruskal.test(Branch.Point.Density...um.3. ~ condition_clean, data = actin_primary)
kruskal
pairwise.wilcox.test(actin_primary$Branch.Point.Density...um.3., actin_primary$condition_clean, p.adjust.method = "BH")


# boxplot for actin_pass Cell.Volume..um.3.
summary_stats <- actin_pass %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p11.5 act","bj p26.5 act","bj p31.5 act",
                                       "bj p36.5 act","bj p38.5 act","bj p40.5 act")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(Cell.Volume..um.3., na.rm = TRUE),
    sd = sd(Cell.Volume..um.3., na.rm = TRUE),
    ymax = max(Cell.Volume..um.3., na.rm = TRUE)
  )

ggplot(mutate(actin_pass, condition_clean = fct_relevel(condition_clean, "bj p11.5 act", "bj p26.5 act", "bj p31.5 act", "bj p36.5 act", "bj p38.5 act", "bj p40.5 act")), aes(condition_clean, Cell.Volume..um.3., fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Passage", y = "Cell volume (μm^3)") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('p11.5', 'p26.5', 'p31.5', 'p36.5', 'p38.5', 'p40.5')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(2000)) +
  geom_text(
    data = summary_stats,
    aes(
      x = condition_clean,
      y = 53500,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 2.5
  ) +
  geom_signif(data = actin_pass, comparisons = list(c("bj p11.5 act","bj p36.5 act")), map_signif_level = TRUE, annotations = "**", y_position = 55000, vjust = 0.7, color = 'black') +
  geom_signif(data = actin_pass, comparisons = list(c("bj p11.5 act","bj p38.5 act")), map_signif_level = TRUE, annotations = "*", vjust = 0.7, y_position = 57000, color = 'black') +
  geom_signif(data = actin_pass, comparisons = list(c("bj p11.5 act","bj p40.5 act")), map_signif_level = TRUE, annotations = "***", y_position = 59000, vjust = 0.7, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/actPassCellSize.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")

# boxplot for actin_pass Skeleton.Density
summary_stats <- actin_pass %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p11.5 act","bj p26.5 act","bj p31.5 act",
                                       "bj p36.5 act","bj p38.5 act","bj p40.5 act")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(Skeleton.Density, na.rm = TRUE),
    sd = sd(Skeleton.Density, na.rm = TRUE),
    ymax = max(Skeleton.Density, na.rm = TRUE)
  )

ggplot(mutate(actin_pass, condition_clean = fct_relevel(condition_clean, "bj p11.5 act", "bj p26.5 act", "bj p31.5 act", "bj p36.5 act", "bj p38.5 act", "bj p40.5 act")), aes(condition_clean, Skeleton.Density, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Passage", y = "Network density") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('p11.5', 'p26.5', 'p31.5', 'p36.5', 'p38.5', 'p40.5')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(1000)) +
  geom_text(
    data = summary_stats,
    aes(
      x = condition_clean,
      y = 1.6,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = actin_pass, comparisons = list(c("bj p11.5 act","bj p40.5 act")), map_signif_level = TRUE, annotations = "ns", y_position = 1.65, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/actPassDensity.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")

# boxplot for actin_pass Branch.Ratio
summary_stats <- actin_pass %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p11.5 act","bj p26.5 act","bj p31.5 act",
                                       "bj p36.5 act","bj p38.5 act","bj p40.5 act")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(Branch.Ratio, na.rm = TRUE),
    sd = sd(Branch.Ratio, na.rm = TRUE),
    ymax = max(Branch.Ratio, na.rm = TRUE)
  )

ggplot(mutate(actin_pass, condition_clean = fct_relevel(condition_clean, "bj p11.5 act", "bj p26.5 act", "bj p31.5 act", "bj p36.5 act", "bj p38.5 act", "bj p40.5 act")), aes(condition_clean, Branch.Ratio, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Passage", y = "Branch ratio (Branch points/Filaments)") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('p11.5', 'p26.5', 'p31.5', 'p36.5', 'p38.5', 'p40.5')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(0.01)) +
  geom_text(
    data = summary_stats,
    aes(
      x = condition_clean,
      y = 1.75,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = actin_pass, comparisons = list(c("bj p11.5 act","bj p36.5 act")), map_signif_level = TRUE, annotations = "**", y_position = 1.76, vjust = 0.7, color = 'black') +
  geom_signif(data = actin_pass, comparisons = list(c("bj p11.5 act","bj p40.5 act")), map_signif_level = TRUE, annotations = "**", y_position = 1.77, vjust = 0.7, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/actPassBranchRatio.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")

# boxplot for actin_pass Branch.Point.Density...um.3.
summary_stats <- actin_pass %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p11.5 act","bj p26.5 act","bj p31.5 act",
                                       "bj p36.5 act","bj p38.5 act","bj p40.5 act")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(Branch.Point.Density...um.3., na.rm = TRUE),
    sd = sd(Branch.Point.Density...um.3., na.rm = TRUE),
    ymax = max(Branch.Point.Density...um.3., na.rm = TRUE)
  )

ggplot(mutate(actin_pass, condition_clean = fct_relevel(condition_clean, "bj p11.5 act", "bj p26.5 act", "bj p31.5 act", "bj p36.5 act", "bj p38.5 act", "bj p40.5 act")), aes(condition_clean, Branch.Point.Density...um.3., fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Passage", y = "Branch point density (Branch points/μm^3)") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('p11.5', 'p26.5', 'p31.5', 'p36.5', 'p38.5', 'p40.5')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(0.2)) +
  geom_text(
    data = summary_stats,
    aes(
      x = condition_clean,
      y = 4.65,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = actin_pass, comparisons = list(c("bj p11.5 act","bj p36.5 act")), map_signif_level = TRUE, annotations = "***", y_position = 4.8, vjust = 0.7, color = 'black') +
  geom_signif(data = actin_pass, comparisons = list(c("bj p11.5 act","bj p38.5 act")), map_signif_level = TRUE, annotations = "***", y_position = 5, vjust = 0.7, color = 'black') +
  geom_signif(data = actin_pass, comparisons = list(c("bj p11.5 act","bj p40.5 act")), map_signif_level = TRUE, annotations = "***", y_position = 5.2, vjust = 0.7, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/actPassBranchDensity.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")

# boxplot for actin_pass extras not plotted
ggplot(mutate(actin_pass, condition_clean = fct_relevel(condition_clean, "bj p11.5 act", "bj p26.5 act", "bj p31.5 act", "bj p36.5 act", "bj p38.5 act", "bj p40.5 act")), aes(condition_clean, Average.Filament.Width..um., fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Passage", y = "Branch point density (Branch points/ROI)") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('p11.5', 'p26.5', 'p31.5', 'p36.5', 'p38.5', 'p40.5')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(0.2)) +
  geom_signif(data = actin_pass, comparisons = list(c("bj p11.5 act","bj p36.5 act")), map_signif_level = TRUE, annotations = "***", y_position = , vjust = 0.7, color = 'black') +
  geom_signif(data = actin_pass, comparisons = list(c("bj p11.5 act","bj p38.5 act")), map_signif_level = TRUE, annotations = "***", y_position = , vjust = 0.7, color = 'black') +
  geom_signif(data = actin_pass, comparisons = list(c("bj p11.5 act","bj p40.5 act")), map_signif_level = TRUE, annotations = "***", y_position = , vjust = 0.7, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/actPassBranchDensity.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")



# boxplot for actin_UV Cell.Volume..um.3.
summary_stats <- actin_UV %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p12.5 -uv act", "bj p12.5 +uv act")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(Cell.Volume..um.3., na.rm = TRUE),
    sd = sd(Cell.Volume..um.3., na.rm = TRUE),
    ymax = max(Cell.Volume..um.3., na.rm = TRUE)
  )

ggplot(mutate(actin_UV, condition_clean = fct_relevel(condition_clean, "bj p12.5 -uv act", "bj p12.5 +uv act")), aes(condition_clean, Cell.Volume..um.3., fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Condition", y = "Cell volume (μm^3)") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('p12.5 -UV', 'p12.5 +UV')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(1000)) +
  geom_text(
    data = summary_stats,
    aes(
      x = condition_clean,
      y = 46000,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = actin_UV, comparisons = list(c("bj p12.5 -uv act", "bj p12.5 +uv act")), annotations = "ns", map_signif_level = TRUE, y_position = 48000, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/actUVCellSize.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")

# boxplot for actin_UV Skeleton.Density
summary_stats <- actin_UV %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p12.5 -uv act", "bj p12.5 +uv act")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(Skeleton.Density, na.rm = TRUE),
    sd = sd(Skeleton.Density, na.rm = TRUE),
    ymax = max(Skeleton.Density, na.rm = TRUE)
  )

ggplot(mutate(actin_UV, condition_clean = fct_relevel(condition_clean, "bj p12.5 -uv act", "bj p12.5 +uv act")), aes(condition_clean, Skeleton.Density, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Condition", y = "Network density") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('p12.5 -UV', 'p12.5 +UV')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(0.1)) +
  geom_text(
    data = summary_stats,
    aes(
      x = condition_clean,
      y = 1.28,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = actin_UV, comparisons = list(c("bj p12.5 -uv act", "bj p12.5 +uv act")), annotations = "ns", map_signif_level = TRUE, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/actUVDensity.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")

# boxplot for actin_UV Branch.Ratio
summary_stats <- actin_UV %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p12.5 -uv act", "bj p12.5 +uv act")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(Branch.Ratio, na.rm = TRUE),
    sd = sd(Branch.Ratio, na.rm = TRUE),
    ymax = max(Branch.Ratio, na.rm = TRUE)
  )

ggplot(mutate(actin_UV, condition_clean = fct_relevel(condition_clean, "bj p12.5 -uv act", "bj p12.5 +uv act")), aes(condition_clean, Branch.Ratio, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Condition", y = "Branch ratio (Branch points/Filaments)") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('p12.5 -UV', 'p12.5 +UV')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(0.005)) +
  geom_text(
    data = summary_stats,
    aes(
      x = condition_clean,
      y = 1.684,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = actin_UV, comparisons = list(c("bj p12.5 -uv act", "bj p12.5 +uv act")), annotations = "ns", map_signif_level = TRUE, y_position = 1.69, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/actUVBranchRatio.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")

# boxplot for actin_UV Branch.Point.Density...um.3.
summary_stats <- actin_UV %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p12.5 -uv act", "bj p12.5 +uv act")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(Branch.Point.Density...um.3., na.rm = TRUE),
    sd = sd(Branch.Point.Density...um.3., na.rm = TRUE),
    ymax = max(Branch.Point.Density...um.3., na.rm = TRUE)
  )

ggplot(mutate(actin_UV, condition_clean = fct_relevel(condition_clean, "bj p12.5 -uv act", "bj p12.5 +uv act")), aes(condition_clean, Branch.Point.Density...um.3., fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Condition", y = "Branch point density (Branch points/μm^3)") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('p12.5 -UV', 'p12.5 +UV')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(0.2)) +
  geom_text(
    data = summary_stats,
    aes(
      x = condition_clean,
      y = 3,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = actin_UV, comparisons = list(c("bj p12.5 -uv act", "bj p12.5 +uv act")), annotations = "ns", map_signif_level = TRUE, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/actUVBranchDensity.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")


# boxplot for actin_etop Cell.Volume..um.3.
summary_stats <- actin_etop %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p8.5 act", "bj p9.5 etop vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(Cell.Volume..um.3., na.rm = TRUE),
    sd = sd(Cell.Volume..um.3., na.rm = TRUE),
    ymax = max(Cell.Volume..um.3., na.rm = TRUE)
  )

ggplot(mutate(actin_etop, condition_clean = fct_relevel(condition_clean, "bj p8.5 act", "bj p9.5 etop vim")), aes(condition_clean, Cell.Volume..um.3., fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Condition", y = "Cell volume (μm^3)") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('p8.5 -etop', 'p9.5 +etop')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(5000)) +
  geom_text(
    data = summary_stats,
    aes(
      x = condition_clean,
      y = 72000,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = actin_etop, comparisons = list(c("bj p8.5 act","bj p9.5 etop vim")), annotations = "*", map_signif_level = TRUE, y_position = 75000, vjust = 0.7, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/actEtopCellSize.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")

# boxplot for actin_etop Skeleton.Density
summary_stats <- actin_etop %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p8.5 act", "bj p9.5 etop vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(Skeleton.Density, na.rm = TRUE),
    sd = sd(Skeleton.Density, na.rm = TRUE),
    ymax = max(Skeleton.Density, na.rm = TRUE)
  )

ggplot(mutate(actin_etop, condition_clean = fct_relevel(condition_clean, "bj p8.5 act", "bj p9.5 etop vim")), aes(condition_clean, Skeleton.Density, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Condition", y = "Network density") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('p8.5 -etop', 'p9.5 +etop')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(0.05)) +
  geom_text(
    data = summary_stats,
    aes(
      x = condition_clean,
      y = 0.7,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = actin_etop, comparisons = list(c("bj p8.5 act","bj p9.5 etop vim")), annotations = "***", map_signif_level = TRUE, vjust = 0.7, y_position = 0.72, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/actEtopDensity.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")

# boxplot for actin_etop Branch.Ratio
summary_stats <- actin_etop %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p8.5 act", "bj p9.5 etop vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(Branch.Ratio, na.rm = TRUE),
    sd = sd(Branch.Ratio, na.rm = TRUE),
    ymax = max(Branch.Ratio, na.rm = TRUE)
  )

ggplot(mutate(actin_etop, condition_clean = fct_relevel(condition_clean, "bj p8.5 act", "bj p9.5 etop vim")), aes(condition_clean, Branch.Ratio, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Condition", y = "Branch ratio (Branch points/Filaments)") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('p8.5 -etop', 'p9.5 +etop')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(0.01)) +
  geom_text(
    data = summary_stats,
    aes(
      x = condition_clean,
      y = 1.735,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = actin_etop, comparisons = list(c("bj p8.5 act","bj p9.5 etop vim")), annotations = "**", map_signif_level = TRUE, vjust = 0.7, y_position = 1.75, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/actEtopBranchRatio.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")

# boxplot for actin_etop Branch.Point.Density...um.3.
summary_stats <- actin_etop %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p8.5 act", "bj p9.5 etop vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(Branch.Point.Density...um.3., na.rm = TRUE),
    sd = sd(Branch.Point.Density...um.3., na.rm = TRUE),
    ymax = max(Branch.Point.Density...um.3., na.rm = TRUE)
  )

ggplot(mutate(actin_etop, condition_clean = fct_relevel(condition_clean, "bj p8.5 act", "bj p9.5 etop vim")), aes(condition_clean, Branch.Point.Density...um.3., fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Condition", y = "Branch point density (Branch points/μm^3)") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('p8.5 -etop', 'p9.5 +etop')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(0.5)) +
  geom_text(
    data = summary_stats,
    aes(
      x = condition_clean,
      y = 10,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = actin_etop, comparisons = list(c("bj p8.5 act","bj p9.5 etop vim")), annotations = "ns", map_signif_level = TRUE, y_position = 10.5, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/actEtopBranchDensity.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")


# boxplot for actin_primary Cell.Volume..um.3.
summary_stats <- actin_primary %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "primary fibro 18yo act", "primary fibro 29yo act", "primary fibro 32yo act", "primary fibro 59yo act")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(Cell.Volume..um.3., na.rm = TRUE),
    sd = sd(Cell.Volume..um.3., na.rm = TRUE),
    ymax = max(Cell.Volume..um.3., na.rm = TRUE)
  )

ggplot(mutate(actin_primary, condition_clean = fct_relevel(condition_clean, "primary fibro 18yo act", "primary fibro 29yo act", "primary fibro 32yo act", "primary fibro 59yo act")), aes(condition_clean, Cell.Volume..um.3., fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Age", y = "Cell volume (μm^3)") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('18 y/o', '29 y/o', '32 y/o', '59 y/o')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(5000)) +
  geom_text(
    data = summary_stats,
    aes(
      x = condition_clean,
      y = 89000,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = actin_primary, comparisons = list(c("primary fibro 18yo act", "primary fibro 32yo act")), annotations = "*", map_signif_level = TRUE, y_position = 95000, vjust = 0.7, color = 'black') +
  geom_signif(data = actin_primary, comparisons = list(c("primary fibro 18yo act", "primary fibro 59yo act")), annotations = "***", map_signif_level = TRUE, vjust = 0.7, color = 'black') 
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/actPrimCellSize.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")

# boxplot for actin_primary Skeleton.Density
summary_stats <- actin_primary %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "primary fibro 18yo act", "primary fibro 29yo act", "primary fibro 32yo act", "primary fibro 59yo act")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(Skeleton.Density, na.rm = TRUE),
    sd = sd(Skeleton.Density, na.rm = TRUE),
    ymax = max(Skeleton.Density, na.rm = TRUE)
  )

ggplot(mutate(actin_primary, condition_clean = fct_relevel(condition_clean, "primary fibro 18yo act", "primary fibro 29yo act", "primary fibro 32yo act", "primary fibro 59yo act")), aes(condition_clean, Skeleton.Density, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Age", y = "Network density") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('18 y/o', '29 y/o', '32 y/o', '59 y/o')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(0.1)) +
  geom_text(
    data = summary_stats,
    aes(
      x = condition_clean,
      y = 1.4,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3.5
  ) +
  geom_signif(data = actin_primary, comparisons = list(c("primary fibro 18yo act", "primary fibro 59yo act")), annotations = "ns", map_signif_level = TRUE, y_position = 1.5, color = 'black') 
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/actPrimDensity.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")

# boxplot for actin_primary Branch.Ratio
summary_stats <- actin_primary %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "primary fibro 18yo act", "primary fibro 29yo act", "primary fibro 32yo act", "primary fibro 59yo act")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(Branch.Ratio, na.rm = TRUE),
    sd = sd(Branch.Ratio, na.rm = TRUE),
    ymax = max(Branch.Ratio, na.rm = TRUE)
  )

ggplot(mutate(actin_primary, condition_clean = fct_relevel(condition_clean, "primary fibro 18yo act", "primary fibro 29yo act", "primary fibro 32yo act", "primary fibro 59yo act")), aes(condition_clean, Branch.Ratio, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Age", y = "Branch ratio (Branch points/Filaments)") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('18 y/o', '29 y/o', '32 y/o', '59 y/o')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(0.005)) +
  geom_text(
    data = summary_stats,
    aes(
      x = condition_clean,
      y = 1.695,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3.5
  ) +
  geom_signif(data = actin_primary, comparisons = list(c("primary fibro 18yo act", "primary fibro 59yo act")), annotations = "ns", map_signif_level = TRUE, y_position = 1.7, color = 'black') 
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/actPrimBranchRatio.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")

# boxplot for actin_primary Branch.Point.Density...um.3.
summary_stats <- actin_primary %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "primary fibro 18yo act", "primary fibro 29yo act", "primary fibro 32yo act", "primary fibro 59yo act")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(Branch.Point.Density...um.3., na.rm = TRUE),
    sd = sd(Branch.Point.Density...um.3., na.rm = TRUE),
    ymax = max(Branch.Point.Density...um.3., na.rm = TRUE)
  )

ggplot(mutate(actin_primary, condition_clean = fct_relevel(condition_clean, "primary fibro 18yo act", "primary fibro 29yo act", "primary fibro 32yo act", "primary fibro 59yo act")), aes(condition_clean, Branch.Point.Density...um.3., fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Age", y = "Branch point density (Branch points/μm^3)") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('18 y/o', '29 y/o', '32 y/o', '59 y/o')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(0.2)) +
  geom_text(
    data = summary_stats,
    aes(
      x = condition_clean,
      y = 3.95,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3.5
  ) +
  geom_signif(data = actin_primary, comparisons = list(c("primary fibro 18yo act", "primary fibro 59yo act")), annotations = "ns", map_signif_level = TRUE, y_position = 4.1, color = 'black') 
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/actPrimBranchDensity.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")

