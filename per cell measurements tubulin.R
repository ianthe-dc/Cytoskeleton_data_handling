library(tidyverse)
library(lme4)
library(car)
library(emmeans)
library(ggprism)
library(forcats)
library(dplyr)
library(ggsignif)
library(DHARMa)

tubulin_df <- read.csv("tubulin cell meas.csv")
tubulin_BJ <- filter(tubulin_df, condition_clean %in% c("bj p11.5 tub", "bj p12.5 -uv tub", "bj p12.5 +uv tub", "bj p26.5 tub", "bj p31.5 tub", "bj p36.5 tub", "bj p38.5 tub", "bj p40.5 tub"))
tubulin_primary <- filter(tubulin_df, condition_clean %in% c("primary fibro 18yo tub", "primary fibro 29yo tub", "primary fibro 32yo tub", "primary fibro 59yo tub"))
tubulin_UV <- filter(tubulin_df, condition_clean %in% c("bj p12.5 -uv tub", "bj p12.5 +uv tub"))
tubulin_pass <- filter(tubulin_df, condition_clean %in% c("bj p11.5 tub", "bj p26.5 tub", "bj p31.5 tub", "bj p36.5 tub", "bj p38.5 tub", "bj p40.5 tub"))


#stats for passage only 
model_pass <- lmer(Skeleton.Density ~ condition_clean + (1|condition_clean:biorep:techrep), tubulin_pass)
summary(model_pass)
Anova(model_pass)
em_pass <- emmeans(model_pass, ~condition_clean)
pairs(em_pass)

simulationOutput <- simulateResiduals(fittedModel = model_pass, plot = T)

kruskal <- kruskal.test(Skeleton.Density ~ condition_clean, data = tubulin_pass)
kruskal
pairwise.wilcox.test(tubulin_pass$Branch.Point.Density...um.3., tubulin_pass$condition_clean, p.adjust.method = "BH")


#stats for tubulin_UV
model_UV <- lmer(Branch.Point.Density...um.3. ~ condition_clean + (1|condition_clean:biorep:techrep), tubulin_UV)
summary(model_UV)
Anova(model_UV)
em_UV <- emmeans(model_UV, ~condition_clean)
pairs(em_UV)

simulationOutput <- simulateResiduals(fittedModel = model_UV, plot = T)

#stats for tubulin_primary
model_prim <- lmer(Cell.Volume..um.3. ~ condition_clean + (1|condition_clean:biorep:techrep), tubulin_primary)
model_prim <- lm(Cell.Volume..um.3. ~ condition_clean + techrep, tubulin_primary)
summary(model_prim)
Anova(model_prim)
em_prim <- emmeans(model_prim, ~condition_clean)
pairs(em_prim)

simulationOutput <- simulateResiduals(fittedModel = model_prim, plot = T)

kruskal <- kruskal.test(Cell.Volume..um.3. ~ condition_clean, data = tubulin_primary)
kruskal
pairwise.wilcox.test(tubulin_primary$Cell.Volume..um.3., tubulin_primary$condition_clean, p.adjust.method = "BH")


# boxplot for tubulin_pass Cell.Volume..um.3.
summary_stats <- tubulin_pass %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p11.5 tub","bj p26.5 tub","bj p31.5 tub",
                                       "bj p36.5 tub","bj p38.5 tub","bj p40.5 tub")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(Cell.Volume..um.3., na.rm = TRUE),
    sd = sd(Cell.Volume..um.3., na.rm = TRUE),
    ymax = max(Cell.Volume..um.3., na.rm = TRUE)
  )

ggplot(mutate(tubulin_pass, condition_clean = fct_relevel(condition_clean, "bj p11.5 tub", "bj p26.5 tub", "bj p31.5 tub", "bj p36.5 tub", "bj p38.5 tub", "bj p40.5 tub")), aes(condition_clean, Cell.Volume..um.3., fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Passage", y = "Cell volume (μm^3)") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('p11.5', 'p26.5', 'p31.5', 'p36.5', 'p38.5', 'p40.5')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(5000)) +
  geom_text(
    data = summary_stats,
    aes(
      x = condition_clean,
      y = 75000,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 2.5
  ) +
  geom_signif(data = tubulin_pass, comparisons = list(c("bj p11.5 tub","bj p38.5 tub")), annotations = "*", map_signif_level = TRUE, y_position = 79000, vjust = 0.7, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/tubPassCellSize.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")

# boxplot for tubulin_pass Skeleton.Density
summary_stats <- tubulin_pass %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p11.5 tub","bj p26.5 tub","bj p31.5 tub",
                                       "bj p36.5 tub","bj p38.5 tub","bj p40.5 tub")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(Skeleton.Density, na.rm = TRUE),
    sd = sd(Skeleton.Density, na.rm = TRUE),
    ymax = max(Skeleton.Density, na.rm = TRUE)
  )

ggplot(mutate(tubulin_pass, condition_clean = fct_relevel(condition_clean, "bj p11.5 tub", "bj p26.5 tub", "bj p31.5 tub", "bj p36.5 tub", "bj p38.5 tub", "bj p40.5 tub")), aes(condition_clean, Skeleton.Density, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Passage", y = "Network density") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('p11.5', 'p26.5', 'p31.5', 'p36.5', 'p38.5', 'p40.5')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(0.05)) +
  geom_text(
    data = summary_stats,
    aes(
      x = condition_clean,
      y = 1.13,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = tubulin_pass, comparisons = list(c("bj p11.5 tub","bj p40.5 tub")), annotations = "ns", map_signif_level = TRUE, y_position = , color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/tubPassDensity.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")

# boxplot for tubulin_pass Branch.Ratio
summary_stats <- tubulin_pass %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p11.5 tub","bj p26.5 tub","bj p31.5 tub",
                                       "bj p36.5 tub","bj p38.5 tub","bj p40.5 tub")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(Branch.Ratio, na.rm = TRUE),
    sd = sd(Branch.Ratio, na.rm = TRUE),
    ymax = max(Branch.Ratio, na.rm = TRUE)
  )

ggplot(mutate(tubulin_pass, condition_clean = fct_relevel(condition_clean, "bj p11.5 tub", "bj p26.5 tub", "bj p31.5 tub", "bj p36.5 tub", "bj p38.5 tub", "bj p40.5 tub")), aes(condition_clean, Branch.Ratio, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Passage", y = "Branch ratio (Branch points/Filaments)") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('p11.5', 'p26.5', 'p31.5', 'p36.5', 'p38.5', 'p40.5')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(0.02)) +
  geom_text(
    data = summary_stats,
    aes(
      x = condition_clean,
      y = 1.68,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = tubulin_pass, comparisons = list(c("bj p11.5 tub","bj p38.5 tub")), annotations = "ns", map_signif_level = TRUE, y_position = 1.7, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/tubPassBranchRatio.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")

# boxplot for tubulin_pass Branch.Point.Density...um.3.
summary_stats <- tubulin_pass %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p11.5 tub","bj p26.5 tub","bj p31.5 tub",
                                       "bj p36.5 tub","bj p38.5 tub","bj p40.5 tub")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(Branch.Point.Density...um.3., na.rm = TRUE),
    sd = sd(Branch.Point.Density...um.3., na.rm = TRUE),
    ymax = max(Branch.Point.Density...um.3., na.rm = TRUE)
  )

ggplot(mutate(tubulin_pass, condition_clean = fct_relevel(condition_clean, "bj p11.5 tub", "bj p26.5 tub", "bj p31.5 tub", "bj p36.5 tub", "bj p38.5 tub", "bj p40.5 tub")), aes(condition_clean, Branch.Point.Density...um.3., fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Passage", y = "Branch point density (Branch points/μm^3)") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('p11.5', 'p26.5', 'p31.5', 'p36.5', 'p38.5', 'p40.5')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(2)) +
  geom_text(
    data = summary_stats,
    aes(
      x = condition_clean,
      y = 26,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = tubulin_pass, comparisons = list(c("bj p11.5 tub","bj p26.5 tub")), annotations = "***", map_signif_level = TRUE, y_position = 28, vjust = 0.7, color = 'black') +
  geom_signif(data = tubulin_pass, comparisons = list(c("bj p11.5 tub","bj p31.5 tub")), annotations = "**", map_signif_level = TRUE, y_position = 29.5, vjust = 0.7, color = 'black') +
  geom_signif(data = tubulin_pass, comparisons = list(c("bj p11.5 tub","bj p38.5 tub")), annotations = "***", map_signif_level = TRUE, y_position = , vjust = 0.7, color = 'black') +
  geom_signif(data = tubulin_pass, comparisons = list(c("bj p11.5 tub","bj p40.5 tub")), annotations = "***", map_signif_level = TRUE, y_position = 32, vjust = 0.7, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/tubPassBranchDensity.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")

# boxplot for tubulin_UV Cell.Volume..um.3.
summary_stats <- tubulin_UV %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p12.5 -uv tub", "bj p12.5 +uv tub")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(Cell.Volume..um.3., na.rm = TRUE),
    sd = sd(Cell.Volume..um.3., na.rm = TRUE),
    ymax = max(Cell.Volume..um.3., na.rm = TRUE)
  )

ggplot(mutate(tubulin_UV, condition_clean = fct_relevel(condition_clean, "bj p12.5 -uv tub", "bj p12.5 +uv tub")), aes(condition_clean, Cell.Volume..um.3., fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Condition", y = "Cell volume (μm^3)") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('p12.5 -UV', 'p12.5 +UV')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(2000)) +
  geom_text(
    data = summary_stats,
    aes(
      x = condition_clean,
      y = 45000,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = tubulin_UV, comparisons = list(c("bj p12.5 -uv tub","bj p12.5 +uv tub")), annotations = "ns", map_signif_level = TRUE, y_position = 47000, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/tubUVCellSize.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")

# boxplot for tubulin_UV Skeleton.Density
summary_stats <- tubulin_UV %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p12.5 -uv tub", "bj p12.5 +uv tub")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(Skeleton.Density, na.rm = TRUE),
    sd = sd(Skeleton.Density, na.rm = TRUE),
    ymax = max(Skeleton.Density, na.rm = TRUE)
  )

ggplot(mutate(tubulin_UV, condition_clean = fct_relevel(condition_clean, "bj p12.5 -uv tub", "bj p12.5 +uv tub")), aes(condition_clean, Skeleton.Density, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Condition", y = "Network density") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('p12.5 -UV', 'p12.5 +UV')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(0.01)) +
  geom_text(
    data = summary_stats,
    aes(
      x = condition_clean,
      y = 0.32,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = tubulin_UV, comparisons = list(c("bj p12.5 -uv tub","bj p12.5 +uv tub")), annotations = "ns", map_signif_level = TRUE, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/tubUVDensity.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")

# boxplot for tubulin_UV Branch.Ratio
summary_stats <- tubulin_UV %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p12.5 -uv tub", "bj p12.5 +uv tub")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(Branch.Ratio, na.rm = TRUE),
    sd = sd(Branch.Ratio, na.rm = TRUE),
    ymax = max(Branch.Ratio, na.rm = TRUE)
  )

ggplot(mutate(tubulin_UV, condition_clean = fct_relevel(condition_clean, "bj p12.5 -uv tub", "bj p12.5 +uv tub")), aes(condition_clean, Branch.Ratio, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Condition", y = "Branch ratio (Branch points/Filaments)") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('p12.5 -UV', 'p12.5 +UV')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(0.01)) +
  geom_text(
    data = summary_stats,
    aes(
      x = condition_clean,
      y = 1.63,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = tubulin_UV, comparisons = list(c("bj p12.5 -uv tub","bj p12.5 +uv tub")), annotations = "ns", map_signif_level = TRUE, y_position = 1.64, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/tubUVBranchRatio.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")

# boxplot for tubulin_UV Branch.Point.Density...um.3.
summary_stats <- tubulin_UV %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p12.5 -uv tub", "bj p12.5 +uv tub")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(Branch.Point.Density...um.3., na.rm = TRUE),
    sd = sd(Branch.Point.Density...um.3., na.rm = TRUE),
    ymax = max(Branch.Point.Density...um.3., na.rm = TRUE)
  )

ggplot(mutate(tubulin_UV, condition_clean = fct_relevel(condition_clean, "bj p12.5 -uv tub", "bj p12.5 +uv tub")), aes(condition_clean, Branch.Point.Density...um.3., fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Condition", y = "Branch point density (Branch points/μm^3)") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('p12.5 -UV', 'p12.5 +UV')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(0.5)) +
  geom_text(
    data = summary_stats,
    aes(
      x = condition_clean,
      y = 14,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = tubulin_UV, comparisons = list(c("bj p12.5 -uv tub","bj p12.5 +uv tub")), annotations = "**", map_signif_level = TRUE, vjust = 0.7, y_position = 14.5, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/tubUVBranchDensity.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")



# boxplot for tubulin_primary Cell.Volume..um.3.
summary_stats <- tubulin_primary %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "primary fibro 18yo tub", "primary fibro 29yo tub", "primary fibro 32yo tub", "primary fibro 59yo tub")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(Cell.Volume..um.3., na.rm = TRUE),
    sd = sd(Cell.Volume..um.3., na.rm = TRUE),
    ymax = max(Cell.Volume..um.3., na.rm = TRUE)
  )

ggplot(mutate(tubulin_primary, condition_clean = fct_relevel(condition_clean, "primary fibro 18yo tub", "primary fibro 29yo tub", "primary fibro 32yo tub", "primary fibro 59yo tub")), aes(condition_clean, Cell.Volume..um.3., fill = condition_clean, color = condition_clean)) + 
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
      y = 100000,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = tubulin_primary, comparisons = list(c("primary fibro 18yo tub", "primary fibro 32yo tub")), annotations = "*", map_signif_level = TRUE, vjust = 0.7, y_position = 105000, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/tubPrimCellSize.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")

# boxplot for tubulin_primary Skeleton.Density
summary_stats <- tubulin_primary %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "primary fibro 18yo tub", "primary fibro 29yo tub", "primary fibro 32yo tub", "primary fibro 59yo tub")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(Skeleton.Density, na.rm = TRUE),
    sd = sd(Skeleton.Density, na.rm = TRUE),
    ymax = max(Skeleton.Density, na.rm = TRUE)
  )

ggplot(mutate(tubulin_primary, condition_clean = fct_relevel(condition_clean, "primary fibro 18yo tub", "primary fibro 29yo tub", "primary fibro 32yo tub", "primary fibro 59yo tub")), aes(condition_clean, Skeleton.Density, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Age", y = "Network density") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('18 y/o', '29 y/o', '32 y/o', '59 y/o')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(0.02)) +
  geom_text(
    data = summary_stats,
    aes(
      x = condition_clean,
      y = 0.55,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = tubulin_primary, comparisons = list(c("primary fibro 18yo tub", "primary fibro 59yo tub")), annotations = "ns", map_signif_level = TRUE, y_position = 0.57, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/tubPrimDensity.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")

# boxplot for tubulin_primary Branch.Ratio
summary_stats <- tubulin_primary %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "primary fibro 18yo tub", "primary fibro 29yo tub", "primary fibro 32yo tub", "primary fibro 59yo tub")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(Branch.Ratio, na.rm = TRUE),
    sd = sd(Branch.Ratio, na.rm = TRUE),
    ymax = max(Branch.Ratio, na.rm = TRUE)
  )

ggplot(mutate(tubulin_primary, condition_clean = fct_relevel(condition_clean, "primary fibro 18yo tub", "primary fibro 29yo tub", "primary fibro 32yo tub", "primary fibro 59yo tub")), aes(condition_clean, Branch.Ratio, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Age", y = "Branch ratio (Branch points/Filaments)") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('18 y/o', '29 y/o', '32 y/o', '59 y/o')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(0.01)) +
  geom_text(
    data = summary_stats,
    aes(
      x = condition_clean,
      y = 1.65,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = tubulin_primary, comparisons = list(c("primary fibro 18yo tub", "primary fibro 59yo tub")), annotations = "ns", map_signif_level = TRUE, y_position = 1.66, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/tubPrimBranchRatio.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")

# boxplot for tubulin_primary Branch.Point.Density...um.3.
summary_stats <- tubulin_primary %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "primary fibro 18yo tub", "primary fibro 29yo tub", "primary fibro 32yo tub", "primary fibro 59yo tub")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(Branch.Point.Density...um.3., na.rm = TRUE),
    sd = sd(Branch.Point.Density...um.3., na.rm = TRUE),
    ymax = max(Branch.Point.Density...um.3., na.rm = TRUE)
  )

ggplot(mutate(tubulin_primary, condition_clean = fct_relevel(condition_clean, "primary fibro 18yo tub", "primary fibro 29yo tub", "primary fibro 32yo tub", "primary fibro 59yo tub")), aes(condition_clean, Branch.Point.Density...um.3., fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Age", y = "Branch point density (Branch points/μm^3)") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('18 y/o', '29 y/o', '32 y/o', '59 y/o')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(1)) +
  geom_text(
    data = summary_stats,
    aes(
      x = condition_clean,
      y = 16.3,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = tubulin_primary, comparisons = list(c("primary fibro 18yo tub", "primary fibro 59yo tub")), annotations = "ns", map_signif_level = TRUE, y_position = 16.8, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/tubPrimBranchDensity.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")
