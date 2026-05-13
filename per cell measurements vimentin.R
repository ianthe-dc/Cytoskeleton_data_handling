library(tidyverse)
library(lme4)
library(car)
library(emmeans)
library(ggprism)
library(forcats)
library(dplyr)
library(ggsignif)
library(DHARMa)

vimentin_df <- read.csv("vimentin cell meas.csv")
vimentin_BJ <- filter(vimentin_df, condition_clean %in% c("bj p8.5 vim", "bj p9.5 etop vim", "bj p11.5 vim", "bj p12.5 -uv vim", "bj p12.5 +uv vim", "bj p26.5 vim", "bj p31.5 vim", "bj p36.5 vim", "bj p38.5 vim", "bj p40.5 vim"))
vimentin_primary <- filter(vimentin_df, condition_clean %in% c("primary fibro 18yo vim", "primary fibro 29yo vim", "primary fibro 32yo vim", "primary fibro 59yo vim"))
vimentin_UV <- filter(vimentin_df, condition_clean %in% c("bj p12.5 -uv vim", "bj p12.5 +uv vim"))
vimentin_etop <- filter(vimentin_df, condition_clean %in% c("bj p8.5 vim", "bj p9.5 etop vim"))
vimentin_pass <- filter(vimentin_df, condition_clean %in% c("bj p11.5 vim", "bj p26.5 vim", "bj p31.5 vim", "bj p36.5 vim", "bj p38.5 vim", "bj p40.5 vim"))


#stats for passage only 
model_pass <- lmer(Branch.Ratio ~ condition_clean + (1|condition_clean:biorep:techrep), vimentin_pass)
summary(model_pass)
Anova(model_pass)
em_pass <- emmeans(model_pass, ~condition_clean)
pairs(em_pass)

simulationOutput <- simulateResiduals(fittedModel = model_pass, plot = T)

kruskal <- kruskal.test(Cell.Volume..um.3. ~ condition_clean, data = vimentin_pass)
kruskal
pairwise.wilcox.test(vimentin_pass$Cell.Volume..um.3., vimentin_pass$condition_clean, p.adjust.method = "BH")


#stats for UV 
model_UV <- lm(Branch.Point.Density...um.3. ~ condition_clean + techrep, vimentin_UV)
summary(model_UV)
Anova(model_UV)
em_UV <- emmeans(model_UV, ~condition_clean)
pairs(em_UV)


simulationOutput <- simulateResiduals(fittedModel = model_UV, plot = T)


#stats for etop only
model_etop <- lmer(Branch.Point.Density...um.3. ~ condition_clean + (1|condition_clean:biorep) + (1|condition_clean:biorep:techrep), vimentin_etop)
summary(model_etop)
Anova(model_etop)
em_etop <- emmeans(model_etop, ~condition_clean)
pairs(em_etop)

simulationOutput <- simulateResiduals(fittedModel = model_etop, plot = T)

kruskal <- kruskal.test(sum_branches ~ condition_clean, data = vimentin_etop)
kruskal


#stats for vimentin primary
model_prim <- lmer(Branch.Point.Density...um.3. ~ condition_clean + (1|condition_clean:biorep:techrep), vimentin_primary)
summary(model_prim)
Anova(model_prim)
em13 <- emmeans(model_prim, ~condition_clean)
pairs(em13)

simulationOutput <- simulateResiduals(fittedModel = model_prim, plot = T)

kruskal <- kruskal.test(Branch.Point.Density...um.3. ~ condition_clean, data = vimentin_primary)
kruskal
pairwise.wilcox.test(vimentin_primary$Branch.Point.Density...um.3., vimentin_primary$condition_clean, p.adjust.method = "BH")


# boxplot for vimentin_pass Cell.Volume..um.3.
summary_stats <- vimentin_pass %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p11.5 vim","bj p26.5 vim","bj p31.5 vim",
                                       "bj p36.5 vim","bj p38.5 vim","bj p40.5 vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(Cell.Volume..um.3., na.rm = TRUE),
    sd = sd(Cell.Volume..um.3., na.rm = TRUE),
    ymax = max(Cell.Volume..um.3., na.rm = TRUE)
  )

ggplot(mutate(vimentin_pass, condition_clean = fct_relevel(condition_clean, "bj p11.5 vim", "bj p26.5 vim", "bj p31.5 vim", "bj p36.5 vim", "bj p38.5 vim", "bj p40.5 vim")), aes(condition_clean, Cell.Volume..um.3., fill = condition_clean, color = condition_clean)) + 
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
      y = 65000,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 2.5
  ) +
  geom_signif(data = vimentin_pass, comparisons = list(c("bj p11.5 vim","bj p31.5 vim")), map_signif_level = TRUE, annotations = "*", y_position = 67000, vjust = 0.7, color = 'black') +
  geom_signif(data = vimentin_pass, comparisons = list(c("bj p11.5 vim","bj p36.5 vim")), map_signif_level = TRUE, annotations = "***", y_position = 70000, vjust = 0.7, color = 'black') +
  geom_signif(data = vimentin_pass, comparisons = list(c("bj p11.5 vim","bj p38.5 vim")), map_signif_level = TRUE, annotations = "***", y_position = 73000, vjust = 0.7, color = 'black') +
  geom_signif(data = vimentin_pass, comparisons = list(c("bj p11.5 vim","bj p40.5 vim")), map_signif_level = TRUE, annotations = "***", y_position = 76000, vjust = 0.7, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/vimPassCellSize.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")

# boxplot for vimentin_pass Skeleton.Density
summary_stats <- vimentin_pass %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p11.5 vim","bj p26.5 vim","bj p31.5 vim",
                                       "bj p36.5 vim","bj p38.5 vim","bj p40.5 vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(Skeleton.Density, na.rm = TRUE),
    sd = sd(Skeleton.Density, na.rm = TRUE),
    ymax = max(Skeleton.Density, na.rm = TRUE)
  )

ggplot(mutate(vimentin_pass, condition_clean = fct_relevel(condition_clean, "bj p11.5 vim", "bj p26.5 vim", "bj p31.5 vim", "bj p36.5 vim", "bj p38.5 vim", "bj p40.5 vim")), aes(condition_clean, Skeleton.Density, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Passage", y = "Network density") +
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
      y = 2.5,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = vimentin_pass, comparisons = list(c("bj p11.5 vim","bj p40.5 vim")), map_signif_level = TRUE, annotations = "ns", y_position = , color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/vimPassDensity.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")

# boxplot for vimentin_pass Branch.Ratio
summary_stats <- vimentin_pass %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p11.5 vim","bj p26.5 vim","bj p31.5 vim",
                                       "bj p36.5 vim","bj p38.5 vim","bj p40.5 vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(Branch.Ratio, na.rm = TRUE),
    sd = sd(Branch.Ratio, na.rm = TRUE),
    ymax = max(Branch.Ratio, na.rm = TRUE)
  )

ggplot(mutate(vimentin_pass, condition_clean = fct_relevel(condition_clean, "bj p11.5 vim", "bj p26.5 vim", "bj p31.5 vim", "bj p36.5 vim", "bj p38.5 vim", "bj p40.5 vim")), aes(condition_clean, Branch.Ratio, fill = condition_clean, color = condition_clean)) + 
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
      y = 1.76,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = vimentin_pass, comparisons = list(c("bj p11.5 vim","bj p26.5 vim")), map_signif_level = TRUE, annotations = "**", y_position = 1.77, vjust = 0.7, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/vimPassBranchRatio.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")

# boxplot for vimentin_pass Branch.Point.Density...um.3.
summary_stats <- vimentin_pass %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p11.5 vim","bj p26.5 vim","bj p31.5 vim",
                                       "bj p36.5 vim","bj p38.5 vim","bj p40.5 vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(Branch.Point.Density...um.3., na.rm = TRUE),
    sd = sd(Branch.Point.Density...um.3., na.rm = TRUE),
    ymax = max(Branch.Point.Density...um.3., na.rm = TRUE)
  )

ggplot(mutate(vimentin_pass, condition_clean = fct_relevel(condition_clean, "bj p11.5 vim", "bj p26.5 vim", "bj p31.5 vim", "bj p36.5 vim", "bj p38.5 vim", "bj p40.5 vim")), aes(condition_clean, Branch.Point.Density...um.3., fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Passage", y = "Branch point density (Branch points/μm^3)") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('p11.5', 'p26.5', 'p31.5', 'p36.5', 'p38.5', 'p40.5')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(0.5)) +
  geom_text(
    data = summary_stats,
    aes(
      x = condition_clean,
      y = 6,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = vimentin_pass, comparisons = list(c("bj p11.5 vim","bj p26.5 vim")), map_signif_level = TRUE, annotations = "*", y_position = , vjust = 0.7, color = 'black') +
  geom_signif(data = vimentin_pass, comparisons = list(c("bj p26.5 vim","bj p40.5 vim")), map_signif_level = TRUE, annotations = "***", y_position = 6.7, vjust = 0.7, color = 'black') 
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/vimPassBranchDensity.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")


# boxplot for vimentin_UV Cell.Volume..um.3.
summary_stats <- vimentin_UV %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p12.5 -uv vim", "bj p12.5 +uv vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(Cell.Volume..um.3., na.rm = TRUE),
    sd = sd(Cell.Volume..um.3., na.rm = TRUE),
    ymax = max(Cell.Volume..um.3., na.rm = TRUE)
  )

ggplot(mutate(vimentin_UV, condition_clean = fct_relevel(condition_clean, "bj p12.5 -uv vim", "bj p12.5 +uv vim")), aes(condition_clean, Cell.Volume..um.3., fill = condition_clean, color = condition_clean)) + 
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
      y = 29000,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = vimentin_UV, comparisons = list(c("bj p12.5 -uv vim","bj p12.5 +uv vim")), annotations = "ns", map_signif_level = TRUE, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/vimUVCellSize.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")

# boxplot for vimentin_UV Skeleton.Density
summary_stats <- vimentin_UV %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p12.5 -uv vim", "bj p12.5 +uv vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(Skeleton.Density, na.rm = TRUE),
    sd = sd(Skeleton.Density, na.rm = TRUE),
    ymax = max(Skeleton.Density, na.rm = TRUE)
  )

ggplot(mutate(vimentin_UV, condition_clean = fct_relevel(condition_clean, "bj p12.5 -uv vim", "bj p12.5 +uv vim")), aes(condition_clean, Skeleton.Density, fill = condition_clean, color = condition_clean)) + 
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
      y = 1.8,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = vimentin_UV, comparisons = list(c("bj p12.5 -uv vim","bj p12.5 +uv vim")), annotations = "ns", map_signif_level = TRUE, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/vimUVDensity.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")

# boxplot for vimentin_UV Branch.Ratio
summary_stats <- vimentin_UV %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p12.5 -uv vim", "bj p12.5 +uv vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(Branch.Ratio, na.rm = TRUE),
    sd = sd(Branch.Ratio, na.rm = TRUE),
    ymax = max(Branch.Ratio, na.rm = TRUE)
  )

ggplot(mutate(vimentin_UV, condition_clean = fct_relevel(condition_clean, "bj p12.5 -uv vim", "bj p12.5 +uv vim")), aes(condition_clean, Branch.Ratio, fill = condition_clean, color = condition_clean)) + 
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
      y = 1.75,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = vimentin_UV, comparisons = list(c("bj p12.5 -uv vim","bj p12.5 +uv vim")), annotations = "ns", map_signif_level = TRUE, y_position = 1.76, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/vimUVBranchRatio.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")

# boxplot for vimentin_UV Branch.Point.Density...um.3.
summary_stats <- vimentin_UV %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p12.5 -uv vim", "bj p12.5 +uv vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(Branch.Point.Density...um.3., na.rm = TRUE),
    sd = sd(Branch.Point.Density...um.3., na.rm = TRUE),
    ymax = max(Branch.Point.Density...um.3., na.rm = TRUE)
  )

ggplot(mutate(vimentin_UV, condition_clean = fct_relevel(condition_clean, "bj p12.5 -uv vim", "bj p12.5 +uv vim")), aes(condition_clean, Branch.Point.Density...um.3., fill = condition_clean, color = condition_clean)) + 
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
      y = 4.2,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = vimentin_UV, comparisons = list(c("bj p12.5 -uv vim","bj p12.5 +uv vim")), annotations = "ns", map_signif_level = TRUE, y_position = 4.4, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/vimUVBranchDensity.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")



# boxplot for vimentin_etop Cell.Volume..um.3.
summary_stats <- vimentin_etop %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p8.5 vim", "bj p9.5 etop vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(Cell.Volume..um.3., na.rm = TRUE),
    sd = sd(Cell.Volume..um.3., na.rm = TRUE),
    ymax = max(Cell.Volume..um.3., na.rm = TRUE)
  )

ggplot(mutate(vimentin_etop, condition_clean = fct_relevel(condition_clean, "bj p8.5 vim", "bj p9.5 etop vim")), aes(condition_clean, Cell.Volume..um.3., fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Condition", y = "Cell volume (μm^3)") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('p8.5 -etop', 'p9.5 +etop')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(2000)) +
  geom_text(
    data = summary_stats,
    aes(
      x = condition_clean,
      y = 53000,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = vimentin_etop, comparisons = list(c("bj p8.5 vim","bj p9.5 etop vim")), annotations = "ns", map_signif_level = TRUE, y_position = 55000, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/vimEtopCellSize.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")

# boxplot for vimentin_etop Skeleton.Density
summary_stats <- vimentin_etop %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p8.5 vim", "bj p9.5 etop vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(Skeleton.Density, na.rm = TRUE),
    sd = sd(Skeleton.Density, na.rm = TRUE),
    ymax = max(Skeleton.Density, na.rm = TRUE)
  )

ggplot(mutate(vimentin_etop, condition_clean = fct_relevel(condition_clean, "bj p8.5 vim", "bj p9.5 etop vim")), aes(condition_clean, Skeleton.Density, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Condition", y = "Network density") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('p8.5 -etop', 'p9.5 +etop')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(0.1)) +
  geom_text(
    data = summary_stats,
    aes(
      x = condition_clean,
      y = 1.78,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = vimentin_etop, comparisons = list(c("bj p8.5 vim","bj p9.5 etop vim")), annotations = "ns", map_signif_level = TRUE, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/vimEtopDensity.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")

# boxplot for vimentin_etop Branch.Ratio
summary_stats <- vimentin_etop %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p8.5 vim", "bj p9.5 etop vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(Branch.Ratio, na.rm = TRUE),
    sd = sd(Branch.Ratio, na.rm = TRUE),
    ymax = max(Branch.Ratio, na.rm = TRUE)
  )

ggplot(mutate(vimentin_etop, condition_clean = fct_relevel(condition_clean, "bj p8.5 vim", "bj p9.5 etop vim")), aes(condition_clean, Branch.Ratio, fill = condition_clean, color = condition_clean)) + 
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
      y = 1.73,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = vimentin_etop, comparisons = list(c("bj p8.5 vim","bj p9.5 etop vim")), annotations = "ns", map_signif_level = TRUE, y_position = 1.74, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/vimEtopBranchRatio.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")

# boxplot for vimentin_etop Branch.Point.Density...um.3.
summary_stats <- vimentin_etop %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p8.5 vim", "bj p9.5 etop vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(Branch.Point.Density...um.3., na.rm = TRUE),
    sd = sd(Branch.Point.Density...um.3., na.rm = TRUE),
    ymax = max(Branch.Point.Density...um.3., na.rm = TRUE)
  )

ggplot(mutate(vimentin_etop, condition_clean = fct_relevel(condition_clean, "bj p8.5 vim", "bj p9.5 etop vim")), aes(condition_clean, Branch.Point.Density...um.3., fill = condition_clean, color = condition_clean)) + 
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
      y = 10.5,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = vimentin_etop, comparisons = list(c("bj p8.5 vim","bj p9.5 etop vim")), annotations = "ns", map_signif_level = TRUE, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/vimEtopBranchDensity.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")


# boxplot for vimentin_primary Cell.Volume..um.3.
summary_stats <- vimentin_primary %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "primary fibro 18yo vim", "primary fibro 29yo vim", "primary fibro 32yo vim", "primary fibro 59yo vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(Cell.Volume..um.3., na.rm = TRUE),
    sd = sd(Cell.Volume..um.3., na.rm = TRUE),
    ymax = max(Cell.Volume..um.3., na.rm = TRUE)
  )

ggplot(mutate(vimentin_primary, condition_clean = fct_relevel(condition_clean, "primary fibro 18yo vim", "primary fibro 29yo vim", "primary fibro 32yo vim", "primary fibro 59yo vim")), aes(condition_clean, Cell.Volume..um.3., fill = condition_clean, color = condition_clean)) + 
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
      y = 81000,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = vimentin_primary, comparisons = list(c("primary fibro 18yo vim", "primary fibro 32yo vim")), annotations = "*", map_signif_level = TRUE, y_position = 84000, vjust = 0.7, color = 'black') +
  geom_signif(data = vimentin_primary, comparisons = list(c("primary fibro 18yo vim", "primary fibro 59yo vim")), annotations = "***", map_signif_level = TRUE, y_position = 87000, vjust = 0.7, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/vimPrimCellSize.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")

# boxplot for vimentin_primary Skeleton.Density
summary_stats <- vimentin_primary %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "primary fibro 18yo vim", "primary fibro 29yo vim", "primary fibro 32yo vim", "primary fibro 59yo vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(Skeleton.Density, na.rm = TRUE),
    sd = sd(Skeleton.Density, na.rm = TRUE),
    ymax = max(Skeleton.Density, na.rm = TRUE)
  )

ggplot(mutate(vimentin_primary, condition_clean = fct_relevel(condition_clean, "primary fibro 18yo vim", "primary fibro 29yo vim", "primary fibro 32yo vim", "primary fibro 59yo vim")), aes(condition_clean, Skeleton.Density, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Age", y = "Network density") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('18 y/o', '29 y/o', '32 y/o', '59 y/o')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(0.05)) +
  geom_text(
    data = summary_stats,
    aes(
      x = condition_clean,
      y = 1.1,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = vimentin_primary, comparisons = list(c("primary fibro 18yo vim", "primary fibro 59yo vim")), annotations = "ns", map_signif_level = TRUE, y_position = 1.15, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/vimPrimDensity.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")

# boxplot for vimentin_primary Branch.Ratio
summary_stats <- vimentin_primary %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "primary fibro 18yo vim", "primary fibro 29yo vim", "primary fibro 32yo vim", "primary fibro 59yo vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(Branch.Ratio, na.rm = TRUE),
    sd = sd(Branch.Ratio, na.rm = TRUE),
    ymax = max(Branch.Ratio, na.rm = TRUE)
  )

ggplot(mutate(vimentin_primary, condition_clean = fct_relevel(condition_clean, "primary fibro 18yo vim", "primary fibro 29yo vim", "primary fibro 32yo vim", "primary fibro 59yo vim")), aes(condition_clean, Branch.Ratio, fill = condition_clean, color = condition_clean)) + 
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
      y = 1.71,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = vimentin_primary, comparisons = list(c("primary fibro 18yo vim", "primary fibro 59yo vim")), annotations = "ns", map_signif_level = TRUE, y_position = 1.72, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/vimPrimBranchRatio.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")

# boxplot for vimentin_primary Branch.Point.Density...um.3.
summary_stats <- vimentin_primary %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "primary fibro 18yo vim", "primary fibro 29yo vim", "primary fibro 32yo vim", "primary fibro 59yo vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(Branch.Point.Density...um.3., na.rm = TRUE),
    sd = sd(Branch.Point.Density...um.3., na.rm = TRUE),
    ymax = max(Branch.Point.Density...um.3., na.rm = TRUE)
  )

ggplot(mutate(vimentin_primary, condition_clean = fct_relevel(condition_clean, "primary fibro 18yo vim", "primary fibro 29yo vim", "primary fibro 32yo vim", "primary fibro 59yo vim")), aes(condition_clean, Branch.Point.Density...um.3., fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Age", y = "Branch point density (Branch points/μm^3)") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('18 y/o', '29 y/o', '32 y/o', '59 y/o')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(0.5)) +
  geom_text(
    data = summary_stats,
    aes(
      x = condition_clean,
      y = 8.5,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = vimentin_primary, comparisons = list(c("primary fibro 18yo vim", "primary fibro 59yo vim")), annotations = "ns", map_signif_level = TRUE, y_position = , color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/vimPrimBranchDensity.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")
