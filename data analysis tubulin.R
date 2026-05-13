# workbook for microtubule analysis

library(tidyverse)
library(lme4)
library(car)
library(emmeans)
library(ggprism)
library(forcats)
library(dplyr)
library(ggsignif)
library(DHARMa)

tubulin_df <- read.csv("tubulin cell level.csv")
tubulin_BJ <- filter(tubulin_df, condition_clean %in% c("bj p11.5 tub", "bj p12.5 -uv tub", "bj p12.5 +uv tub", "bj p26.5 tub", "bj p31.5 tub", "bj p36.5 tub", "bj p38.5 tub", "bj p40.5 tub"))
tubulin_primary <- tubulin_df[167:201,]
tubulin_UV <- filter(tubulin_df, condition_clean %in% c("bj p12.5 -uv tub", "bj p12.5 +uv tub"))
tubulin_pass <- filter(tubulin_df, condition_clean %in% c("bj p11.5 tub", "bj p26.5 tub", "bj p31.5 tub", "bj p36.5 tub", "bj p38.5 tub", "bj p40.5 tub"))

#stats for all tubulin BJ datasets
model <- lmer(filament_count ~ condition_clean + (1|condition_clean:biorep) +(1|condition_clean:biorep:techrep), tubulin_BJ)
summary(model)
Anova(model)
em1 <- emmeans(model, ~condition_clean)
pairs(em1)

#stats for passage only sum_length_2d
model_pass <- lm(sum_length_2d ~ condition_clean + techrep, tubulin_pass)
summary(model_pass)
Anova(model_pass)
em_pass <- emmeans(model_pass, ~condition_clean)
pairs(em_pass)

simulationOutput <- simulateResiduals(fittedModel = model_pass, plot = T)

kruskal <- kruskal.test(sum_length_2d ~ condition_clean, data = tubulin_pass)
kruskal
pairwise.wilcox.test(tubulin_pass$sum_length_2d, tubulin_pass$condition_clean, p.adjust.method = "BH")

#stats for passage only mean_length_2d
model2_pass <- lmer(mean_length_2d ~ condition_clean + (1|condition_clean:biorep) + (1|condition_clean:biorep:techrep), tubulin_pass)
summary(model2_pass)
Anova(model2_pass)
em_pass2 <- emmeans(model2_pass, ~condition_clean)
pairs(em_pass2)

simulationOutput <- simulateResiduals(fittedModel = model2_pass, plot = T)

#stats for passage only sum_length_3d
model3_pass <- lm(sum_length_3d ~ condition_clean + techrep, tubulin_pass)
summary(model3_pass)
Anova(model3_pass)
em_pass3 <- emmeans(model3_pass, ~condition_clean)
pairs(em_pass3)

simulationOutput <- simulateResiduals(fittedModel = model3_pass, plot = T)

kruskal <- kruskal.test(sum_length_3d ~ condition_clean, data = tubulin_pass)
kruskal
pairwise.wilcox.test(tubulin_pass$sum_length_3d, tubulin_pass$condition_clean, p.adjust.method = "BH")

#stats for passage only mean_length_3d
model4_pass <- lmer(mean_length_3d ~ condition_clean + (1|condition_clean:biorep:techrep), tubulin_pass)
summary(model4_pass)
Anova(model4_pass)
em_pass4 <- emmeans(model4_pass, ~condition_clean)
pairs(em_pass4)

simulationOutput <- simulateResiduals(fittedModel = model4_pass, plot = T)

#stats for passage only mean_avg_angle
model5_pass <- lmer(mean_avg_angle ~ condition_clean + (1|condition_clean:biorep) + (1|condition_clean:biorep:techrep), tubulin_pass)
summary(model5_pass)
Anova(model5_pass)
em_pass5 <- emmeans(model5_pass, ~condition_clean)
pairs(em_pass5)

simulationOutput <- simulateResiduals(fittedModel = model5_pass, plot = T)

#stats for passage only mean_acc_curv
model6_pass <- lmer(mean_acc_curv ~ condition_clean + (1|condition_clean:biorep) + (1|condition_clean:biorep:techrep), tubulin_pass)
summary(model6_pass)
Anova(model6_pass)
em_pass6 <- emmeans(model6_pass, ~condition_clean)
pairs(em_pass6)

simulationOutput <- simulateResiduals(fittedModel = model6_pass, plot = T)

#stats for passage only mean_net_curv
model6_pass <- lmer(mean_net_curv ~ condition_clean + (1|condition_clean:biorep) + (1|condition_clean:biorep:techrep), tubulin_pass)
summary(model6_pass)
Anova(model6_pass)
em_pass6 <- emmeans(model6_pass, ~condition_clean)
pairs(em_pass6)

simulationOutput <- simulateResiduals(fittedModel = model6_pass, plot = T)

#stats for passage only filament_count
model7_pass <- lm(filament_count ~ condition_clean + techrep + biorep, tubulin_pass)
summary(model7_pass)
Anova(model7_pass)
em_pass7 <- emmeans(model7_pass, ~condition_clean)
pairs(em_pass7)

simulationOutput <- simulateResiduals(fittedModel = model7_pass, plot = T)

kruskal <- kruskal.test(filament_count ~ condition_clean, data = tubulin_pass)
kruskal
pairwise.wilcox.test(tubulin_pass$filament_count, tubulin_pass$condition_clean, p.adjust.method = "BH")

#stats for passage only sum_branches
model8_pass <- lm(sum_branches ~ condition_clean + techrep + biorep, tubulin_pass)
summary(model8_pass)
Anova(model8_pass)
em_pass8 <- emmeans(model8_pass, ~condition_clean)
pairs(em_pass8)

simulationOutput <- simulateResiduals(fittedModel = model8_pass, plot = T)

kruskal <- kruskal.test(sum_branches ~ condition_clean, data = tubulin_pass)
kruskal
pairwise.wilcox.test(tubulin_pass$sum_branches, tubulin_pass$condition_clean, p.adjust.method = "BH")

#stats for passage only mean_branches
model9_pass <- lmer(mean_branches ~ condition_clean + (1|condition_clean:biorep) + (1|condition_clean:biorep:techrep), tubulin_pass)
summary(model9_pass)
Anova(model9_pass)
em_pass9 <- emmeans(model9_pass, ~condition_clean)
pairs(em_pass9)

simulationOutput <- simulateResiduals(fittedModel = model9_pass, plot = T)

kruskal <- kruskal.test(mean_branches ~ condition_clean, data = tubulin_pass)
kruskal
pairwise.wilcox.test(tubulin_pass$mean_branches, tubulin_pass$condition_clean, p.adjust.method = "BH")

#stats for passage only mean_deviation
model10_pass <- lmer(mean_deviation ~ condition_clean + (1|condition_clean:biorep), tubulin_pass)
summary(model10_pass)
Anova(model10_pass)
em_pass10 <- emmeans(model10_pass, ~condition_clean)
pairs(em_pass10)

simulationOutput <- simulateResiduals(fittedModel = model10_pass, plot = T)

#stats for tubulin_UV
model_UV <- lmer(mean_avg_angle ~ condition_clean + (1|condition_clean:biorep:techrep), tubulin_UV)
summary(model_UV)
Anova(model_UV)
em_UV <- emmeans(model_UV, ~condition_clean)
pairs(em_UV)

simulationOutput <- simulateResiduals(fittedModel = model_UV, plot = T)

#stats for tubulin_primary
model_prim <- lm(sum_length_3d ~ condition_clean, tubulin_primary)
summary(model_prim)
Anova(model_prim)
em_prim <- emmeans(model_prim, ~condition_clean)
pairs(em_prim)

simulationOutput <- simulateResiduals(fittedModel = model_prim, plot = T)

# boxplot for all tubulin 
ggplot(mutate(tubulin_BJ, condition_clean = fct_relevel(condition_clean, "bj p8.5 tub", "bj p9.5 etop tub", "bj p11.5 tub", "bj p12.5 -uv tub", "bj p12.5 +uv tub", "bj p26.5 tub", "bj p31.5 tub", "bj p36.5 tub", "bj p38.5 tub", "bj p40.5 tub")), aes(condition_clean, sum_length_2d, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(size = 0.5, position = position_jitterdodge()) +
  labs(x = "Condition", y = "Total 2D filament length per cell (μm)") +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('p11.5', 'p12.5 -UV', 'p12.5 +UV', 'p26.5', 'p31.5', 'p36.5', 'p38.5', 'p40.5')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(5000))
ggsave("tubBJsum2D.tiff", device = "tiff", dpi = 300, width = 177, height = 150, units = "mm")


# boxplot for all tubulin mean_length_2D
ggplot(tubulin_BJ, aes(condition_clean, mean_length_2d, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(size = 0.5, position = position_jitterdodge()) +
  labs(x = "Condition", y = "Mean 2D filament length per cell (μm)") +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('p11.5', 'p12.5 -UV', 'p12.5 +UV', 'p26.5', 'p31.5', 'p36.5', 'p38.5', 'p40.5')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(0.01))
ggsave("tubBJmean2D.tiff", device = "tiff", dpi = 300, width = 177, height = 150, units = "mm")

# boxplot for all tubulin sum_length_3D
ggplot(mutate(tubulin_BJ, condition_clean = fct_relevel(condition_clean, "bj p8.5 tub", "bj p9.5 etop tub", "bj p11.5 tub", "bj p12.5 -uv tub", "bj p12.5 +uv tub", "bj p26.5 tub", "bj p31.5 tub", "bj p36.5 tub", "bj p38.5 tub", "bj p40.5 tub")), aes(condition_clean, sum_length_3d, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(size = 0.5, position = position_jitterdodge()) +
  labs(x = "Condition", y = "Total 3D filament length per cell (μm)") +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('p11.5', 'p12.5 -UV', 'p12.5 +UV', 'p26.5', 'p31.5', 'p36.5', 'p38.5', 'p40.5')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(200000))
ggsave("tubBJsum3D.tiff", device = "tiff", dpi = 300, width = 177, height = 150, units = "mm")

# boxplot for all tubulin mean_length_3D
ggplot(mutate(tubulin_BJ, condition_clean = fct_relevel(condition_clean, "bj p8.5 tub", "bj p9.5 etop tub", "bj p11.5 tub", "bj p12.5 -uv tub", "bj p12.5 +uv tub", "bj p26.5 tub", "bj p31.5 tub", "bj p36.5 tub", "bj p38.5 tub", "bj p40.5 tub")), aes(condition_clean, mean_length_3d, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(size = 0.5, position = position_jitterdodge()) +
  labs(x = "Condition", y = "Mean 3D filament length per cell (μm)") +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('p11.5', 'p12.5 -UV', 'p12.5 +UV', 'p26.5', 'p31.5', 'p36.5', 'p38.5', 'p40.5')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(0.5))
ggsave("tubBJmean3D.tiff", device = "tiff", dpi = 300, width = 177, height = 150, units = "mm")

# boxplot for tubulin_pass sum_length_2d
summary_stats <- tubulin_pass %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p11.5 tub","bj p26.5 tub","bj p31.5 tub",
                                       "bj p36.5 tub","bj p38.5 tub","bj p40.5 tub")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(sum_length_2d, na.rm = TRUE),
    sd = sd(sum_length_2d, na.rm = TRUE),
    ymax = max(sum_length_2d, na.rm = TRUE)
  )

ggplot(mutate(tubulin_pass, condition_clean = fct_relevel(condition_clean, "bj p11.5 tub", "bj p26.5 tub", "bj p31.5 tub", "bj p36.5 tub", "bj p38.5 tub", "bj p40.5 tub")), aes(condition_clean, sum_length_2d, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Passage", y = "Total 2D filament length per cell (μm)") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('p11.5', 'p26.5', 'p31.5', 'p36.5', 'p38.5', 'p40.5')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(10000)) +
  geom_text(
    data = summary_stats,
    aes(
      x = condition_clean,
      y = 290000,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 2.4
  ) +
  geom_signif(data = tubulin_pass, comparisons = list(c("bj p11.5 tub","bj p26.5 tub")), annotations = "***", map_signif_level = TRUE, y_position = 3.1e+05, vjust = 0.7, color = 'black') +
  geom_signif(data = tubulin_pass, comparisons = list(c("bj p11.5 tub","bj p31.5 tub")), annotations = "***", map_signif_level = TRUE, y_position = 3.3e+05, vjust = 0.7, color = 'black') +
  geom_signif(data = tubulin_pass, comparisons = list(c("bj p11.5 tub","bj p36.5 tub")), annotations = "***", map_signif_level = TRUE, y_position = 3.5e+05, vjust = 0.7, color = 'black') +
  geom_signif(data = tubulin_pass, comparisons = list(c("bj p11.5 tub","bj p38.5 tub")), annotations = "***", map_signif_level = TRUE, y_position = 3.7e+05, vjust = 0.7, color = 'black') +
  geom_signif(data = tubulin_pass, comparisons = list(c("bj p11.5 tub","bj p40.5 tub")), annotations = "***", map_signif_level = TRUE, y_position = 3.9e+05, vjust = 0.7, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/tubPassSum2D.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")

# boxplot for tubulin_pass mean_length_2d
summary_stats <- tubulin_pass %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p11.5 tub","bj p26.5 tub","bj p31.5 tub",
                                       "bj p36.5 tub","bj p38.5 tub","bj p40.5 tub")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_length_2d, na.rm = TRUE),
    sd = sd(mean_length_2d, na.rm = TRUE),
    ymax = max(mean_length_2d, na.rm = TRUE)
  )

ggplot(mutate(tubulin_pass, condition_clean = fct_relevel(condition_clean, "bj p11.5 tub", "bj p26.5 tub", "bj p31.5 tub", "bj p36.5 tub", "bj p38.5 tub", "bj p40.5 tub")), aes(condition_clean, mean_length_2d, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Passage", y = "Mean 2D filament length per cell (μm)") +
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
      y = 0.93,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = tubulin_pass, comparisons = list(c("bj p11.5 tub","bj p40.5 tub")), annotations = "ns", map_signif_level = TRUE, y_position = 0.94, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/tubPassMean2D.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")

# boxplot for tubulin_pass sum_length_3d
summary_stats <- tubulin_pass %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p11.5 tub","bj p26.5 tub","bj p31.5 tub",
                                       "bj p36.5 tub","bj p38.5 tub","bj p40.5 tub")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(sum_length_3d, na.rm = TRUE),
    sd = sd(sum_length_3d, na.rm = TRUE),
    ymax = max(sum_length_3d, na.rm = TRUE)
  )

ggplot(mutate(tubulin_pass, condition_clean = fct_relevel(condition_clean, "bj p11.5 tub", "bj p26.5 tub", "bj p31.5 tub", "bj p36.5 tub", "bj p38.5 tub", "bj p40.5 tub")), aes(condition_clean, sum_length_3d, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Passage", y = "Total 3D filament length per cell (μm)") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('p11.5', 'p26.5', 'p31.5', 'p36.5', 'p38.5', 'p40.5')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(250000)) +
  geom_text(
    data = summary_stats,
    aes(
      x = condition_clean,
      y = 5000000,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 2.3
  ) +
  geom_signif(data = tubulin_pass, comparisons = list(c("bj p11.5 tub","bj p26.5 tub")), annotations = "***", map_signif_level = TRUE, y_position = 5.4e+06, vjust = 0.7, color = 'black') +
  geom_signif(data = tubulin_pass, comparisons = list(c("bj p11.5 tub","bj p31.5 tub")), annotations = "***", map_signif_level = TRUE, y_position = 5.7e+06, vjust = 0.7, color = 'black') +
  geom_signif(data = tubulin_pass, comparisons = list(c("bj p11.5 tub","bj p36.5 tub")), annotations = "**", map_signif_level = TRUE, y_position = 6e+06, vjust = 0.7, color = 'black') +
  geom_signif(data = tubulin_pass, comparisons = list(c("bj p11.5 tub","bj p38.5 tub")), annotations = "***", map_signif_level = TRUE, y_position = 6.3e+06, vjust = 0.7, color = 'black') +
  geom_signif(data = tubulin_pass, comparisons = list(c("bj p11.5 tub","bj p40.5 tub")), annotations = "***", map_signif_level = TRUE, y_position = 6.6e+06, vjust = 0.7, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/tubPassSum3D.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")

# boxplot for tubulin_pass mean_length_3d
summary_stats <- tubulin_pass %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p11.5 tub","bj p26.5 tub","bj p31.5 tub",
                                       "bj p36.5 tub","bj p38.5 tub","bj p40.5 tub")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_length_3d, na.rm = TRUE),
    sd = sd(mean_length_3d, na.rm = TRUE),
    ymax = max(mean_length_3d, na.rm = TRUE)
  )

ggplot(mutate(tubulin_pass, condition_clean = fct_relevel(condition_clean, "bj p11.5 tub", "bj p26.5 tub", "bj p31.5 tub", "bj p36.5 tub", "bj p38.5 tub", "bj p40.5 tub")), aes(condition_clean, mean_length_3d, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Passage", y = "Mean 3D filament length per cell (μm)") +
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
      y = 12.5,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = tubulin_pass, comparisons = list(c("bj p11.5 tub","bj p26.5 tub")), annotations = "**", map_signif_level = TRUE, y_position = 12.7, vjust = 0.7, color = 'black') +
  geom_signif(data = tubulin_pass, comparisons = list(c("bj p11.5 tub","bj p38.5 tub")), annotations = "***", map_signif_level = TRUE, y_position = 13.1, vjust = 0.7, color = 'black') 
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/tubPassMean3D.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")

# boxplot for tubulin_pass mean_avg_angle
summary_stats <- tubulin_pass %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p11.5 tub","bj p26.5 tub","bj p31.5 tub",
                                       "bj p36.5 tub","bj p38.5 tub","bj p40.5 tub")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_avg_angle, na.rm = TRUE),
    sd = sd(mean_avg_angle, na.rm = TRUE),
    ymax = max(mean_avg_angle, na.rm = TRUE)
  )

ggplot(mutate(tubulin_pass, condition_clean = fct_relevel(condition_clean, "bj p11.5 tub", "bj p26.5 tub", "bj p31.5 tub", "bj p36.5 tub", "bj p38.5 tub", "bj p40.5 tub")), aes(condition_clean, mean_avg_angle, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Passage", y = "Mean average filament angle (deg)") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('p11.5', 'p26.5', 'p31.5', 'p36.5', 'p38.5', 'p40.5')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(1)) +
  geom_text(
    data = summary_stats,
    aes(
      x = condition_clean,
      y = 90,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = tubulin_pass, comparisons = list(c("bj p11.5 tub","bj p40.5 tub")), annotations = "ns", map_signif_level = TRUE, y_position = 91, color = 'black') 
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/tubPassMnAvgAng.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")

# boxplot for tubulin_pass mean_acc_curv
summary_stats <- tubulin_pass %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p11.5 tub","bj p26.5 tub","bj p31.5 tub",
                                       "bj p36.5 tub","bj p38.5 tub","bj p40.5 tub")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_acc_curv, na.rm = TRUE),
    sd = sd(mean_acc_curv, na.rm = TRUE),
    ymax = max(mean_acc_curv, na.rm = TRUE)
  )

ggplot(mutate(tubulin_pass, condition_clean = fct_relevel(condition_clean, "bj p11.5 tub", "bj p26.5 tub", "bj p31.5 tub", "bj p36.5 tub", "bj p38.5 tub", "bj p40.5 tub")), aes(condition_clean, mean_acc_curv, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Passage", y = "Mean accumulated filament curvature /cell") +
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
      y = 4.4,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = tubulin_pass, comparisons = list(c("bj p11.5 tub","bj p40.5 tub")), annotations = "ns", map_signif_level = TRUE, y_position = 4.42, color = 'black') 
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/tubPassMnAccCurv.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")

# boxplot for tubulin_pass mean_net_curv
summary_stats <- tubulin_pass %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p11.5 tub","bj p26.5 tub","bj p31.5 tub",
                                       "bj p36.5 tub","bj p38.5 tub","bj p40.5 tub")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_net_curv, na.rm = TRUE),
    sd = sd(mean_net_curv, na.rm = TRUE),
    ymax = max(mean_net_curv, na.rm = TRUE)
  )

ggplot(mutate(tubulin_pass, condition_clean = fct_relevel(condition_clean, "bj p11.5 tub", "bj p26.5 tub", "bj p31.5 tub", "bj p36.5 tub", "bj p38.5 tub", "bj p40.5 tub")), aes(condition_clean, mean_net_curv, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Passage", y = "Mean end-to-end filament curvature /cell") +
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
      y = 1.57,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = tubulin_pass, comparisons = list(c("bj p11.5 tub","bj p40.5 tub")), annotations = "ns", map_signif_level = TRUE, y_position = 1.59, color = 'black') 
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/tubPassMnNetCurv.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")

# boxplot for tubulin_pass filament_count
summary_stats <- tubulin_pass %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p11.5 tub","bj p26.5 tub","bj p31.5 tub",
                                       "bj p36.5 tub","bj p38.5 tub","bj p40.5 tub")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(filament_count, na.rm = TRUE),
    sd = sd(filament_count, na.rm = TRUE),
    ymax = max(filament_count, na.rm = TRUE)
  )

ggplot(mutate(tubulin_pass, condition_clean = fct_relevel(condition_clean, "bj p11.5 tub", "bj p26.5 tub", "bj p31.5 tub", "bj p36.5 tub", "bj p38.5 tub", "bj p40.5 tub")), aes(condition_clean, filament_count, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Passage", y = "Filament count /cell") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('p11.5', 'p26.5', 'p31.5', 'p36.5', 'p38.5', 'p40.5')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(20000)) +
  geom_text(
    data = summary_stats,
    aes(
      x = condition_clean,
      y = 400000,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 2.3
  ) +
  geom_signif(data = tubulin_pass, comparisons = list(c("bj p11.5 tub","bj p26.5 tub")), annotations = "***", map_signif_level = TRUE, y_position = 4.4e+05, vjust = 0.7, color = 'black') +
  geom_signif(data = tubulin_pass, comparisons = list(c("bj p11.5 tub","bj p31.5 tub")), annotations = "***", map_signif_level = TRUE, y_position = 4.7e+05, vjust = 0.7, color = 'black') +
  geom_signif(data = tubulin_pass, comparisons = list(c("bj p11.5 tub","bj p36.5 tub")), annotations = "**", map_signif_level = TRUE, y_position = 5e+05, vjust = 0.7, color = 'black') +
  geom_signif(data = tubulin_pass, comparisons = list(c("bj p11.5 tub","bj p38.5 tub")), annotations = "***", map_signif_level = TRUE, y_position = 5.3e+05, vjust = 0.7, color = 'black') +
  geom_signif(data = tubulin_pass, comparisons = list(c("bj p11.5 tub","bj p40.5 tub")), annotations = "***", map_signif_level = TRUE, y_position = 5.6e+05, vjust = 0.7, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/tubPassFilCount.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")

# boxplot for tubulin_pass sum_branches
summary_stats <- tubulin_pass %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p11.5 tub","bj p26.5 tub","bj p31.5 tub",
                                       "bj p36.5 tub","bj p38.5 tub","bj p40.5 tub")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(sum_branches, na.rm = TRUE),
    sd = sd(sum_branches, na.rm = TRUE),
    ymax = max(sum_branches, na.rm = TRUE)
  )

ggplot(mutate(tubulin_pass, condition_clean = fct_relevel(condition_clean, "bj p11.5 tub", "bj p26.5 tub", "bj p31.5 tub", "bj p36.5 tub", "bj p38.5 tub", "bj p40.5 tub")), aes(condition_clean, sum_branches, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Passage", y = "Total filament branches /cell") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('p11.5', 'p26.5', 'p31.5', 'p36.5', 'p38.5', 'p40.5')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(10000)) +
  geom_text(
    data = summary_stats,
    aes(
      x = condition_clean,
      y = 370000,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 2.4
  ) +
  geom_signif(data = tubulin_pass, comparisons = list(c("bj p11.5 tub","bj p26.5 tub")), annotations = "**", map_signif_level = TRUE, y_position = 400000, vjust = 0.7, color = 'black') +
  geom_signif(data = tubulin_pass, comparisons = list(c("bj p11.5 tub","bj p31.5 tub")), annotations = "**", map_signif_level = TRUE, y_position = 420000, vjust = 0.7, color = 'black') +
  geom_signif(data = tubulin_pass, comparisons = list(c("bj p11.5 tub","bj p36.5 tub")), annotations = "**", map_signif_level = TRUE, y_position = 440000, vjust = 0.7, color = 'black') +
  geom_signif(data = tubulin_pass, comparisons = list(c("bj p11.5 tub","bj p38.5 tub")), annotations = "***", map_signif_level = TRUE, y_position = 460000, vjust = 0.7, color = 'black') +
  geom_signif(data = tubulin_pass, comparisons = list(c("bj p11.5 tub","bj p40.5 tub")), annotations = "***", map_signif_level = TRUE, y_position = 480000, vjust = 0.7, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/tubPassBranchCount.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")

# boxplot for tubulin_pass mean_branches
summary_stats <- tubulin_pass %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p11.5 tub","bj p26.5 tub","bj p31.5 tub",
                                       "bj p36.5 tub","bj p38.5 tub","bj p40.5 tub")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_branches, na.rm = TRUE),
    sd = sd(mean_branches, na.rm = TRUE),
    ymax = max(mean_branches, na.rm = TRUE)
  )

ggplot(mutate(tubulin_pass, condition_clean = fct_relevel(condition_clean, "bj p11.5 tub", "bj p26.5 tub", "bj p31.5 tub", "bj p36.5 tub", "bj p38.5 tub", "bj p40.5 tub")), aes(condition_clean, mean_branches, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Passage", y = "Mean filament branches /cell") +
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
      y = 0.84,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = tubulin_pass, comparisons = list(c("bj p11.5 tub","bj p40.5 tub")), annotations = "ns", map_signif_level = TRUE, y_position = 0.86, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/tubPassMeanBranch.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")

# boxplot for tubulin_pass mean_deviation
summary_stats <- tubulin_pass %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p11.5 tub","bj p26.5 tub","bj p31.5 tub",
                                       "bj p36.5 tub","bj p38.5 tub","bj p40.5 tub")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_deviation, na.rm = TRUE),
    sd = sd(mean_deviation, na.rm = TRUE),
    ymax = max(mean_deviation, na.rm = TRUE)
  )

ggplot(mutate(tubulin_pass, condition_clean = fct_relevel(condition_clean, "bj p11.5 tub", "bj p26.5 tub", "bj p31.5 tub", "bj p36.5 tub", "bj p38.5 tub", "bj p40.5 tub")), aes(condition_clean, mean_deviation, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Passage", y = "Mean filament deviation (μm)") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('p11.5', 'p26.5', 'p31.5', 'p36.5', 'p38.5', 'p40.5')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(0.1)) +
  geom_text(
    data = summary_stats,
    aes(
      x = condition_clean,
      y = 2.3,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = tubulin_pass, comparisons = list(c("bj p11.5 tub","bj p26.5 tub")), annotations = "*", map_signif_level = TRUE, y_position = 2.4, vjust = 0.7, color = 'black') +
  geom_signif(data = tubulin_pass, comparisons = list(c("bj p11.5 tub","bj p38.5 tub")), annotations = "*", map_signif_level = TRUE, y_position = 2.5, vjust = 0.7, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/tubPassMeanDev.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")

# dotplot for tubulin_pass sum_length_3d vs sum_branches
ggplot(data = tubulin_pass, aes(x = sum_length_3d, y = sum_branches, fill = condition_clean, color = condition_clean)) + 
  geom_point(size = 0.8) + 
  labs(x = "Total 3D filament length (μm)", y = "Total branches /cell") +
  theme_bw() +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE), x = guide_axis(minor.ticks = TRUE)) +
  scale_y_continuous(minor_breaks = scales::breaks_width(5000)) 
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/tubPassSum2D.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")

# boxplot for tubulin_UV sum_length_2d
summary_stats <- tubulin_UV %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p12.5 -uv tub", "bj p12.5 +uv tub")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(sum_length_2d, na.rm = TRUE),
    sd = sd(sum_length_2d, na.rm = TRUE),
    ymax = max(sum_length_2d, na.rm = TRUE)
  )

ggplot(mutate(tubulin_UV, condition_clean = fct_relevel(condition_clean, "bj p12.5 -uv tub", "bj p12.5 +uv tub")), aes(condition_clean, sum_length_2d, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Condition", y = "Total 2D filament length per cell (μm)") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('p12.5 -UV', 'p12.5 +UV')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(10000)) +
  geom_text(
    data = summary_stats,
    aes(
      x = condition_clean,
      y = 170000,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = tubulin_UV, comparisons = list(c("bj p12.5 -uv tub","bj p12.5 +uv tub")), annotations = "ns", map_signif_level = TRUE, y_position = 180000, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/tubUVsum2d.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")

# boxplot for tubulin_UV mean_length_2d
summary_stats <- tubulin_UV %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p12.5 -uv tub", "bj p12.5 +uv tub")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_length_2d, na.rm = TRUE),
    sd = sd(mean_length_2d, na.rm = TRUE),
    ymax = max(mean_length_2d, na.rm = TRUE)
  )

ggplot(mutate(tubulin_UV, condition_clean = fct_relevel(condition_clean, "bj p12.5 -uv tub", "bj p12.5 +uv tub")), aes(condition_clean, mean_length_2d, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Condition", y = "Mean 2D filament length per cell (μm)") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('p12.5 -UV', 'p12.5 +UV')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(0.002)) +
  geom_text(
    data = summary_stats,
    aes(
      x = condition_clean,
      y = 0.75,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = tubulin_UV, comparisons = list(c("bj p12.5 -uv tub","bj p12.5 +uv tub")), annotations = "ns", map_signif_level = TRUE, y_position = 0.753, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/tubUVmean2d.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")

# boxplot for tubulin_UV sum_length_3d
summary_stats <- tubulin_UV %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p12.5 -uv tub", "bj p12.5 +uv tub")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(sum_length_3d, na.rm = TRUE),
    sd = sd(sum_length_3d, na.rm = TRUE),
    ymax = max(sum_length_3d, na.rm = TRUE)
  )

ggplot(mutate(tubulin_UV, condition_clean = fct_relevel(condition_clean, "bj p12.5 -uv tub", "bj p12.5 +uv tub")), aes(condition_clean, sum_length_3d, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Condition", y = "Total 3D filament length per cell (μm)") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('p12.5 -UV', 'p12.5 +UV')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(100000)) +
  geom_text(
    data = summary_stats,
    aes(
      x = condition_clean,
      y = 1800000,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = tubulin_UV, comparisons = list(c("bj p12.5 -uv tub","bj p12.5 +uv tub")), annotations = "ns", map_signif_level = TRUE, y_position = 1900000, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/tubUVsum3d.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")

# boxplot for tubulin_UV mean_length_3d
summary_stats <- tubulin_UV %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p12.5 -uv tub", "bj p12.5 +uv tub")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_length_3d, na.rm = TRUE),
    sd = sd(mean_length_3d, na.rm = TRUE),
    ymax = max(mean_length_3d, na.rm = TRUE)
  )

ggplot(mutate(tubulin_UV, condition_clean = fct_relevel(condition_clean, "bj p12.5 -uv tub", "bj p12.5 +uv tub")), aes(condition_clean, mean_length_3d, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Condition", y = "Mean 3D filament length per cell (μm)") +
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
      y = 11.3,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = tubulin_UV, comparisons = list(c("bj p12.5 -uv tub","bj p12.5 +uv tub")), annotations = "ns", map_signif_level = TRUE, y_position = 11.6, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/tubUVmean3d.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")

# boxplot for tubulin_UV mean_deviation
summary_stats <- tubulin_UV %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p12.5 -uv tub", "bj p12.5 +uv tub")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_deviation, na.rm = TRUE),
    sd = sd(mean_deviation, na.rm = TRUE),
    ymax = max(mean_deviation, na.rm = TRUE)
  )

ggplot(mutate(tubulin_UV, condition_clean = fct_relevel(condition_clean, "bj p12.5 -uv tub", "bj p12.5 +uv tub")), aes(condition_clean, mean_deviation, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Condition", y = "Mean filament deviation (μm)") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('p12.5 -UV', 'p12.5 +UV')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(0.05)) +
  geom_text(
    data = summary_stats,
    aes(
      x = condition_clean,
      y = 2.05,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = tubulin_UV, comparisons = list(c("bj p12.5 -uv tub","bj p12.5 +uv tub")), annotations = "ns", map_signif_level = TRUE, y_position = 2.1, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/tubUVmeanDev.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")

# boxplot for tubulin_UV mean_avg_angle
summary_stats <- tubulin_UV %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p12.5 -uv tub", "bj p12.5 +uv tub")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_avg_angle, na.rm = TRUE),
    sd = sd(mean_avg_angle, na.rm = TRUE),
    ymax = max(mean_avg_angle, na.rm = TRUE)
  )

ggplot(mutate(tubulin_UV, condition_clean = fct_relevel(condition_clean, "bj p12.5 -uv tub", "bj p12.5 +uv tub")), aes(condition_clean, mean_avg_angle, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Condition", y = "Mean average filament angle (deg)") +
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
      y = 84.6,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = tubulin_UV, comparisons = list(c("bj p12.5 -uv tub","bj p12.5 +uv tub")), annotations = "ns", map_signif_level = TRUE, y_position = 85, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/tubUVmeanAvgAng.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")

# boxplot for tubulin_UV mean_acc_curv
summary_stats <- tubulin_UV %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p12.5 -uv tub", "bj p12.5 +uv tub")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_acc_curv, na.rm = TRUE),
    sd = sd(mean_acc_curv, na.rm = TRUE),
    ymax = max(mean_acc_curv, na.rm = TRUE)
  )

ggplot(mutate(tubulin_UV, condition_clean = fct_relevel(condition_clean, "bj p12.5 -uv tub", "bj p12.5 +uv tub")), aes(condition_clean, mean_acc_curv, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Condition", y = "Mean accumulated filament curvature /cell") +
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
      y = 4.2,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = tubulin_UV, comparisons = list(c("bj p12.5 -uv tub","bj p12.5 +uv tub")), annotations = "ns", map_signif_level = TRUE, y_position = 4.21, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/tubUVmeanAccCurv.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")

# boxplot for tubulin_UV mean_net_curv
summary_stats <- tubulin_UV %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p12.5 -uv tub", "bj p12.5 +uv tub")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_net_curv, na.rm = TRUE),
    sd = sd(mean_net_curv, na.rm = TRUE),
    ymax = max(mean_net_curv, na.rm = TRUE)
  )

ggplot(mutate(tubulin_UV, condition_clean = fct_relevel(condition_clean, "bj p12.5 -uv tub", "bj p12.5 +uv tub")), aes(condition_clean, mean_net_curv, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Condition", y = "Mean end-to-end filament curvature /cell") +
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
      y = 1.46,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = tubulin_UV, comparisons = list(c("bj p12.5 -uv tub","bj p12.5 +uv tub")), annotations = "ns", map_signif_level = TRUE, y_position = 1.47, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/tubUVmeanNetCurv.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")

# boxplot for tubulin_UV filament_count
summary_stats <- tubulin_UV %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p12.5 -uv tub", "bj p12.5 +uv tub")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(filament_count, na.rm = TRUE),
    sd = sd(filament_count, na.rm = TRUE),
    ymax = max(filament_count, na.rm = TRUE)
  )

ggplot(mutate(tubulin_UV, condition_clean = fct_relevel(condition_clean, "bj p12.5 -uv tub", "bj p12.5 +uv tub")), aes(condition_clean, filament_count, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Condition", y = "Number of filaments /cell") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('p12.5 -UV', 'p12.5 +UV')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(10000)) +
  geom_text(
    data = summary_stats,
    aes(
      x = condition_clean,
      y = 240000,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = tubulin_UV, comparisons = list(c("bj p12.5 -uv tub","bj p12.5 +uv tub")), annotations = "ns", map_signif_level = TRUE, y_position = 250000, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/tubUVFilCount.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")

# boxplot for tubulin_UV sum_branches
summary_stats <- tubulin_UV %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p12.5 -uv tub", "bj p12.5 +uv tub")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(sum_branches, na.rm = TRUE),
    sd = sd(sum_branches, na.rm = TRUE),
    ymax = max(sum_branches, na.rm = TRUE)
  )

ggplot(mutate(tubulin_UV, condition_clean = fct_relevel(condition_clean, "bj p12.5 -uv tub", "bj p12.5 +uv tub")), aes(condition_clean, sum_branches, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Condition", y = "Number of branches /cell") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('p12.5 -UV', 'p12.5 +UV')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(10000)) +
  geom_text(
    data = summary_stats,
    aes(
      x = condition_clean,
      y = 160000,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = tubulin_UV, comparisons = list(c("bj p12.5 -uv tub","bj p12.5 +uv tub")), annotations = "ns", map_signif_level = TRUE, y_position = 170000, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/tubUVSumBranch.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")

# boxplot for tubulin_UV mean_branches
summary_stats <- tubulin_UV %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p12.5 -uv tub", "bj p12.5 +uv tub")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_branches, na.rm = TRUE),
    sd = sd(mean_branches, na.rm = TRUE),
    ymax = max(mean_branches, na.rm = TRUE)
  )

ggplot(mutate(tubulin_UV, condition_clean = fct_relevel(condition_clean, "bj p12.5 -uv tub", "bj p12.5 +uv tub")), aes(condition_clean, mean_branches, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Condition", y = "Mean branches per filament /cell") +
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
      y = 0.74,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = tubulin_UV, comparisons = list(c("bj p12.5 -uv tub","bj p12.5 +uv tub")), annotations = "ns", map_signif_level = TRUE, y_position = 0.75, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/tubUVMeanBranch.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")

# boxplot for tubulin_primary sum_length_2d
summary_stats <- tubulin_primary %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "primary fibro 18yo tub", "primary fibro 29yo tub", "primary fibro 32yo tub", "primary fibro 59yo tub")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(sum_length_2d, na.rm = TRUE),
    sd = sd(sum_length_2d, na.rm = TRUE),
    ymax = max(sum_length_2d, na.rm = TRUE)
  )

ggplot(mutate(tubulin_primary, condition_clean = fct_relevel(condition_clean, "primary fibro 18yo tub", "primary fibro 29yo tub", "primary fibro 32yo tub", "primary fibro 59yo tub")), aes(condition_clean, sum_length_2d, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Age", y = "Total 2D filament length per cell (μm)") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('18 y/o', '29 y/o', '32 y/o', '59 y/o')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(10000)) +
  geom_text(
    data = summary_stats,
    aes(
      x = condition_clean,
      y = 340000,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = tubulin_primary, comparisons = list(c("primary fibro 18yo tub", "primary fibro 59yo tub")), annotations = "*", map_signif_level = TRUE, vjust = 0.7, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/tubPrimLength2d.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")

# boxplot for tubulin_primary mean_length_2d
summary_stats <- tubulin_primary %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "primary fibro 18yo tub", "primary fibro 29yo tub", "primary fibro 32yo tub", "primary fibro 59yo tub")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_length_2d, na.rm = TRUE),
    sd = sd(mean_length_2d, na.rm = TRUE),
    ymax = max(mean_length_2d, na.rm = TRUE)
  )

ggplot(mutate(tubulin_primary, condition_clean = fct_relevel(condition_clean, "primary fibro 18yo tub", "primary fibro 29yo tub", "primary fibro 32yo tub", "primary fibro 59yo tub")), aes(condition_clean, mean_length_2d, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Age", y = "Mean 2D filament length per cell (μm)") +
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
      y = 0.83,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = tubulin_primary, comparisons = list(c("primary fibro 18yo tub", "primary fibro 59yo tub")), annotations = "ns", map_signif_level = TRUE, y_position = 0.84, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/tubPrimMean2D.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")

# boxplot for tubulin_primary sum_length_3d
summary_stats <- tubulin_primary %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "primary fibro 18yo tub", "primary fibro 29yo tub", "primary fibro 32yo tub", "primary fibro 59yo tub")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(sum_length_3d, na.rm = TRUE),
    sd = sd(sum_length_3d, na.rm = TRUE),
    ymax = max(sum_length_3d, na.rm = TRUE)
  )

ggplot(mutate(tubulin_primary, condition_clean = fct_relevel(condition_clean, "primary fibro 18yo tub", "primary fibro 29yo tub", "primary fibro 32yo tub", "primary fibro 59yo tub")), aes(condition_clean, sum_length_3d, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Age", y = "Total 3D filament length per cell (μm)") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('18 y/o', '29 y/o', '32 y/o', '59 y/o')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(250000)) +
  geom_text(
    data = summary_stats,
    aes(
      x = condition_clean,
      y = 6800000,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = tubulin_primary, comparisons = list(c("primary fibro 18yo tub", "primary fibro 59yo tub")), annotations = "*", map_signif_level = TRUE, vjust = 0.7, y_position = , color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/tubPrimSum3D.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")

# boxplot for tubulin_primary mean_length_3d
summary_stats <- tubulin_primary %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "primary fibro 18yo tub", "primary fibro 29yo tub", "primary fibro 32yo tub", "primary fibro 59yo tub")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_length_3d, na.rm = TRUE),
    sd = sd(mean_length_3d, na.rm = TRUE),
    ymax = max(mean_length_3d, na.rm = TRUE)
  )

ggplot(mutate(tubulin_primary, condition_clean = fct_relevel(condition_clean, "primary fibro 18yo tub", "primary fibro 29yo tub", "primary fibro 32yo tub", "primary fibro 59yo tub")), aes(condition_clean, mean_length_3d, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Age", y = "Mean 3D filament length per cell (μm)") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('18 y/o', '29 y/o', '32 y/o', '59 y/o')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(0.25)) +
  geom_text(
    data = summary_stats,
    aes(
      x = condition_clean,
      y = 13.5,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = tubulin_primary, comparisons = list(c("primary fibro 18yo tub", "primary fibro 59yo tub")), annotations = "**", map_signif_level = TRUE, vjust = 0.7, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/tubPrimmean3d.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")
  
# boxplot for tubulin_primary mean_avg_angle
summary_stats <- tubulin_primary %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "primary fibro 18yo tub", "primary fibro 29yo tub", "primary fibro 32yo tub", "primary fibro 59yo tub")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_avg_angle, na.rm = TRUE),
    sd = sd(mean_avg_angle, na.rm = TRUE),
    ymax = max(mean_avg_angle, na.rm = TRUE)
  )

ggplot(mutate(tubulin_primary, condition_clean = fct_relevel(condition_clean, "primary fibro 18yo tub", "primary fibro 29yo tub", "primary fibro 32yo tub", "primary fibro 59yo tub")), aes(condition_clean, mean_avg_angle, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Age", y = "Mean average filament angle (deg)") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('18 y/o', '29 y/o', '32 y/o', '59 y/o')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(0.25)) +
  geom_text(
    data = summary_stats,
    aes(
      x = condition_clean,
      y = 85.8,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = tubulin_primary, comparisons = list(c("primary fibro 18yo tub", "primary fibro 59yo tub")), annotations = "*", map_signif_level = TRUE, vjust = 0.7, y_position = 86.2, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/tubPrimmeanAvgAng.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")

# boxplot for tubulin_primary mean_acc_curv
summary_stats <- tubulin_primary %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "primary fibro 18yo tub", "primary fibro 29yo tub", "primary fibro 32yo tub", "primary fibro 59yo tub")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_acc_curv, na.rm = TRUE),
    sd = sd(mean_acc_curv, na.rm = TRUE),
    ymax = max(mean_acc_curv, na.rm = TRUE)
  )

ggplot(mutate(tubulin_primary, condition_clean = fct_relevel(condition_clean, "primary fibro 18yo tub", "primary fibro 29yo tub", "primary fibro 32yo tub", "primary fibro 59yo tub")), aes(condition_clean, mean_acc_curv, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Age", y = "Mean accumulated filament curvature /cell") +
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
      y = 4.21,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = tubulin_primary, comparisons = list(c("primary fibro 18yo tub", "primary fibro 59yo tub")), annotations = "ns", map_signif_level = TRUE, y_position = 4.22, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/tubPrimmeanAccCurv.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")

# boxplot for tubulin_primary mean_net_curv
summary_stats <- tubulin_primary %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "primary fibro 18yo tub", "primary fibro 29yo tub", "primary fibro 32yo tub", "primary fibro 59yo tub")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_net_curv, na.rm = TRUE),
    sd = sd(mean_net_curv, na.rm = TRUE),
    ymax = max(mean_net_curv, na.rm = TRUE)
  )

ggplot(mutate(tubulin_primary, condition_clean = fct_relevel(condition_clean, "primary fibro 18yo tub", "primary fibro 29yo tub", "primary fibro 32yo tub", "primary fibro 59yo tub")), aes(condition_clean, mean_net_curv, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Age", y = "Mean end-to-end filament curvature /cell") +
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
      y = 1.5,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = tubulin_primary, comparisons = list(c("primary fibro 18yo tub", "primary fibro 59yo tub")), annotations = "*", map_signif_level = TRUE, vjust = 0.7, y_position = 1.52, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/tubPrimmeanNetCurv.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")


# boxplot for tubulin_primary filament_count
summary_stats <- tubulin_primary %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "primary fibro 18yo tub", "primary fibro 29yo tub", "primary fibro 32yo tub", "primary fibro 59yo tub")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(filament_count, na.rm = TRUE),
    sd = sd(filament_count, na.rm = TRUE),
    ymax = max(filament_count, na.rm = TRUE)
  )

ggplot(mutate(tubulin_primary, condition_clean = fct_relevel(condition_clean, "primary fibro 18yo tub", "primary fibro 29yo tub", "primary fibro 32yo tub", "primary fibro 59yo tub")), aes(condition_clean, filament_count, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Age", y = "Number of filaments /cell") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('18 y/o', '29 y/o', '32 y/o', '59 y/o')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(20000)) +
  geom_text(
    data = summary_stats,
    aes(
      x = condition_clean,
      y = 480000,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = tubulin_primary, comparisons = list(c("primary fibro 18yo tub", "primary fibro 59yo tub")), annotations = "*", map_signif_level = TRUE, vjust = 0.7, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/tubPrimFilCount.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")

# boxplot for tubulin_primary sum_branches
summary_stats <- tubulin_primary %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "primary fibro 18yo tub", "primary fibro 29yo tub", "primary fibro 32yo tub", "primary fibro 59yo tub")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(sum_branches, na.rm = TRUE),
    sd = sd(sum_branches, na.rm = TRUE),
    ymax = max(sum_branches, na.rm = TRUE)
  )

ggplot(mutate(tubulin_primary, condition_clean = fct_relevel(condition_clean, "primary fibro 18yo tub", "primary fibro 29yo tub", "primary fibro 32yo tub", "primary fibro 59yo tub")), aes(condition_clean, sum_branches, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Age", y = "Total filament branches /cell") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('18 y/o', '29 y/o', '32 y/o', '59 y/o')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(10000)) +
  geom_text(
    data = summary_stats,
    aes(
      x = condition_clean,
      y = 350000,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = tubulin_primary, comparisons = list(c("primary fibro 18yo tub", "primary fibro 59yo tub")), annotations = "*", map_signif_level = TRUE, vjust = 0.7, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/tubPrimBranchCount.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")

# boxplot for tubulin_primary mean_branches
summary_stats <- tubulin_primary %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "primary fibro 18yo tub", "primary fibro 29yo tub", "primary fibro 32yo tub", "primary fibro 59yo tub")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_branches, na.rm = TRUE),
    sd = sd(mean_branches, na.rm = TRUE),
    ymax = max(mean_branches, na.rm = TRUE)
  )

ggplot(mutate(tubulin_primary, condition_clean = fct_relevel(condition_clean, "primary fibro 18yo tub", "primary fibro 29yo tub", "primary fibro 32yo tub", "primary fibro 59yo tub")), aes(condition_clean, mean_branches, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Age", y = "Mean filament branches /cell") +
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
      y = 0.77,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = tubulin_primary, comparisons = list(c("primary fibro 18yo tub", "primary fibro 59yo tub")), annotations = "**", map_signif_level = TRUE, vjust = 0.7, y_position = 0.78, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/tubPrimMeanBranch.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")

# boxplot for tubulin_primary mean_deviation
summary_stats <- tubulin_primary %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "primary fibro 18yo tub", "primary fibro 29yo tub", "primary fibro 32yo tub", "primary fibro 59yo tub")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_deviation, na.rm = TRUE),
    sd = sd(mean_deviation, na.rm = TRUE),
    ymax = max(mean_deviation, na.rm = TRUE)
  )

ggplot(mutate(tubulin_primary, condition_clean = fct_relevel(condition_clean, "primary fibro 18yo tub", "primary fibro 29yo tub", "primary fibro 32yo tub", "primary fibro 59yo tub")), aes(condition_clean, mean_deviation, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Age", y = "Mean filament deviatiom (μm)") +
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
      y = 2.53,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = tubulin_primary, comparisons = list(c("primary fibro 18yo tub", "primary fibro 59yo tub")), annotations = "ns", map_signif_level = TRUE, y_position = 2.6, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/tubPrimMeanDev.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")

# doplot for tubulin_primary sum_length_3d vs sum_branches
ggplot(data = tubulin_primary, aes(x = sum_length_3d, y = sum_branches, fill = condition_clean, color = condition_clean)) + 
  geom_point(size = 0.8) + 
  labs(x = "Total 3D filament length (μm)", y = "Total branches /cell") +
  theme_bw() +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, y = guide_axis(minor.ticks = TRUE), x = guide_axis(minor.ticks = TRUE)) +
  scale_y_continuous(minor_breaks = scales::breaks_width(5000)) 
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/tubPrimSum2D.tiff", device = "tiff", dpi = 300, width = 177, height = 150, units = "mm")