library(tidyverse)
library(lme4)
library(car)
library(emmeans)
library(ggprism)
library(forcats)
library(dplyr)
library(ggsignif)
library(DHARMa)

vimentin_df <- read.csv("vimentin cell level.csv")
vimentin_BJ <- vimentin_df[0:253,]
vimentin_primary <- vimentin_df[254:314,]
vimentin_UV <- filter(vimentin_df, condition_clean %in% c("bj p12.5 -uv vim", "bj p12.5 +uv vim"))
vimentin_etop <- filter(vimentin_df, condition_clean %in% c("bj p8.5 vim", "bj p9.5 etop vim"))
vimentin_pass <- filter(vimentin_df, condition_clean %in% c("bj p11.5 vim", "bj p26.5 vim", "bj p31.5 vim", "bj p36.5 vim", "bj p38.5 vim", "bj p40.5 vim"))

#stats for all vimentin BJ datasets
model <- lmer(sum_length_2d ~ condition_clean + (1|condition_clean:biorep) +(1|condition_clean:biorep:techrep), vimentin_BJ)
summary(model)
Anova(model)
em1 <- emmeans(model, ~condition_clean)
pairs(em1)

#stats for UV only
model_UV <- lmer(sum_length_2d ~ condition_clean + (1|condition_clean:biorep) +(1|condition_clean:biorep:techrep), vimentin_UV)
help('isSingular')
summary(model_UV)
Anova(model_UV)
em_UV <- emmeans(model_UV, ~condition_clean)
pairs(em_UV)

#stats for etop only
model_etop <- lmer(mean_branches ~ condition_clean + (1|condition_clean:biorep:techrep), vimentin_etop)
summary(model_etop)
Anova(model_etop)
em_etop <- emmeans(model_etop, ~condition_clean)
pairs(em_etop)

simulationOutput <- simulateResiduals(fittedModel = model_etop, plot = T)

kruskal <- kruskal.test(sum_branches ~ condition_clean, data = vimentin_etop)
kruskal

#stats for passage only sum_length_2d
model_pass <- lmer(sum_length_2d ~ condition_clean + (1|condition_clean:biorep) + (1|condition_clean:biorep:techrep), vimentin_pass)
summary(model_pass)
Anova(model_pass)
em_pass <- emmeans(model_pass, ~condition_clean)
pairs(em_pass)

simulationOutput <- simulateResiduals(fittedModel = model_pass, plot = T)

kruskal <- kruskal.test(sum_length_2d ~ condition_clean, data = vimentin_pass)
kruskal
pairwise.wilcox.test(vimentin_pass$sum_length_2d, vimentin_pass$condition_clean, p.adjust.method = "BH")

#stats for passage only mean_length_2d
model2_pass <- lmer(mean_length_2d ~ condition_clean + (1|condition_clean:biorep:techrep), vimentin_pass)
summary(model2_pass)
Anova(model2_pass)
em_pass2 <- emmeans(model2_pass, ~condition_clean)
pairs(em_pass2)

simulationOutput <- simulateResiduals(fittedModel = model2_pass, plot = T)

#stats for passage only sum_length_3d
model3_pass <- lmer(sum_length_3d ~ condition_clean + (1|condition_clean:biorep) + (1|condition_clean:biorep:techrep), vimentin_pass)
summary(model3_pass)
Anova(model3_pass)
em_pass3 <- emmeans(model3_pass, ~condition_clean)
pairs(em_pass3)

simulationOutput <- simulateResiduals(fittedModel = model3_pass, plot = T)

kruskal <- kruskal.test(sum_length_3d ~ condition_clean, data = vimentin_pass)
kruskal
pairwise.wilcox.test(vimentin_pass$sum_length_3d, vimentin_pass$condition_clean, p.adjust.method = "BH")

#stats for passage only mean_length_3d
model4_pass <- lmer(mean_length_3d ~ condition_clean + (1|condition_clean:biorep:techrep), vimentin_pass)
summary(model4_pass)
Anova(model4_pass)
em_pass4 <- emmeans(model4_pass, ~condition_clean)
pairs(em_pass4)

simulationOutput <- simulateResiduals(fittedModel = model4_pass, plot = T)

#stats for passage only mean_avg_angle
model5_pass <- lmer(mean_avg_angle ~ condition_clean + (1|condition_clean:biorep) + (1|condition_clean:biorep:techrep), vimentin_pass)
summary(model5_pass)
Anova(model5_pass)
em_pass5 <- emmeans(model5_pass, ~condition_clean)
pairs(em_pass5)

simulationOutput <- simulateResiduals(fittedModel = model5_pass, plot = T)

#stats for passage only mean_acc_curv
model6_pass <- lmer(mean_net_curv ~ condition_clean + (1|condition_clean:biorep:techrep), vimentin_pass)
summary(model6_pass)
Anova(model6_pass)
em_pass6 <- emmeans(model6_pass, ~condition_clean)
pairs(em_pass6)

simulationOutput <- simulateResiduals(fittedModel = model6_pass, plot = T)

#stats for passage only filament_count
model7_pass <- lmer(filament_count ~ condition_clean + (1|condition_clean:biorep) + (1|condition_clean:biorep:techrep), vimentin_pass)
summary(model7_pass)
Anova(model7_pass)
em_pass7 <- emmeans(model7_pass, ~condition_clean)
pairs(em_pass7)

simulationOutput <- simulateResiduals(fittedModel = model7_pass, plot = T)

kruskal <- kruskal.test(filament_count ~ condition_clean, data = vimentin_pass)
kruskal
pairwise.wilcox.test(vimentin_pass$filament_count, vimentin_pass$condition_clean, p.adjust.method = "BH")

#stats for passage only sum_branches
model8_pass <- lmer(sum_branches ~ condition_clean + (1|condition_clean:biorep) + (1|condition_clean:biorep:techrep), vimentin_pass)
summary(model8_pass)
Anova(model8_pass)
em_pass8 <- emmeans(model8_pass, ~condition_clean)
pairs(em_pass8)

simulationOutput <- simulateResiduals(fittedModel = model8_pass, plot = T)

kruskal <- kruskal.test(sum_branches ~ condition_clean, data = vimentin_pass)
kruskal
pairwise.wilcox.test(vimentin_pass$sum_branches, vimentin_pass$condition_clean, p.adjust.method = "BH")

#stats for passage only mean_branches
model9_pass <- lmer(mean_branches ~ condition_clean + (1|condition_clean:biorep) + (1|condition_clean:biorep:techrep), vimentin_pass)
summary(model9_pass)
Anova(model9_pass)
em_pass9 <- emmeans(model9_pass, ~condition_clean)
pairs(em_pass9)

simulationOutput <- simulateResiduals(fittedModel = model9_pass, plot = T)

#stats for passage only mean_deviation
model10_pass <- lmer(mean_deviation ~ condition_clean + (1|condition_clean:biorep:techrep), vimentin_pass)
summary(model10_pass)
Anova(model10_pass)
em_pass10 <- emmeans(model10_pass, ~condition_clean)
pairs(em_pass10)

simulationOutput <- simulateResiduals(fittedModel = model10_pass, plot = T)

#stats for UV 
model_UV <- lmer(mean_branches ~ condition_clean + (1|condition_clean:biorep:techrep), vimentin_UV)
summary(model_UV)
Anova(model_UV)
em_UV <- emmeans(model_UV, ~condition_clean)
pairs(em_UV)


simulationOutput <- simulateResiduals(fittedModel = model_UV, plot = T)

#stats for vimentin primary
model_prim <- lmer(mean_branches ~ condition_clean + (1|condition_clean:biorep:techrep), vimentin_primary)
summary(model_prim)
Anova(model_prim)
em13 <- emmeans(model_prim, ~condition_clean)
pairs(em13)

simulationOutput <- simulateResiduals(fittedModel = model_prim, plot = T)

kruskal <- kruskal.test(sum_branches ~ condition_clean, data = vimentin_primary)
kruskal
pairwise.wilcox.test(vimentin_primary$sum_branches, vimentin_primary$condition_clean, p.adjust.method = "BH")

# boxplot for all vimentin sum_length_2D
ggplot(mutate(vimentin_BJ, condition_clean = fct_relevel(condition_clean, "bj p8.5 vim", "bj p9.5 etop vim", "bj p11.5 vim", "bj p12.5 -uv vim", "bj p12.5 +uv vim", "bj p26.5 vim", "bj p31.5 vim", "bj p36.5 vim", "bj p38.5 vim", "bj p40.5 vim")), aes(condition_clean, sum_length_2d, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(size = 0.5, position = position_jitterdodge()) +
  labs(x = "Condition", y = "Total 2D filament length per cell (μm)") +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('p8.5', 'p9.5 +etop', 'p11.5', 'p12.5 -UV', 'p12.5 +UV', 'p26.5', 'p31.5', 'p36.5', 'p38.5', 'p40.5')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(5000))
ggsave("vimBJsum2D.tiff", device = "tiff", dpi = 300, width = 177, height = 150, units = "mm")

# boxplot for all vimentin mean_length_2D
ggplot(mutate(vimentin_BJ, condition_clean = fct_relevel(condition_clean, "bj p8.5 vim", "bj p9.5 etop vim", "bj p11.5 vim", "bj p12.5 -uv vim", "bj p12.5 +uv vim", "bj p26.5 vim", "bj p31.5 vim", "bj p36.5 vim", "bj p38.5 vim", "bj p40.5 vim")), aes(condition_clean, mean_length_2d, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(size = 0.5, position = position_jitterdodge()) +
  labs(x = "Condition", y = "Mean 2D filament length per cell (μm)") +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('p8.5', 'p9.5 +etop', 'p11.5', 'p12.5 -UV', 'p12.5 +UV', 'p26.5', 'p31.5', 'p36.5', 'p38.5', 'p40.5')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(0.01))
ggsave("vimBJmean2D.tiff", device = "tiff", dpi = 300, width = 177, height = 150, units = "mm")

# boxplot for all vimentin sum_length_3D
ggplot(mutate(vimentin_BJ, condition_clean = fct_relevel(condition_clean, "bj p8.5 vim", "bj p9.5 etop vim", "bj p11.5 vim", "bj p12.5 -uv vim", "bj p12.5 +uv vim", "bj p26.5 vim", "bj p31.5 vim", "bj p36.5 vim", "bj p38.5 vim", "bj p40.5 vim")), aes(condition_clean, sum_length_3d, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(size = 0.5, position = position_jitterdodge()) +
  labs(x = "Condition", y = "Total 3D filament length per cell (μm)") +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('p8.5', 'p9.5 +etop', 'p11.5', 'p12.5 -UV', 'p12.5 +UV', 'p26.5', 'p31.5', 'p36.5', 'p38.5', 'p40.5')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(50000))
ggsave("vimBJsum3D.tiff", device = "tiff", dpi = 300, width = 177, height = 150, units = "mm")

# boxplot for all vimentin mean_length_3D
ggplot(mutate(vimentin_BJ, condition_clean = fct_relevel(condition_clean, "bj p8.5 vim", "bj p9.5 etop vim", "bj p11.5 vim", "bj p12.5 -uv vim", "bj p12.5 +uv vim", "bj p26.5 vim", "bj p31.5 vim", "bj p36.5 vim", "bj p38.5 vim", "bj p40.5 vim")), aes(condition_clean, mean_length_3d, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(size = 0.5, position = position_jitterdodge()) +
  labs(x = "Condition", y = "Mean 3D filament length per cell (μm)") +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('p8.5', 'p9.5 +etop', 'p11.5', 'p12.5 -UV', 'p12.5 +UV', 'p26.5', 'p31.5', 'p36.5', 'p38.5', 'p40.5')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(0.5))
ggsave("vimBJmean3D.tiff", device = "tiff", dpi = 300, width = 177, height = 150, units = "mm")


# boxplot for vimentin_pass sum_length_2d
summary_stats <- vimentin_pass %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p11.5 vim","bj p26.5 vim","bj p31.5 vim",
                                       "bj p36.5 vim","bj p38.5 vim","bj p40.5 vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(sum_length_2d, na.rm = TRUE),
    sd = sd(sum_length_2d, na.rm = TRUE),
    ymax = max(sum_length_2d, na.rm = TRUE)
  )

ggplot(mutate(vimentin_pass, condition_clean = fct_relevel(condition_clean, "bj p11.5 vim", "bj p26.5 vim", "bj p31.5 vim", "bj p36.5 vim", "bj p38.5 vim", "bj p40.5 vim")), aes(condition_clean, sum_length_2d, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Passage", y = "Total 2D filament length per cell (μm)") +
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
      y = 52000,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 2.5
  ) +
  geom_signif(data = vimentin_pass, comparisons = list(c("bj p11.5 vim","bj p40.5 vim")), map_signif_level = TRUE, annotations = "***", vjust = 0.7, y_position = 62000, color = 'black') +
  geom_signif(data = vimentin_pass, comparisons = list(c("bj p11.5 vim","bj p38.5 vim")), map_signif_level = TRUE, annotations = "***", vjust = 0.7, y_position = 60000, color = 'black') +
  geom_signif(data = vimentin_pass, comparisons = list(c("bj p11.5 vim","bj p36.5 vim")), map_signif_level = TRUE, annotations = "***", vjust = 0.7, y_position = 58000, color = 'black') +
  geom_signif(data = vimentin_pass, comparisons = list(c("bj p11.5 vim","bj p31.5 vim")), map_signif_level = TRUE, annotations = "**", vjust = 0.7, y_position = 56000, color = 'black') +
  geom_signif(data = vimentin_pass, comparisons = list(c("bj p11.5 vim","bj p26.5 vim")), map_signif_level = TRUE, annotations = "***", vjust = 0.7, y_position = 54000, color = 'black')

ggsave("vimPassSum2D.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")


# boxplot for vimentin_pass mean_length_2d
summary_stats <- vimentin_pass %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p11.5 vim","bj p26.5 vim","bj p31.5 vim",
                                       "bj p36.5 vim","bj p38.5 vim","bj p40.5 vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_length_2d, na.rm = TRUE),
    sd = sd(mean_length_2d, na.rm = TRUE),
    ymax = max(mean_length_2d, na.rm = TRUE)
  )

ggplot(mutate(vimentin_pass, condition_clean = fct_relevel(condition_clean, "bj p11.5 vim", "bj p26.5 vim", "bj p31.5 vim", "bj p36.5 vim", "bj p38.5 vim", "bj p40.5 vim")), aes(condition_clean, mean_length_2d, fill = condition_clean, color = condition_clean)) + 
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
      y = 0.985,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = vimentin_pass, comparisons = list(c("bj p11.5 vim","bj p36.5 vim")), map_signif_level = TRUE, annotations = "**", vjust = 0.7, color = 'black') +
  geom_signif(data = vimentin_pass, comparisons = list(c("bj p11.5 vim","bj p40.5 vim")), map_signif_level = TRUE, annotations = "*", y_position = 1.025, vjust = 0.7, color = 'black') 
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/vimPassMean2D.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")


# boxplot for vimentin_pass sum_length_3d
summary_stats <- vimentin_pass %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p11.5 vim","bj p26.5 vim","bj p31.5 vim",
                                       "bj p36.5 vim","bj p38.5 vim","bj p40.5 vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(sum_length_3d, na.rm = TRUE),
    sd = sd(sum_length_3d, na.rm = TRUE),
    ymax = max(sum_length_3d, na.rm = TRUE)
  )

ggplot(mutate(vimentin_pass, condition_clean = fct_relevel(condition_clean, "bj p11.5 vim", "bj p26.5 vim", "bj p31.5 vim", "bj p36.5 vim", "bj p38.5 vim", "bj p40.5 vim")), aes(condition_clean, sum_length_3d, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Passage", y = "Total 3D filament length per cell (μm)") +
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
      y = 430000,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 2.3
  ) +
  geom_signif(data = vimentin_pass, comparisons = list(c("bj p11.5 vim","bj p40.5 vim")), map_signif_level = TRUE, annotations = "***", vjust = 0.7, y_position = 560000, color = 'black') +
  geom_signif(data = vimentin_pass, comparisons = list(c("bj p11.5 vim","bj p38.5 vim")), map_signif_level = TRUE, annotations = "***", vjust = 0.7, y_position = 540000, color = 'black') +
  geom_signif(data = vimentin_pass, comparisons = list(c("bj p11.5 vim","bj p36.5 vim")), map_signif_level = TRUE, annotations = "***", vjust = 0.7, y_position = 520000, color = 'black') +
  geom_signif(data = vimentin_pass, comparisons = list(c("bj p11.5 vim","bj p31.5 vim")), map_signif_level = TRUE, annotations = "*", vjust = 0.7, y_position = 500000, color = 'black') +
  geom_signif(data = vimentin_pass, comparisons = list(c("bj p11.5 vim","bj p26.5 vim")), map_signif_level = TRUE, annotations = "***", vjust = 0.7, y_position = 480000, color = 'black')
ggsave("vimPassSum3D.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")


# boxplot for vimentin_pass mean_length_3d
summary_stats <- vimentin_pass %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p11.5 vim","bj p26.5 vim","bj p31.5 vim",
                                       "bj p36.5 vim","bj p38.5 vim","bj p40.5 vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_length_3d, na.rm = TRUE),
    sd = sd(mean_length_3d, na.rm = TRUE),
    ymax = max(mean_length_3d, na.rm = TRUE)
  )

ggplot(mutate(vimentin_pass, condition_clean = fct_relevel(condition_clean, "bj p11.5 vim", "bj p26.5 vim", "bj p31.5 vim", "bj p36.5 vim", "bj p38.5 vim", "bj p40.5 vim")), aes(condition_clean, mean_length_3d, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Passage", y = "Mean 3D filament length per cell (μm)") +
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
      y = 9.5,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = vimentin_pass, comparisons = list(c("bj p11.5 vim","bj p40.5 vim")), map_signif_level = TRUE, annotations = "ns", y_position = 9.8, color = 'black') 
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/vimPassMean3D.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")


# boxplot for vimentin_pass mean_avg_angle
summary_stats <- vimentin_pass %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p11.5 vim","bj p26.5 vim","bj p31.5 vim",
                                       "bj p36.5 vim","bj p38.5 vim","bj p40.5 vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_avg_angle, na.rm = TRUE),
    sd = sd(mean_avg_angle, na.rm = TRUE),
    ymax = max(mean_avg_angle, na.rm = TRUE)
  )

ggplot(mutate(vimentin_pass, condition_clean = fct_relevel(condition_clean, "bj p11.5 vim", "bj p26.5 vim", "bj p31.5 vim", "bj p36.5 vim", "bj p38.5 vim", "bj p40.5 vim")), aes(condition_clean, mean_avg_angle, fill = condition_clean, color = condition_clean)) + 
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
      y = 86,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = vimentin_pass, comparisons = list(c("bj p11.5 vim","bj p40.5 vim")), map_signif_level = TRUE, annotations = "ns", y_position = 87, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/vimPassMeanAvgAng.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")


# boxplot for vimentin_pass mean_acc_curv
summary_stats <- vimentin_pass %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p11.5 vim","bj p26.5 vim","bj p31.5 vim",
                                       "bj p36.5 vim","bj p38.5 vim","bj p40.5 vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_acc_curv, na.rm = TRUE),
    sd = sd(mean_acc_curv, na.rm = TRUE),
    ymax = max(mean_acc_curv, na.rm = TRUE)
  )

ggplot(mutate(vimentin_pass, condition_clean = fct_relevel(condition_clean, "bj p11.5 vim", "bj p26.5 vim", "bj p31.5 vim", "bj p36.5 vim", "bj p38.5 vim", "bj p40.5 vim")), aes(condition_clean, mean_acc_curv, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Passage", y = "Mean accumulated filament curvature /cell") +
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
      y = 4.225,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = vimentin_pass, comparisons = list(c("bj p11.5 vim","bj p40.5 vim")), map_signif_level = TRUE, annotations = "ns", y_position = 4.24, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/vimPassMeanAccCurv.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")


# boxplot for vimentin_pass mean_net_curv
summary_stats <- vimentin_pass %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p11.5 vim","bj p26.5 vim","bj p31.5 vim",
                                       "bj p36.5 vim","bj p38.5 vim","bj p40.5 vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_net_curv, na.rm = TRUE),
    sd = sd(mean_net_curv, na.rm = TRUE),
    ymax = max(mean_net_curv, na.rm = TRUE)
  )

ggplot(mutate(vimentin_pass, condition_clean = fct_relevel(condition_clean, "bj p11.5 vim", "bj p26.5 vim", "bj p31.5 vim", "bj p36.5 vim", "bj p38.5 vim", "bj p40.5 vim")), aes(condition_clean, mean_net_curv, fill = condition_clean, color = condition_clean)) + 
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
      y = 1.54,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = vimentin_pass, comparisons = list(c("bj p11.5 vim","bj p40.5 vim")), map_signif_level = TRUE, annotations = "ns", y_position = 1.55, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/vimPassMeanNetCurv.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")


# boxplot for vimentin_pass filament_count
summary_stats <- vimentin_pass %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p11.5 vim","bj p26.5 vim","bj p31.5 vim",
                                       "bj p36.5 vim","bj p38.5 vim","bj p40.5 vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(filament_count, na.rm = TRUE),
    sd = sd(filament_count, na.rm = TRUE),
    ymax = max(filament_count, na.rm = TRUE)
  )

ggplot(mutate(vimentin_pass, condition_clean = fct_relevel(condition_clean, "bj p11.5 vim", "bj p26.5 vim", "bj p31.5 vim", "bj p36.5 vim", "bj p38.5 vim", "bj p40.5 vim")), aes(condition_clean, filament_count, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Passage", y = "Filament count /cell") +
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
      y = 68000,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 2.5
  ) +
  geom_signif(data = vimentin_pass, comparisons = list(c("bj p11.5 vim","bj p40.5 vim")), map_signif_level = TRUE, annotations = "***", vjust = 0.7, y_position = 80000, color = 'black') +
  geom_signif(data = vimentin_pass, comparisons = list(c("bj p11.5 vim","bj p38.5 vim")), map_signif_level = TRUE, annotations = "***", vjust = 0.7, y_position = 77500, color = 'black') +
  geom_signif(data = vimentin_pass, comparisons = list(c("bj p11.5 vim","bj p36.5 vim")), map_signif_level = TRUE, annotations = "***", vjust = 0.7, y_position = 75000, color = 'black') +
  geom_signif(data = vimentin_pass, comparisons = list(c("bj p11.5 vim","bj p31.5 vim")), map_signif_level = TRUE, annotations = "*", vjust = 0.7, y_position = 72500, color = 'black') +
  geom_signif(data = vimentin_pass, comparisons = list(c("bj p11.5 vim","bj p26.5 vim")), map_signif_level = TRUE, annotations = "***", vjust = 0.7, y_position = 70000, color = 'black')
ggsave("vimPassFilCount.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")


# boxplot for vimentin_pass sum_branches
summary_stats <- vimentin_pass %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p11.5 vim","bj p26.5 vim","bj p31.5 vim",
                                       "bj p36.5 vim","bj p38.5 vim","bj p40.5 vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(sum_branches, na.rm = TRUE),
    sd = sd(sum_branches, na.rm = TRUE),
    ymax = max(sum_branches, na.rm = TRUE)
  )

ggplot(mutate(vimentin_pass, condition_clean = fct_relevel(condition_clean, "bj p11.5 vim", "bj p26.5 vim", "bj p31.5 vim", "bj p36.5 vim", "bj p38.5 vim", "bj p40.5 vim")), aes(condition_clean, sum_branches, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Passage", y = "Total filament branches /cell") +
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
      y = 55000,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 2.5
  ) +
  geom_signif(data = vimentin_pass, comparisons = list(c("bj p11.5 vim","bj p40.5 vim")), map_signif_level = TRUE, annotations = "***", vjust = 0.7, y_position = 65000, color = 'black') +
  geom_signif(data = vimentin_pass, comparisons = list(c("bj p11.5 vim","bj p38.5 vim")), map_signif_level = TRUE, annotations = "***", vjust = 0.7, y_position = 63000, color = 'black') +
  geom_signif(data = vimentin_pass, comparisons = list(c("bj p11.5 vim","bj p36.5 vim")), map_signif_level = TRUE, annotations = "***", vjust = 0.7, y_position = 61000, color = 'black') +
  geom_signif(data = vimentin_pass, comparisons = list(c("bj p11.5 vim","bj p31.5 vim")), map_signif_level = TRUE, annotations = "**", vjust = 0.7, y_position = 59000, color = 'black') +
  geom_signif(data = vimentin_pass, comparisons = list(c("bj p11.5 vim","bj p26.5 vim")), map_signif_level = TRUE, annotations = "***", vjust = 0.7, y_position = 57000, color = 'black')
ggsave("vimPassBranchCount.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")


# boxplot for vimentin_pass mean_branches
summary_stats <- vimentin_pass %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p11.5 vim","bj p26.5 vim","bj p31.5 vim",
                                       "bj p36.5 vim","bj p38.5 vim","bj p40.5 vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_branches, na.rm = TRUE),
    sd = sd(mean_branches, na.rm = TRUE),
    ymax = max(mean_branches, na.rm = TRUE)
  )

ggplot(mutate(vimentin_pass, condition_clean = fct_relevel(condition_clean, "bj p11.5 vim", "bj p26.5 vim", "bj p31.5 vim", "bj p36.5 vim", "bj p38.5 vim", "bj p40.5 vim")), aes(condition_clean, mean_branches, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Passage", y = "Mean filament branches /cell") +
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
      y = 0.91,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = vimentin_pass, comparisons = list(c("bj p11.5 vim","bj p26.5 vim")), annotations = "**", map_signif_level = TRUE, vjust = 0.7, y_position = 0.92, color = 'black') +
  geom_signif(data = vimentin_pass, comparisons = list(c("bj p11.5 vim","bj p36.5 vim")), annotations = "*", map_signif_level = TRUE, vjust = 0.7, y_position = 0.93, color = 'black') +
  geom_signif(data = vimentin_pass, comparisons = list(c("bj p11.5 vim","bj p38.5 vim")), annotations = "*", map_signif_level = TRUE, vjust = 0.7, y_position = 0.94, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/vimPassMeanBranch.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")


# boxplot for vimentin_pass mean_deviation
summary_stats <- vimentin_pass %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p11.5 vim","bj p26.5 vim","bj p31.5 vim",
                                       "bj p36.5 vim","bj p38.5 vim","bj p40.5 vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_deviation, na.rm = TRUE),
    sd = sd(mean_deviation, na.rm = TRUE),
    ymax = max(mean_deviation, na.rm = TRUE)
  )

ggplot(mutate(vimentin_pass, condition_clean = fct_relevel(condition_clean, "bj p11.5 vim", "bj p26.5 vim", "bj p31.5 vim", "bj p36.5 vim", "bj p38.5 vim", "bj p40.5 vim")), aes(condition_clean, mean_deviation, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Passage", y = "Mean filament deviation (μm)") +
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
      y = 1.55,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = vimentin_pass, comparisons = list(c("bj p11.5 vim","bj p40.5 vim")), annotations = "ns", map_signif_level = TRUE, y_position = 1.6, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/vimPassMeanDev.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")

# doplot for vimentin_pass sum_length_3d vs sum_branches
ggplot(data = vimentin_pass, aes(x = sum_length_3d, y = sum_branches, fill = condition_clean, color = condition_clean)) + 
  geom_point(size = 0.8) + 
  labs(x = "Total 3D filament length (μm)", y = "Total branches /cell") +
  theme_bw() +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE), x = guide_axis(minor.ticks = TRUE)) +
  scale_y_continuous(minor_breaks = scales::breaks_width(5000)) 
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/vimPassSum2D.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")


# boxplot for vimentin_UV sum_length_2d
summary_stats <- vimentin_UV %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p12.5 -uv vim", "bj p12.5 +uv vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(sum_length_2d, na.rm = TRUE),
    sd = sd(sum_length_2d, na.rm = TRUE),
    ymax = max(sum_length_2d, na.rm = TRUE)
  )

ggplot(mutate(vimentin_UV, condition_clean = fct_relevel(condition_clean, "bj p12.5 -uv vim", "bj p12.5 +uv vim")), aes(condition_clean, sum_length_2d, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Condition", y = "Total 2D filament length per cell (μm)") +
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
      y = 27500,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = vimentin_UV, comparisons = list(c("bj p12.5 -uv vim","bj p12.5 +uv vim")), annotations = "ns", map_signif_level = TRUE, y_position = 29000, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/vimUVsum2d.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")


# boxplot for vimentin_UV mean_length_2d
summary_stats <- vimentin_UV %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p12.5 -uv vim", "bj p12.5 +uv vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_length_2d, na.rm = TRUE),
    sd = sd(mean_length_2d, na.rm = TRUE),
    ymax = max(mean_length_2d, na.rm = TRUE)
  )

ggplot(mutate(vimentin_UV, condition_clean = fct_relevel(condition_clean, "bj p12.5 -uv vim", "bj p12.5 +uv vim")), aes(condition_clean, mean_length_2d, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Condition", y = "Mean 2D filament length per cell (μm)") +
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
      y = 0.9,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = vimentin_UV, comparisons = list(c("bj p12.5 -uv vim","bj p12.5 +uv vim")), annotations = "ns", map_signif_level = TRUE, y_position = 0.91, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/vimUVmean2d.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")


# boxplot for vimentin_UV sum_length_3d
summary_stats <- vimentin_UV %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p12.5 -uv vim", "bj p12.5 +uv vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(sum_length_3d, na.rm = TRUE),
    sd = sd(sum_length_3d, na.rm = TRUE),
    ymax = max(sum_length_3d, na.rm = TRUE)
  )

ggplot(mutate(vimentin_UV, condition_clean = fct_relevel(condition_clean, "bj p12.5 -uv vim", "bj p12.5 +uv vim")), aes(condition_clean, sum_length_3d, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Condition", y = "Total 3D filament length per cell (μm)") +
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
      y = 225000,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = vimentin_UV, comparisons = list(c("bj p12.5 -uv vim","bj p12.5 +uv vim")), annotations = "ns", map_signif_level = TRUE, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/vimUVsum3d.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")


# boxplot for vimentin_UV mean_length_3d
summary_stats <- vimentin_UV %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p12.5 -uv vim", "bj p12.5 +uv vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_length_3d, na.rm = TRUE),
    sd = sd(mean_length_3d, na.rm = TRUE),
    ymax = max(mean_length_3d, na.rm = TRUE)
  )

ggplot(mutate(vimentin_UV, condition_clean = fct_relevel(condition_clean, "bj p12.5 -uv vim", "bj p12.5 +uv vim")), aes(condition_clean, mean_length_3d, fill = condition_clean, color = condition_clean)) + 
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
      y = 8.8,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = vimentin_UV, comparisons = list(c("bj p12.5 -uv vim","bj p12.5 +uv vim")), annotations = "ns", map_signif_level = TRUE, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/vimUVmean3d.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")


# boxplot for vimentin_UV mean_avg_angle
summary_stats <- vimentin_UV %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p12.5 -uv vim", "bj p12.5 +uv vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_avg_angle, na.rm = TRUE),
    sd = sd(mean_avg_angle, na.rm = TRUE),
    ymax = max(mean_avg_angle, na.rm = TRUE)
  )

ggplot(mutate(vimentin_UV, condition_clean = fct_relevel(condition_clean, "bj p12.5 -uv vim", "bj p12.5 +uv vim")), aes(condition_clean, mean_avg_angle, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Passage", y = "Mean average filament angle (deg)") +
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
      y = 85.7,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = vimentin_UV, comparisons = list(c("bj p12.5 -uv vim", "bj p12.5 +uv vim")), map_signif_level = TRUE, annotations = "ns", color = 'black') 
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/vimUVMeanAvgAng.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")


# boxplot for vimentin_UV mean_deviation
summary_stats <- vimentin_UV %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p12.5 -uv vim", "bj p12.5 +uv vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_deviation, na.rm = TRUE),
    sd = sd(mean_deviation, na.rm = TRUE),
    ymax = max(mean_deviation, na.rm = TRUE)
  )

ggplot(mutate(vimentin_UV, condition_clean = fct_relevel(condition_clean, "bj p12.5 -uv vim", "bj p12.5 +uv vim")), aes(condition_clean, mean_deviation, fill = condition_clean, color = condition_clean)) + 
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
      y = 1.43,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = vimentin_UV, comparisons = list(c("bj p12.5 -uv vim","bj p12.5 +uv vim")), annotations = "ns", map_signif_level = TRUE, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/vimUVmeanDev.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")


# boxplot for vimentin_UV mean_acc_curv
summary_stats <- vimentin_UV %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p12.5 -uv vim", "bj p12.5 +uv vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_acc_curv, na.rm = TRUE),
    sd = sd(mean_acc_curv, na.rm = TRUE),
    ymax = max(mean_acc_curv, na.rm = TRUE)
  )

ggplot(mutate(vimentin_UV, condition_clean = fct_relevel(condition_clean, "bj p12.5 -uv vim", "bj p12.5 +uv vim")), aes(condition_clean, mean_acc_curv, fill = condition_clean, color = condition_clean)) + 
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
  geom_signif(data = vimentin_UV, comparisons = list(c("bj p12.5 -uv vim","bj p12.5 +uv vim")), annotations = "ns", map_signif_level = TRUE, y_position = 4.21, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/vimUVmeanAccCurv.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")


# boxplot for vimentin_UV mean_net_curv
summary_stats <- vimentin_UV %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p12.5 -uv vim", "bj p12.5 +uv vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_net_curv, na.rm = TRUE),
    sd = sd(mean_net_curv, na.rm = TRUE),
    ymax = max(mean_net_curv, na.rm = TRUE)
  )

ggplot(mutate(vimentin_UV, condition_clean = fct_relevel(condition_clean, "bj p12.5 -uv vim", "bj p12.5 +uv vim")), aes(condition_clean, mean_net_curv, fill = condition_clean, color = condition_clean)) + 
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
      y = 1.44,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = vimentin_UV, comparisons = list(c("bj p12.5 -uv vim","bj p12.5 +uv vim")), annotations = "ns", map_signif_level = TRUE, y_position = 1.45, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/vimUVmeanNetCurv.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")


# boxplot for vimentin_UV filament_count
summary_stats <- vimentin_UV %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p12.5 -uv vim", "bj p12.5 +uv vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(filament_count, na.rm = TRUE),
    sd = sd(filament_count, na.rm = TRUE),
    ymax = max(filament_count, na.rm = TRUE)
  )

ggplot(mutate(vimentin_UV, condition_clean = fct_relevel(condition_clean, "bj p12.5 -uv vim", "bj p12.5 +uv vim")), aes(condition_clean, filament_count, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Condition", y = "Number of filaments /cell") +
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
      y = 30300,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = vimentin_UV, comparisons = list(c("bj p12.5 -uv vim","bj p12.5 +uv vim")), annotations = "ns", map_signif_level = TRUE, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/vimUVFilCount.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")


# boxplot for vimentin_UV sum_branches
summary_stats <- vimentin_UV %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p12.5 -uv vim", "bj p12.5 +uv vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(sum_branches, na.rm = TRUE),
    sd = sd(sum_branches, na.rm = TRUE),
    ymax = max(sum_branches, na.rm = TRUE)
  )

ggplot(mutate(vimentin_UV, condition_clean = fct_relevel(condition_clean, "bj p12.5 -uv vim", "bj p12.5 +uv vim")), aes(condition_clean, sum_branches, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Condition", y = "Number of branches /cell") +
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
      y = 24000,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = vimentin_UV, comparisons = list(c("bj p12.5 -uv vim","bj p12.5 +uv vim")), annotations = "ns", map_signif_level = TRUE, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/vimUVSumBranch.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")


# boxplot for vimentin_UV mean_branches
summary_stats <- vimentin_UV %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p12.5 -uv vim", "bj p12.5 +uv vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_branches, na.rm = TRUE),
    sd = sd(mean_branches, na.rm = TRUE),
    ymax = max(mean_branches, na.rm = TRUE)
  )

ggplot(mutate(vimentin_UV, condition_clean = fct_relevel(condition_clean, "bj p12.5 -uv vim", "bj p12.5 +uv vim")), aes(condition_clean, mean_branches, fill = condition_clean, color = condition_clean)) + 
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
      y = 0.89,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = vimentin_UV, comparisons = list(c("bj p12.5 -uv vim","bj p12.5 +uv vim")), annotations = "ns", map_signif_level = TRUE, y_position = 0.9, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/vimUVMeanBranch.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")


# boxplot for vimentin_etop sum_length_2d
summary_stats <- vimentin_etop %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p8.5 vim", "bj p9.5 etop vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(sum_length_2d, na.rm = TRUE),
    sd = sd(sum_length_2d, na.rm = TRUE),
    ymax = max(sum_length_2d, na.rm = TRUE)
  )

ggplot(mutate(vimentin_etop, condition_clean = fct_relevel(condition_clean, "bj p8.5 vim", "bj p9.5 etop vim")), aes(condition_clean, sum_length_2d, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Condition", y = "Total 2D filament length per cell (μm)") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('p8.5 -etop', 'p9.5 +etop')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(10000)) +
  geom_text(
    data = summary_stats,
    aes(
      x = condition_clean,
      y = 130000,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = vimentin_etop, comparisons = list(c("bj p8.5 vim","bj p9.5 etop vim")), annotations = "***", map_signif_level = TRUE, vjust = 0.7, y_position = , color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/vimEtopSum2d.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")


# boxplot for vimentin_etop mean_length_2d
summary_stats <- vimentin_etop %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p8.5 vim", "bj p9.5 etop vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_length_2d, na.rm = TRUE),
    sd = sd(mean_length_2d, na.rm = TRUE),
    ymax = max(mean_length_2d, na.rm = TRUE)
  )

ggplot(mutate(vimentin_etop, condition_clean = fct_relevel(condition_clean, "bj p8.5 vim", "bj p9.5 etop vim")), aes(condition_clean, mean_length_2d, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Condition", y = "Mean 2D filament length per cell (μm)") +
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
      y = 0.96,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = vimentin_etop, comparisons = list(c("bj p8.5 vim","bj p9.5 etop vim")), annotations = "*", map_signif_level = TRUE, vjust = 0.7, y_position = 0.97, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/vimEtopMean2d.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")


# boxplot for vimentin_etop sum_length_3d
summary_stats <- vimentin_etop %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p8.5 vim", "bj p9.5 etop vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(sum_length_3d, na.rm = TRUE),
    sd = sd(sum_length_3d, na.rm = TRUE),
    ymax = max(sum_length_3d, na.rm = TRUE)
  )

ggplot(mutate(vimentin_etop, condition_clean = fct_relevel(condition_clean, "bj p8.5 vim", "bj p9.5 etop vim")), aes(condition_clean, sum_length_3d, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Condition", y = "Total 3D filament length per cell (μm)") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('p8.5 -etop', 'p9.5 +etop')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(50000)) +
  geom_text(
    data = summary_stats,
    aes(
      x = condition_clean,
      y = 1500000,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = vimentin_etop, comparisons = list(c("bj p8.5 vim","bj p9.5 etop vim")), annotations = "*", map_signif_level = TRUE, vjust = 0.7, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/vimEtopSum3d.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")


# boxplot for vimentin_etop mean_length_3d
summary_stats <- vimentin_etop %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p8.5 vim", "bj p9.5 etop vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_length_3d, na.rm = TRUE),
    sd = sd(mean_length_3d, na.rm = TRUE),
    ymax = max(mean_length_3d, na.rm = TRUE)
  )

ggplot(mutate(vimentin_etop, condition_clean = fct_relevel(condition_clean, "bj p8.5 vim", "bj p9.5 etop vim")), aes(condition_clean, mean_length_3d, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Condition", y = "Mean 3D filament length per cell (μm)") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('p8.5 -etop', 'p9.5 +etop')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(0.2)) +
  geom_text(
    data = summary_stats,
    aes(
      x = condition_clean,
      y = 11.7,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = vimentin_etop, comparisons = list(c("bj p8.5 vim","bj p9.5 etop vim")), annotations = "ns", map_signif_level = TRUE, y_position = 12, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/vimEtopMean3d.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")


# boxplot for vimentin_etop mean_avg_angle
summary_stats <- vimentin_etop %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p8.5 vim", "bj p9.5 etop vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_avg_angle, na.rm = TRUE),
    sd = sd(mean_avg_angle, na.rm = TRUE),
    ymax = max(mean_avg_angle, na.rm = TRUE)
  )

ggplot(mutate(vimentin_etop, condition_clean = fct_relevel(condition_clean, "bj p8.5 vim", "bj p9.5 etop vim")), aes(condition_clean, mean_avg_angle, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Passage", y = "Mean average filament angle (deg)") +
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
      y = 83.8,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = vimentin_etop, comparisons = list(c("bj p8.5 vim", "bj p9.5 etop vim")), map_signif_level = TRUE, annotations = "ns", color = 'black') 
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/vimEtopMeanAvgAng.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")


# boxplot for vimentin_etop mean_deviation
summary_stats <- vimentin_etop %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p8.5 vim", "bj p9.5 etop vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_deviation, na.rm = TRUE),
    sd = sd(mean_deviation, na.rm = TRUE),
    ymax = max(mean_deviation, na.rm = TRUE)
  )

ggplot(mutate(vimentin_etop, condition_clean = fct_relevel(condition_clean, "bj p8.5 vim", "bj p9.5 etop vim")), aes(condition_clean, mean_deviation, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Condition", y = "Mean filament deviation (μm)") +
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
      y = 2,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = vimentin_etop, comparisons = list(c("bj p8.5 vim", "bj p9.5 etop vim")), annotations = "ns", map_signif_level = TRUE, y_position = 2.08, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/vimEtopmeanDev.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")


# boxplot for vimentin_etop mean_acc_curv
summary_stats <- vimentin_etop %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p8.5 vim", "bj p9.5 etop vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_acc_curv, na.rm = TRUE),
    sd = sd(mean_acc_curv, na.rm = TRUE),
    ymax = max(mean_acc_curv, na.rm = TRUE)
  )

ggplot(mutate(vimentin_etop, condition_clean = fct_relevel(condition_clean, "bj p8.5 vim", "bj p9.5 etop vim")), aes(condition_clean, mean_acc_curv, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Condition", y = "Mean accumulated filament curvature /cell") +
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
      y = 4.24,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = vimentin_etop, comparisons = list(c("bj p8.5 vim","bj p9.5 etop vim")), annotations = "ns", map_signif_level = TRUE, y_position = 4.25, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/vimEtopMeanAccCurv.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")


# boxplot for vimentin_etop mean_net_curv
summary_stats <- vimentin_etop %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p8.5 vim", "bj p9.5 etop vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_net_curv, na.rm = TRUE),
    sd = sd(mean_net_curv, na.rm = TRUE),
    ymax = max(mean_net_curv, na.rm = TRUE)
  )

ggplot(mutate(vimentin_etop, condition_clean = fct_relevel(condition_clean, "bj p8.5 vim", "bj p9.5 etop vim")), aes(condition_clean, mean_net_curv, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Condition", y = "Mean end-to-end filament curvature /cell") +
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
      y = 1.55,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = vimentin_etop, comparisons = list(c("bj p8.5 vim","bj p9.5 etop vim")), annotations = "ns", map_signif_level = TRUE, y_position = 1.56, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/vimEtopMeanNetCurv.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")


# boxplot for vimentin_etop filament_count
summary_stats <- vimentin_etop %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p8.5 vim", "bj p9.5 etop vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(filament_count, na.rm = TRUE),
    sd = sd(filament_count, na.rm = TRUE),
    ymax = max(filament_count, na.rm = TRUE)
  )

ggplot(mutate(vimentin_etop, condition_clean = fct_relevel(condition_clean, "bj p8.5 vim", "bj p9.5 etop vim")), aes(condition_clean, filament_count, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Condition", y = "Number of filaments /cell") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('p8.5 -etop', 'p9.5 +etop')) +
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
  geom_signif(data = vimentin_etop, comparisons = list(c("bj p8.5 vim","bj p9.5 etop vim")), annotations = "***", map_signif_level = TRUE, vjust = 0.7, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/vimEtopFilCount.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")


# boxplot for vimentin_etop sum_branches
summary_stats <- vimentin_etop %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p8.5 vim", "bj p9.5 etop vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(sum_branches, na.rm = TRUE),
    sd = sd(sum_branches, na.rm = TRUE),
    ymax = max(sum_branches, na.rm = TRUE)
  )

ggplot(mutate(vimentin_etop, condition_clean = fct_relevel(condition_clean, "bj p8.5 vim", "bj p9.5 etop vim")), aes(condition_clean, sum_branches, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Condition", y = "Number of branches /cell") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('p8.5 -etop', 'p9.5 +etop')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(10000)) +
  geom_text(
    data = summary_stats,
    aes(
      x = condition_clean,
      y = 140000,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = vimentin_etop, comparisons = list(c("bj p8.5 vim","bj p9.5 etop vim")), annotations = "***", map_signif_level = TRUE, vjust = 0.7, y_position = 148000, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/vimEtopSumBranch.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")


# boxplot for vimentin_etop mean_branches
summary_stats <- vimentin_etop %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p8.5 vim", "bj p9.5 etop vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_branches, na.rm = TRUE),
    sd = sd(mean_branches, na.rm = TRUE),
    ymax = max(mean_branches, na.rm = TRUE)
  )

ggplot(mutate(vimentin_etop, condition_clean = fct_relevel(condition_clean, "bj p8.5 vim", "bj p9.5 etop vim")), aes(condition_clean, mean_branches, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Condition", y = "Mean branches per filament /cell") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('p8.5 -etop', 'p9.5 +etop')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(0.02)) +
  geom_text(
    data = summary_stats,
    aes(
      x = condition_clean,
      y = 0.88,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = vimentin_etop, comparisons = list(c("bj p8.5 vim","bj p9.5 etop vim")), annotations = "**", map_signif_level = TRUE, vjust = 0.7, y_position = 0.9, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/vimEtopMeanBranch.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")

# doplot for vimentin_etop sum_length_3d vs sum_branches
ggplot(data = vimentin_etop, aes(x = sum_length_3d, y = sum_branches, fill = condition_clean, color = condition_clean)) + 
  geom_point(size = 0.8) + 
  labs(x = "Total 3D filament length (μm)", y = "Total branches /cell") +
  theme_bw() +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE), x = guide_axis(minor.ticks = TRUE)) +
  scale_y_continuous(minor_breaks = scales::breaks_width(5000)) +
  scale_x_continuous(minor_breaks = scales::breaks_width(100000)) 
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/vimPassSum2D.tiff", device = "tiff", dpi = 300, width = 177, height = 150, units = "mm")


# boxplot for vimentin_primary sum_length_2d
summary_stats <- vimentin_primary %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "primary fibro 18yo vim", "primary fibro 29yo vim", "primary fibro 32yo vim", "primary fibro 59yo vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(sum_length_2d, na.rm = TRUE),
    sd = sd(sum_length_2d, na.rm = TRUE),
    ymax = max(sum_length_2d, na.rm = TRUE)
  )

ggplot(mutate(vimentin_primary, condition_clean = fct_relevel(condition_clean, "primary fibro 18yo vim", "primary fibro 29yo vim", "primary fibro 32yo vim", "primary fibro 59yo vim")), aes(condition_clean, sum_length_2d, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Age", y = "Total 2D filament length per cell (μm)") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('18 y/o', '29 y/o', '32 y/o', '59 y/o')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(25000)) +
  geom_text(
    data = summary_stats,
    aes(
      x = condition_clean,
      y = 210000,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = vimentin_primary, comparisons = list(c("primary fibro 18yo vim", "primary fibro 29yo vim")), annotations = "**", map_signif_level = TRUE, y_position = 230000, vjust = 0.7, color = 'black') +
  geom_signif(data = vimentin_primary, comparisons = list(c("primary fibro 18yo vim", "primary fibro 32yo vim")), annotations = "***", map_signif_level = TRUE, y_position = 245000, vjust = 0.7, color = 'black') +
  geom_signif(data = vimentin_primary, comparisons = list(c("primary fibro 18yo vim", "primary fibro 59yo vim")), annotations = "***", map_signif_level = TRUE, y_position = 260000, vjust = 0.7, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/vimPrimSum2D.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")


# boxplot for vimentin_primary mean_length_2d
summary_stats <- vimentin_primary %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "primary fibro 18yo vim", "primary fibro 29yo vim", "primary fibro 32yo vim", "primary fibro 59yo vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_length_2d, na.rm = TRUE),
    sd = sd(mean_length_2d, na.rm = TRUE),
    ymax = max(mean_length_2d, na.rm = TRUE)
  )

ggplot(mutate(vimentin_primary, condition_clean = fct_relevel(condition_clean, "primary fibro 18yo vim", "primary fibro 29yo vim", "primary fibro 32yo vim", "primary fibro 59yo vim")), aes(condition_clean, mean_length_2d, fill = condition_clean, color = condition_clean)) + 
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
      y = 0.96,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = vimentin_primary, comparisons = list(c("primary fibro 18yo vim", "primary fibro 59yo vim")), annotations = "ns", map_signif_level = TRUE, y_position = 0.97, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/vimPrimMean2D.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")


# boxplot for vimentin_primary sum_length_3d
summary_stats <- vimentin_primary %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "primary fibro 18yo vim", "primary fibro 29yo vim", "primary fibro 32yo vim", "primary fibro 59yo vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(sum_length_3d, na.rm = TRUE),
    sd = sd(sum_length_3d, na.rm = TRUE),
    ymax = max(sum_length_3d, na.rm = TRUE)
  )

ggplot(mutate(vimentin_primary, condition_clean = fct_relevel(condition_clean, "primary fibro 18yo vim", "primary fibro 29yo vim", "primary fibro 32yo vim", "primary fibro 59yo vim")), aes(condition_clean, sum_length_3d, fill = condition_clean, color = condition_clean)) + 
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
      y = 3000000,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = vimentin_primary, comparisons = list(c("primary fibro 18yo vim", "primary fibro 29yo vim")), annotations = "***", map_signif_level = TRUE, y_position = 3200000, vjust = 0.7, color = 'black') +
  geom_signif(data = vimentin_primary, comparisons = list(c("primary fibro 18yo vim", "primary fibro 32yo vim")), annotations = "***", map_signif_level = TRUE, y_position = 3400000, vjust = 0.7, color = 'black') +
  geom_signif(data = vimentin_primary, comparisons = list(c("primary fibro 18yo vim", "primary fibro 59yo vim")), annotations = "***", map_signif_level = TRUE, y_position = 3600000, vjust = 0.7, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/vimPrimSum3D.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")


# boxplot for vimentin_primary mean_length_3d
summary_stats <- vimentin_primary %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "primary fibro 18yo vim", "primary fibro 29yo vim", "primary fibro 32yo vim", "primary fibro 59yo vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_length_3d, na.rm = TRUE),
    sd = sd(mean_length_3d, na.rm = TRUE),
    ymax = max(mean_length_3d, na.rm = TRUE)
  )

ggplot(mutate(vimentin_primary, condition_clean = fct_relevel(condition_clean, "primary fibro 18yo vim", "primary fibro 29yo vim", "primary fibro 32yo vim", "primary fibro 59yo vim")), aes(condition_clean, mean_length_3d, fill = condition_clean, color = condition_clean)) + 
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
      y = 11,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = vimentin_primary, comparisons = list(c("primary fibro 18yo vim", "primary fibro 59yo vim")), annotations = "ns", map_signif_level = TRUE, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/vimPrimmean3d.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")


# boxplot for vimentin_primary mean_avg_angle
summary_stats <- vimentin_primary %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "primary fibro 18yo vim", "primary fibro 29yo vim", "primary fibro 32yo vim", "primary fibro 59yo vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_avg_angle, na.rm = TRUE),
    sd = sd(mean_avg_angle, na.rm = TRUE),
    ymax = max(mean_avg_angle, na.rm = TRUE)
  )

ggplot(mutate(vimentin_primary, condition_clean = fct_relevel(condition_clean, "primary fibro 18yo vim", "primary fibro 29yo vim", "primary fibro 32yo vim", "primary fibro 59yo vim")), aes(condition_clean, mean_avg_angle, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Age", y = "Mean average filament angle (deg)") +
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
      y = 86,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = vimentin_primary, comparisons = list(c("primary fibro 18yo vim", "primary fibro 59yo vim")), annotations = "ns", map_signif_level = TRUE, y_position = 87, color = 'black') 
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/vimPrimmeanAvgAng.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")


# boxplot for vimentin_primary mean_acc_curv
summary_stats <- vimentin_primary %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "primary fibro 18yo vim", "primary fibro 29yo vim", "primary fibro 32yo vim", "primary fibro 59yo vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_acc_curv, na.rm = TRUE),
    sd = sd(mean_acc_curv, na.rm = TRUE),
    ymax = max(mean_acc_curv, na.rm = TRUE)
  )

ggplot(mutate(vimentin_primary, condition_clean = fct_relevel(condition_clean, "primary fibro 18yo vim", "primary fibro 29yo vim", "primary fibro 32yo vim", "primary fibro 59yo vim")), aes(condition_clean, mean_acc_curv, fill = condition_clean, color = condition_clean)) + 
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
      y = 4.23,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = vimentin_primary, comparisons = list(c("primary fibro 18yo vim", "primary fibro 59yo vim")), annotations = "ns", map_signif_level = TRUE, y_position = 4.24, color = 'black') 
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/vimPrimmeanAccCurv.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")


# boxplot for vimentin_primary mean_net_curv
summary_stats <- vimentin_primary %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "primary fibro 18yo vim", "primary fibro 29yo vim", "primary fibro 32yo vim", "primary fibro 59yo vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_net_curv, na.rm = TRUE),
    sd = sd(mean_net_curv, na.rm = TRUE),
    ymax = max(mean_net_curv, na.rm = TRUE)
  )

ggplot(mutate(vimentin_primary, condition_clean = fct_relevel(condition_clean, "primary fibro 18yo vim", "primary fibro 29yo vim", "primary fibro 32yo vim", "primary fibro 59yo vim")), aes(condition_clean, mean_net_curv, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Age", y = "Mean end-to-end filament curvature /cell") +
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
      y = 1.485,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = vimentin_primary, comparisons = list(c("primary fibro 18yo vim", "primary fibro 59yo vim")), annotations = "ns", map_signif_level = TRUE, y_position = 1.51, color = 'black') 
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/vimPrimmeanNetCurv.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")


# boxplot for vimentin_primary filament_count
summary_stats <- vimentin_primary %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "primary fibro 18yo vim", "primary fibro 29yo vim", "primary fibro 32yo vim", "primary fibro 59yo vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(filament_count, na.rm = TRUE),
    sd = sd(filament_count, na.rm = TRUE),
    ymax = max(filament_count, na.rm = TRUE)
  )

ggplot(mutate(vimentin_primary, condition_clean = fct_relevel(condition_clean, "primary fibro 18yo vim", "primary fibro 29yo vim", "primary fibro 32yo vim", "primary fibro 59yo vim")), aes(condition_clean, filament_count, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size= 0.7) + 
  labs(x = "Age", y = "Filament count /cell") +
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
      y = 280000,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = vimentin_primary, comparisons = list(c("primary fibro 18yo vim", "primary fibro 59yo vim")), annotations = "***", map_signif_level = TRUE, vjust = 0.7, y_position = 348000, color = 'black') +
  geom_signif(data = vimentin_primary, comparisons = list(c("primary fibro 18yo vim", "primary fibro 32yo vim")), annotations = "***", map_signif_level = TRUE, vjust = 0.7, y_position = 334000, color = 'black') +
  geom_signif(data = vimentin_primary, comparisons = list(c("primary fibro 18yo vim", "primary fibro 29yo vim")), annotations = "**", map_signif_level = TRUE, vjust = 0.7, y_position = 320000, color = 'black') 
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/vimPrimFilCount.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")


# boxplot for vimentin_primary sum_branches
summary_stats <- vimentin_primary %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "primary fibro 18yo vim", "primary fibro 29yo vim", "primary fibro 32yo vim", "primary fibro 59yo vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(sum_branches, na.rm = TRUE),
    sd = sd(sum_branches, na.rm = TRUE),
    ymax = max(sum_branches, na.rm = TRUE)
  )

ggplot(mutate(vimentin_primary, condition_clean = fct_relevel(condition_clean, "primary fibro 18yo vim", "primary fibro 29yo vim", "primary fibro 32yo vim", "primary fibro 59yo vim")), aes(condition_clean, sum_branches, fill = condition_clean, color = condition_clean)) + 
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
      y = 210000,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = vimentin_primary, comparisons = list(c("primary fibro 18yo vim", "primary fibro 59yo vim")), annotations = "***", map_signif_level = TRUE, vjust = 0.7, y_position = 260000, color = 'black') +
  geom_signif(data = vimentin_primary, comparisons = list(c("primary fibro 18yo vim", "primary fibro 32yo vim")), annotations = "***", map_signif_level = TRUE, vjust = 0.7, y_position = 250000, color = 'black') +
  geom_signif(data = vimentin_primary, comparisons = list(c("primary fibro 18yo vim", "primary fibro 29yo vim")), annotations = "**", map_signif_level = TRUE, vjust = 0.7, y_position = 240000, color = 'black') 
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/vimPrimBranchCount.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")


# boxplot for vimentin_primary mean_branches
summary_stats <- vimentin_primary %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "primary fibro 18yo vim", "primary fibro 29yo vim", "primary fibro 32yo vim", "primary fibro 59yo vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_branches, na.rm = TRUE),
    sd = sd(mean_branches, na.rm = TRUE),
    ymax = max(mean_branches, na.rm = TRUE)
  )

ggplot(mutate(vimentin_primary, condition_clean = fct_relevel(condition_clean, "primary fibro 18yo vim", "primary fibro 29yo vim", "primary fibro 32yo vim", "primary fibro 59yo vim")), aes(condition_clean, mean_branches, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Age", y = "Mean branches per filament /cell") +
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
      y = 0.85,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = vimentin_primary, comparisons = list(c("primary fibro 18yo vim", "primary fibro 59yo vim")), annotations = "ns", map_signif_level = TRUE, y_position = 0.86, color = 'black') 
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/vimPrimMeanBranch.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")


# boxplot for vimentin_primary mean_deviation
summary_stats <- vimentin_primary %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "primary fibro 18yo vim", "primary fibro 29yo vim", "primary fibro 32yo vim", "primary fibro 59yo vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_deviation, na.rm = TRUE),
    sd = sd(mean_deviation, na.rm = TRUE),
    ymax = max(mean_deviation, na.rm = TRUE)
  )

ggplot(mutate(vimentin_primary, condition_clean = fct_relevel(condition_clean, "primary fibro 18yo vim", "primary fibro 29yo vim", "primary fibro 32yo vim", "primary fibro 59yo vim")), aes(condition_clean, mean_deviation, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Age", y = "Mean filament deviation (μm)") +
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
      y = 1.9,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = vimentin_primary, comparisons = list(c("primary fibro 18yo vim", "primary fibro 59yo vim")), annotations = "ns", map_signif_level = TRUE, color = 'black') 
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/vimPrimMeanDev.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")


# doplot for vimentin_primary sum_length_3d vs sum_branches
ggplot(data = vimentin_primary, aes(x = sum_length_3d, y = sum_branches, fill = condition_clean, color = condition_clean)) + 
  geom_point(size = 0.8) + 
  labs(x = "Total 3D filament length (μm)", y = "Total branches /cell") +
  theme_bw() +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, y = guide_axis(minor.ticks = TRUE), x = guide_axis(minor.ticks = TRUE)) +
  scale_y_continuous(minor_breaks = scales::breaks_width(5000)) 
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/vimPrimSum2D.tiff", device = "tiff", dpi = 300, width = 177, height = 150, units = "mm")
