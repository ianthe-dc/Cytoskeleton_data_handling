library(tidyverse)
library(lme4)
library(car)
library(emmeans)
library(ggprism)
library(forcats)
library(dplyr)
library(ggsignif)
library(DHARMa)

actin_df <- read.csv("actin cell level.csv")
actin_BJ <- filter(actin_df, condition_clean %in% c("bj p8.5 act", "bj p9.5 etop vim", "bj p11.5 act", "bj p12.5 -uv act", "bj p12.5 +uv act", "bj p26.5 act", "bj p31.5 act", "bj p36.5 act", "bj p38.5 act", "bj p40.5 act"))
actin_primary <- actin_df[259:310,]
actin_UV <- filter(actin_df, condition_clean %in% c("bj p12.5 -uv act", "bj p12.5 +uv act"))
actin_etop <- filter(actin_df, condition_clean %in% c("bj p8.5 act", "bj p9.5 etop vim"))
actin_pass <- filter(actin_df, condition_clean %in% c("bj p11.5 act", "bj p26.5 act", "bj p31.5 act", "bj p36.5 act", "bj p38.5 act", "bj p40.5 act"))

#stats for all actin BJ datasets
model <- lmer(sum_length_2d ~ condition_clean + (1|condition_clean:biorep) +(1|condition_clean:biorep:techrep), actin_BJ)
summary(model)
Anova(model)
em1 <- emmeans(model, ~condition_clean)
pairs(em1)

#stats for UV only
model_UV <- lmer(sum_length_2d ~ condition_clean + (1|condition_clean:biorep) +(1|condition_clean:biorep:techrep), actin_UV)
summary(model_UV)
Anova(model_UV)
em_UV <- emmeans(model_UV, ~condition_clean)
pairs(em_UV)


#stats for passage only sum_length_2d
model_pass <- lmer(sum_length_2d ~ condition_clean + (1|condition_clean:biorep:techrep), actin_pass)
summary(model_pass)
Anova(model_pass)
em_pass <- emmeans(model_pass, ~condition_clean)
pairs(em_pass)

simulationOutput <- simulateResiduals(fittedModel = model_pass, plot = T)

kruskal <- kruskal.test(sum_length_2d ~ condition_clean, data = actin_pass)
kruskal

#stats for passage only mean_length_2d
model2_pass <- lmer(mean_length_2d ~ condition_clean + (1|condition_clean:biorep:techrep), actin_pass)
summary(model2_pass)
Anova(model2_pass)
em_pass2 <- emmeans(model2_pass, ~condition_clean)
pairs(em_pass2)

simulationOutput <- simulateResiduals(fittedModel = model2_pass, plot = T)

#stats for passage only sum_length_3d
model3_pass <- lmer(sum_length_3d ~ condition_clean + (1|condition_clean:biorep:techrep), actin_pass)
summary(model3_pass)
Anova(model3_pass)
em_pass3 <- emmeans(model3_pass, ~condition_clean)
pairs(em_pass3)

simulationOutput <- simulateResiduals(fittedModel = model3_pass, plot = T)

kruskal <- kruskal.test(sum_length_3d ~ condition_clean, data = actin_pass)
kruskal
pairwise.wilcox.test(actin_pass$sum_length_3d, actin_pass$condition_clean, p.adjust.method = "BH")

#stats for passage only mean_length_3d
model4_pass <- lm(mean_length_3d ~ condition_clean + techrep, actin_pass)
summary(model4_pass)
Anova(model4_pass)
em_pass4 <- emmeans(model4_pass, ~condition_clean)
pairs(em_pass4)

simulationOutput <- simulateResiduals(fittedModel = model4_pass, plot = T)

#stats for passage only mean_avg_angle
model5_pass <- lm(mean_avg_angle ~ condition_clean + techrep, actin_pass)
summary(model5_pass)
Anova(model5_pass)
em_pass5 <- emmeans(model5_pass, ~condition_clean)
pairs(em_pass5)

simulationOutput <- simulateResiduals(fittedModel = model5_pass, plot = T)

kruskal <- kruskal.test(mean_avg_angle ~ condition_clean, data = actin_pass)
kruskal
pairwise.wilcox.test(actin_pass$mean_avg_angle, actin_pass$condition_clean, p.adjust.method = "BH")

#stats for passage only mean_acc_curv
model6_pass <- lmer(mean_acc_curv ~ condition_clean + (1|condition_clean:biorep:techrep), actin_pass)
summary(model6_pass)
Anova(model6_pass)
em_pass6 <- emmeans(model6_pass, ~condition_clean)
pairs(em_pass6)

simulationOutput <- simulateResiduals(fittedModel = model6_pass, plot = T)

#stats for passage only mean_net_curv
model6_pass <- lmer(mean_net_curv ~ condition_clean + (1|condition_clean:biorep:techrep), actin_pass)
summary(model6_pass)
Anova(model6_pass)
em_pass6 <- emmeans(model6_pass, ~condition_clean)
pairs(em_pass6)

simulationOutput <- simulateResiduals(fittedModel = model6_pass, plot = T)

#stats for passage only filament_count
model7_pass <- lmer(filament_count ~ condition_clean + (1|condition_clean:biorep:techrep), actin_pass)
summary(model7_pass)
Anova(model7_pass)
em_pass7 <- emmeans(model7_pass, ~condition_clean)
pairs(em_pass7)

simulationOutput <- simulateResiduals(fittedModel = model7_pass, plot = T)

kruskal <- kruskal.test(filament_count ~ condition_clean, data = actin_pass)
kruskal
pairwise.wilcox.test(actin_pass$filament_count, actin_pass$condition_clean, p.adjust.method = "BH")

#stats for passage only sum_branches
model8_pass <- lmer(sum_branches ~ condition_clean + (1|condition_clean:biorep:techrep), actin_pass)
summary(model8_pass)
Anova(model8_pass)
em_pass8 <- emmeans(model8_pass, ~condition_clean)
pairs(em_pass8)

simulationOutput <- simulateResiduals(fittedModel = model8_pass, plot = T)

kruskal <- kruskal.test(sum_branches ~ condition_clean, data = actin_pass)
kruskal
pairwise.wilcox.test(actin_pass$sum_branches, actin_pass$condition_clean, p.adjust.method = "BH")

#stats for passage only mean_branches
model9_pass <- lmer(mean_branches ~ condition_clean + (1|condition_clean:biorep), actin_pass)
summary(model9_pass)
Anova(model9_pass)
em_pass9 <- emmeans(model9_pass, ~condition_clean)
pairs(em_pass9)

simulationOutput <- simulateResiduals(fittedModel = model9_pass, plot = T)

#stats for passage only mean_deviation
model10_pass <- lmer(mean_deviation ~ condition_clean + (1|condition_clean:biorep:techrep), actin_pass)
summary(model10_pass)
Anova(model10_pass)
em_pass10 <- emmeans(model10_pass, ~condition_clean)
pairs(em_pass10)

simulationOutput <- simulateResiduals(fittedModel = model10_pass, plot = T)

#stats for UV 
model_UV <- lmer(sum_length_3d ~ condition_clean + (1|condition_clean:biorep:techrep), actin_UV)
summary(model_UV)
Anova(model_UV)
em_UV <- emmeans(model_UV, ~condition_clean)
pairs(em_UV)

simulationOutput <- simulateResiduals(fittedModel = model_UV, plot = T)


#stats for etop only
model_etop <- lmer(mean_branches ~ condition_clean + (1|condition_clean:biorep:techrep), actin_etop)
summary(model_etop)
Anova(model_etop)
em_etop <- emmeans(model_etop, ~condition_clean)
pairs(em_etop)

simulationOutput <- simulateResiduals(fittedModel = model_etop, plot = T)

kruskal <- kruskal.test(sum_branches ~ condition_clean, data = actin_etop)
kruskal

#stats for actin primary
model13 <- lmer(sum_length_2d ~ condition_clean + (1|condition_clean:biorep), actin_primary)
summary(model13)
Anova(model13)
em13 <- emmeans(model13, ~condition_clean)
pairs(em13)

simulationOutput <- simulateResiduals(fittedModel = model13, plot = T)

kruskal <- kruskal.test(sum_length_2d ~ condition_clean, data = actin_etop)
kruskal
pairwise.wilcox.test(actin_primary$sum_length_2d, actin_primary$condition_clean, p.adjust.method = "BH")

# boxplot for all actin sum_length_2D
ggplot(mutate(actin_BJ, condition_clean = fct_relevel(condition_clean, "bj p8.5 act", "bj p9.5 etop vim", "bj p11.5 act", "bj p12.5 -uv act", "bj p12.5 +uv act", "bj p26.5 act", "bj p31.5 act", "bj p36.5 act", "bj p38.5 act", "bj p40.5 act")), aes(condition_clean, sum_length_2d, fill = condition_clean, color = condition_clean)) + 
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

# boxplot for all actin mean_length_2D
ggplot(mutate(actin_BJ, condition_clean = fct_relevel(condition_clean, "bj p8.5 act", "bj p9.5 etop vim", "bj p11.5 act", "bj p12.5 -uv act", "bj p12.5 +uv act", "bj p26.5 act", "bj p31.5 act", "bj p36.5 act", "bj p38.5 act", "bj p40.5 act")), aes(condition_clean, mean_length_2d, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(size = 0.5, position = position_jitterdodge()) +
  labs(x = "Condition", y = "Mean 2D filament length per cell (μm)") +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('p8.5', 'p9.5 +etop', 'p11.5', 'p12.5 -UV', 'p12.5 +UV', 'p26.5', 'p31.5', 'p36.5', 'p38.5', 'p40.5')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(0.1))
ggsave("vimBJmean2D.tiff", device = "tiff", dpi = 300, width = 177, height = 150, units = "mm")

# boxplot for all actin sum_length_3D
ggplot(mutate(actin_BJ, condition_clean = fct_relevel(condition_clean, "bj p8.5 act", "bj p9.5 etop vim", "bj p11.5 act", "bj p12.5 -uv act", "bj p12.5 +uv act", "bj p26.5 act", "bj p31.5 act", "bj p36.5 act", "bj p38.5 act", "bj p40.5 act")), aes(condition_clean, sum_length_3d, fill = condition_clean, color = condition_clean)) + 
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

# boxplot for all actin mean_length_3D
ggplot(mutate(actin_BJ, condition_clean = fct_relevel(condition_clean, "bj p8.5 act", "bj p9.5 etop vim", "bj p11.5 act", "bj p12.5 -uv act", "bj p12.5 +uv act", "bj p26.5 act", "bj p31.5 act", "bj p36.5 act", "bj p38.5 act", "bj p40.5 act")), aes(condition_clean, mean_length_3d, fill = condition_clean, color = condition_clean)) + 
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


# boxplot for actin_pass sum_length_2d
summary_stats <- actin_pass %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p11.5 act","bj p26.5 act","bj p31.5 act",
                                       "bj p36.5 act","bj p38.5 act","bj p40.5 act")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(sum_length_2d, na.rm = TRUE),
    sd = sd(sum_length_2d, na.rm = TRUE),
    ymax = max(sum_length_2d, na.rm = TRUE)
  )

ggplot(mutate(actin_pass, condition_clean = fct_relevel(condition_clean, "bj p11.5 act", "bj p26.5 act", "bj p31.5 act", "bj p36.5 act", "bj p38.5 act", "bj p40.5 act")), aes(condition_clean, sum_length_2d, fill = condition_clean, color = condition_clean)) + 
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
      y = 66000,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 2.5
  ) +
  geom_signif(data = actin_pass, comparisons = list(c("bj p11.5 act","bj p40.5 act")), map_signif_level = TRUE, annotations = "ns", y_position = 68000, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/actPassSum2D.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")


# boxplot for actin_pass mean_length_2d
summary_stats <- actin_pass %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p11.5 act","bj p26.5 act","bj p31.5 act",
                                       "bj p36.5 act","bj p38.5 act","bj p40.5 act")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_length_2d, na.rm = TRUE),
    sd = sd(mean_length_2d, na.rm = TRUE),
    ymax = max(mean_length_2d, na.rm = TRUE)
  )

ggplot(mutate(actin_pass, condition_clean = fct_relevel(condition_clean, "bj p11.5 act", "bj p26.5 act", "bj p31.5 act", "bj p36.5 act", "bj p38.5 act", "bj p40.5 act")), aes(condition_clean, mean_length_2d, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Passage", y = "Mean 2D filament length per cell (μm)") +
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
      y = 1.85,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 2.9
  ) +
  geom_signif(data = actin_pass, comparisons = list(c("bj p11.5 act","bj p36.5 act")), map_signif_level = TRUE, annotations = "*", y_position = 1.93, vjust = 0.7, color = 'black') +
  geom_signif(data = actin_pass, comparisons = list(c("bj p11.5 act","bj p38.5 act")), map_signif_level = TRUE, annotations = "*", y_position = 1.98, vjust = 0.7, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/actPassMean2D.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")


# boxplot for actin_pass sum_length_3d
summary_stats <- actin_pass %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p11.5 act","bj p26.5 act","bj p31.5 act",
                                       "bj p36.5 act","bj p38.5 act","bj p40.5 act")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(sum_length_3d, na.rm = TRUE),
    sd = sd(sum_length_3d, na.rm = TRUE),
    ymax = max(sum_length_3d, na.rm = TRUE)
  )

ggplot(mutate(actin_pass, condition_clean = fct_relevel(condition_clean, "bj p11.5 act", "bj p26.5 act", "bj p31.5 act", "bj p36.5 act", "bj p38.5 act", "bj p40.5 act")), aes(condition_clean, sum_length_3d, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Passage", y = "Total 3D filament length per cell (μm)") +
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
      y = 500000,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 2.3
  ) +
  geom_signif(data = actin_pass, comparisons = list(c("bj p11.5 act","bj p31.5 act")), map_signif_level = TRUE, annotations = "*", y_position = 5.5e+05, vjust = 0.7, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/actPassSum3D.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")


# boxplot for actin_pass mean_length_3d
summary_stats <- actin_pass %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p11.5 act","bj p26.5 act","bj p31.5 act",
                                       "bj p36.5 act","bj p38.5 act","bj p40.5 act")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_length_3d, na.rm = TRUE),
    sd = sd(mean_length_3d, na.rm = TRUE),
    ymax = max(mean_length_3d, na.rm = TRUE)
  )

ggplot(mutate(actin_pass, condition_clean = fct_relevel(condition_clean, "bj p11.5 act", "bj p26.5 act", "bj p31.5 act", "bj p36.5 act", "bj p38.5 act", "bj p40.5 act")), aes(condition_clean, mean_length_3d, fill = condition_clean, color = condition_clean)) + 
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
      y = 9.4,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = actin_pass, comparisons = list(c("bj p11.5 act","bj p31.5 act")), map_signif_level = TRUE, annotations = "***", y_position = 10, color = 'black') +
  geom_signif(data = actin_pass, comparisons = list(c("bj p31.5 act","bj p36.5 act")), map_signif_level = TRUE, annotations = "***", y_position = 9.8, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/actPassMean3D.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")


# boxplot for actin_pass mean_avg_angle
summary_stats <- actin_pass %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p11.5 act","bj p26.5 act","bj p31.5 act",
                                       "bj p36.5 act","bj p38.5 act","bj p40.5 act")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_avg_angle, na.rm = TRUE),
    sd = sd(mean_avg_angle, na.rm = TRUE),
    ymax = max(mean_avg_angle, na.rm = TRUE)
  )

ggplot(mutate(actin_pass, condition_clean = fct_relevel(condition_clean, "bj p11.5 act", "bj p26.5 act", "bj p31.5 act", "bj p36.5 act", "bj p38.5 act", "bj p40.5 act")), aes(condition_clean, mean_avg_angle, fill = condition_clean, color = condition_clean)) + 
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
      y = 92.4,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = actin_pass, comparisons = list(c("bj p11.5 act","bj p36.5 act")), map_signif_level = TRUE, annotations = "**", y_position = 93.4, vjust = 0.7, color = 'black') +
  geom_signif(data = actin_pass, comparisons = list(c("bj p11.5 act","bj p38.5 act")), map_signif_level = TRUE, annotations = "**", y_position = 94.4, vjust = 0.7, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/actPassMeanAvgAng.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")


# boxplot for actin_pass mean_acc_curv
summary_stats <- actin_pass %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p11.5 act","bj p26.5 act","bj p31.5 act",
                                       "bj p36.5 act","bj p38.5 act","bj p40.5 act")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_acc_curv, na.rm = TRUE),
    sd = sd(mean_acc_curv, na.rm = TRUE),
    ymax = max(mean_acc_curv, na.rm = TRUE)
  )

ggplot(mutate(actin_pass, condition_clean = fct_relevel(condition_clean, "bj p11.5 act", "bj p26.5 act", "bj p31.5 act", "bj p36.5 act", "bj p38.5 act", "bj p40.5 act")), aes(condition_clean, mean_acc_curv, fill = condition_clean, color = condition_clean)) + 
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
      y = 4.26,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = actin_pass, comparisons = list(c("bj p11.5 act","bj p36.5 act")), map_signif_level = TRUE, annotations = "*", y_position = 4.28, vjust = 0.7, color = 'black') +
  geom_signif(data = actin_pass, comparisons = list(c("bj p11.5 act","bj p38.5 act")), map_signif_level = TRUE, annotations = "*", y_position = 4.3, vjust = 0.7, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/actPassMeanAccCurv.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")


# boxplot for actin_pass mean_net_curv
summary_stats <- actin_pass %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p11.5 act","bj p26.5 act","bj p31.5 act",
                                       "bj p36.5 act","bj p38.5 act","bj p40.5 act")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_net_curv, na.rm = TRUE),
    sd = sd(mean_net_curv, na.rm = TRUE),
    ymax = max(mean_net_curv, na.rm = TRUE)
  )

ggplot(mutate(actin_pass, condition_clean = fct_relevel(condition_clean, "bj p11.5 act", "bj p26.5 act", "bj p31.5 act", "bj p36.5 act", "bj p38.5 act", "bj p40.5 act")), aes(condition_clean, mean_net_curv, fill = condition_clean, color = condition_clean)) + 
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
      y = 1.615,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = actin_pass, comparisons = list(c("bj p11.5 act","bj p36.5 act")), map_signif_level = TRUE, annotations = "**", y_position = 1.63, vjust = 0.7, color = 'black') +
  geom_signif(data = actin_pass, comparisons = list(c("bj p11.5 act","bj p38.5 act")), map_signif_level = TRUE, annotations = "**", y_position = 1.65, vjust = 0.7, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/actPassMeanNetCurv.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")


# boxplot for actin_pass filament_count
summary_stats <- actin_pass %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p11.5 act","bj p26.5 act","bj p31.5 act",
                                       "bj p36.5 act","bj p38.5 act","bj p40.5 act")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(filament_count, na.rm = TRUE),
    sd = sd(filament_count, na.rm = TRUE),
    ymax = max(filament_count, na.rm = TRUE)
  )

ggplot(mutate(actin_pass, condition_clean = fct_relevel(condition_clean, "bj p11.5 act", "bj p26.5 act", "bj p31.5 act", "bj p36.5 act", "bj p38.5 act", "bj p40.5 act")), aes(condition_clean, filament_count, fill = condition_clean, color = condition_clean)) + 
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
      y = 83000,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 2.5
  ) +
  geom_signif(data = actin_pass, comparisons = list(c("bj p11.5 act","bj p40.5 act")), map_signif_level = TRUE, annotations = "ns", y_position = 85000, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/actPassFilCount.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")


# boxplot for actin_pass sum_branches
summary_stats <- actin_pass %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p11.5 act","bj p26.5 act","bj p31.5 act",
                                       "bj p36.5 act","bj p38.5 act","bj p40.5 act")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(sum_branches, na.rm = TRUE),
    sd = sd(sum_branches, na.rm = TRUE),
    ymax = max(sum_branches, na.rm = TRUE)
  )

ggplot(mutate(actin_pass, condition_clean = fct_relevel(condition_clean, "bj p11.5 act", "bj p26.5 act", "bj p31.5 act", "bj p36.5 act", "bj p38.5 act", "bj p40.5 act")), aes(condition_clean, sum_branches, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Passage", y = "Total filament branches /cell") +
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
      y = 58500,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 2.5
  ) +
  geom_signif(data = actin_pass, comparisons = list(c("bj p11.5 act","bj p40.5 act")), map_signif_level = TRUE, annotations = "ns", y_position = 60000, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/actPassBranchCount.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")


# boxplot for actin_pass mean_branches
summary_stats <- actin_pass %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p11.5 act","bj p26.5 act","bj p31.5 act",
                                       "bj p36.5 act","bj p38.5 act","bj p40.5 act")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_branches, na.rm = TRUE),
    sd = sd(mean_branches, na.rm = TRUE),
    ymax = max(mean_branches, na.rm = TRUE)
  )

ggplot(mutate(actin_pass, condition_clean = fct_relevel(condition_clean, "bj p11.5 act", "bj p26.5 act", "bj p31.5 act", "bj p36.5 act", "bj p38.5 act", "bj p40.5 act")), aes(condition_clean, mean_branches, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Passage", y = "Mean filament branches per cell") +
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
      y = 0.905,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = actin_pass, comparisons = list(c("bj p11.5 act","bj p36.5 act")), annotations = "**", map_signif_level = TRUE, y_position = 0.913, vjust = 0.7, color = 'black') +
  geom_signif(data = actin_pass, comparisons = list(c("bj p11.5 act","bj p40.5 act")), annotations = "**", map_signif_level = TRUE, y_position = 0.923, vjust = 0.7, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/actPassMeanBranch.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")


# boxplot for actin_pass mean_deviation
summary_stats <- actin_pass %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p11.5 act","bj p26.5 act","bj p31.5 act",
                                       "bj p36.5 act","bj p38.5 act","bj p40.5 act")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_deviation, na.rm = TRUE),
    sd = sd(mean_deviation, na.rm = TRUE),
    ymax = max(mean_deviation, na.rm = TRUE)
  )

ggplot(mutate(actin_pass, condition_clean = fct_relevel(condition_clean, "bj p11.5 act", "bj p26.5 act", "bj p31.5 act", "bj p36.5 act", "bj p38.5 act", "bj p40.5 act")), aes(condition_clean, mean_deviation, fill = condition_clean, color = condition_clean)) + 
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
      y = 1.49,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = actin_pass, comparisons = list(c("bj p11.5 act","bj p31.5 act")), annotations = "**", map_signif_level = TRUE, y_position = 1.56, vjust = 0.7, color = 'black') +
  geom_signif(data = actin_pass, comparisons = list(c("bj p31.5 act","bj p36.5 act")), annotations = "**", map_signif_level = TRUE, y_position = 1.52, vjust = 0.7, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/actPassMeanDev.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")

# doplot for actin_pass sum_length_3d vs sum_branches
ggplot(data = actin_pass, aes(x = sum_length_3d, y = sum_branches, fill = condition_clean, color = condition_clean)) + 
  geom_point(size = 0.8) + 
  labs(x = "Total 3D filament length (μm)", y = "Total branches /cell") +
  theme_bw() +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE), x = guide_axis(minor.ticks = TRUE)) +
  scale_y_continuous(minor_breaks = scales::breaks_width(5000)) 
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/vimPassSum2D.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")


# boxplot for actin_UV sum_length_2d
summary_stats <- actin_UV %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p12.5 -uv act", "bj p12.5 +uv act")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(sum_length_2d, na.rm = TRUE),
    sd = sd(sum_length_2d, na.rm = TRUE),
    ymax = max(sum_length_2d, na.rm = TRUE)
  )

ggplot(mutate(actin_UV, condition_clean = fct_relevel(condition_clean, "bj p12.5 -uv act", "bj p12.5 +uv act")), aes(condition_clean, sum_length_2d, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Condition", y = "Total 2D filament length (μm)") +
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
      y = 36000,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = actin_UV, comparisons = list(c("bj p12.5 -uv act", "bj p12.5 +uv act")), annotations = "*", map_signif_level = TRUE, y_position = 37500, vjust = 0.7, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/actUVsum2d.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")


# boxplot for actin_UV mean_length_2d
summary_stats <- actin_UV %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p12.5 -uv act", "bj p12.5 +uv act")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_length_2d, na.rm = TRUE),
    sd = sd(mean_length_2d, na.rm = TRUE),
    ymax = max(mean_length_2d, na.rm = TRUE)
  )

ggplot(mutate(actin_UV, condition_clean = fct_relevel(condition_clean, "bj p12.5 -uv act", "bj p12.5 +uv act")), aes(condition_clean, mean_length_2d, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Condition", y = "Mean 2D filament length (μm)") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('p12.5 -UV', 'p12.5 +UV')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(0.02)) +
  geom_text(
    data = summary_stats,
    aes(
      x = condition_clean,
      y = 1.45,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = actin_UV, comparisons = list(c("bj p12.5 -uv act", "bj p12.5 +uv act")), annotations = "ns", map_signif_level = TRUE, y_position = 1.48, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/actUVmean2d.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")

# boxplot for actin_UV sum_length_3d
summary_stats <- actin_UV %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p12.5 -uv act", "bj p12.5 +uv act")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(sum_length_3d, na.rm = TRUE),
    sd = sd(sum_length_3d, na.rm = TRUE),
    ymax = max(sum_length_3d, na.rm = TRUE)
  )

ggplot(mutate(actin_UV, condition_clean = fct_relevel(condition_clean, "bj p12.5 -uv act", "bj p12.5 +uv act")), aes(condition_clean, sum_length_3d, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Condition", y = "Total 3D filament length (μm)") +
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
      y = 290000,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = actin_UV, comparisons = list(c("bj p12.5 -uv act", "bj p12.5 +uv act")), annotations = "ns", map_signif_level = TRUE, y_position = 300000, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/actUVsum3d.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")


# boxplot for actin_UV mean_length_3d
summary_stats <- actin_UV %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p12.5 -uv act", "bj p12.5 +uv act")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_length_3d, na.rm = TRUE),
    sd = sd(mean_length_3d, na.rm = TRUE),
    ymax = max(mean_length_3d, na.rm = TRUE)
  )

ggplot(mutate(actin_UV, condition_clean = fct_relevel(condition_clean, "bj p12.5 -uv act", "bj p12.5 +uv act")), aes(condition_clean, mean_length_3d, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Condition", y = "Mean 3D filament length (μm)") +
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
      y = 8.55,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = actin_UV, comparisons = list(c("bj p12.5 -uv act", "bj p12.5 +uv act")), annotations = "ns", map_signif_level = TRUE, y_position = 8.65, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/actUVmean3d.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")


# boxplot for actin_UV mean_avg_angle
summary_stats <- actin_UV %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p12.5 -uv act", "bj p12.5 +uv act")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_avg_angle, na.rm = TRUE),
    sd = sd(mean_avg_angle, na.rm = TRUE),
    ymax = max(mean_avg_angle, na.rm = TRUE)
  )

ggplot(mutate(actin_UV, condition_clean = fct_relevel(condition_clean, "bj p12.5 -uv act", "bj p12.5 +uv act")), aes(condition_clean, mean_avg_angle, fill = condition_clean, color = condition_clean)) + 
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
      y = 87.7,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = actin_UV, comparisons = list(c("bj p12.5 -uv act", "bj p12.5 +uv act")), map_signif_level = TRUE, annotations = "ns", y_position = 88.2, color = 'black') 
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/actUVMeanAvgAng.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")


# boxplot for actin_UV mean_deviation
summary_stats <- actin_UV %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p12.5 -uv act", "bj p12.5 +uv act")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_deviation, na.rm = TRUE),
    sd = sd(mean_deviation, na.rm = TRUE),
    ymax = max(mean_deviation, na.rm = TRUE)
  )

ggplot(mutate(actin_UV, condition_clean = fct_relevel(condition_clean, "bj p12.5 -uv act", "bj p12.5 +uv act")), aes(condition_clean, mean_deviation, fill = condition_clean, color = condition_clean)) + 
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
      y = 1.42,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = actin_UV, comparisons = list(c("bj p12.5 -uv act", "bj p12.5 +uv act")), annotations = "ns", map_signif_level = TRUE, y_position = 1.44,  color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/actUVmeanDev.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")


# boxplot for actin_UV mean_acc_curv
summary_stats <- actin_UV %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p12.5 -uv act", "bj p12.5 +uv act")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_acc_curv, na.rm = TRUE),
    sd = sd(mean_acc_curv, na.rm = TRUE),
    ymax = max(mean_acc_curv, na.rm = TRUE)
  )

ggplot(mutate(actin_UV, condition_clean = fct_relevel(condition_clean, "bj p12.5 -uv act", "bj p12.5 +uv act")), aes(condition_clean, mean_acc_curv, fill = condition_clean, color = condition_clean)) + 
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
      y = 4.18,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = actin_UV, comparisons = list(c("bj p12.5 -uv act", "bj p12.5 +uv act")), annotations = "ns", map_signif_level = TRUE, y_position = 4.19, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/actUVmeanAccCurv.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")

# boxplot for actin_UV mean_net_curv
summary_stats <- actin_UV %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p12.5 -uv act", "bj p12.5 +uv act")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_net_curv, na.rm = TRUE),
    sd = sd(mean_net_curv, na.rm = TRUE),
    ymax = max(mean_net_curv, na.rm = TRUE)
  )

ggplot(mutate(actin_UV, condition_clean = fct_relevel(condition_clean, "bj p12.5 -uv act", "bj p12.5 +uv act")), aes(condition_clean, mean_net_curv, fill = condition_clean, color = condition_clean)) + 
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
      y = 1.55,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = actin_UV, comparisons = list(c("bj p12.5 -uv act", "bj p12.5 +uv act")), annotations = "ns", map_signif_level = TRUE, y_position = 1.57, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/actUVmeanNetCurv.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")

# boxplot for actin_UV filament_count
summary_stats <- actin_UV %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p12.5 -uv act", "bj p12.5 +uv act")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(filament_count, na.rm = TRUE),
    sd = sd(filament_count, na.rm = TRUE),
    ymax = max(filament_count, na.rm = TRUE)
  )

ggplot(mutate(actin_UV, condition_clean = fct_relevel(condition_clean, "bj p12.5 -uv act", "bj p12.5 +uv act")), aes(condition_clean, filament_count, fill = condition_clean, color = condition_clean)) + 
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
      y = 39000,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = actin_UV, comparisons = list(c("bj p12.5 -uv act", "bj p12.5 +uv act")), annotations = "ns", map_signif_level = TRUE, y_position = , color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/actUVFilCount.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")

# boxplot for actin_UV sum_branches
summary_stats <- actin_UV %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p12.5 -uv act", "bj p12.5 +uv act")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(sum_branches, na.rm = TRUE),
    sd = sd(sum_branches, na.rm = TRUE),
    ymax = max(sum_branches, na.rm = TRUE)
  )

ggplot(mutate(actin_UV, condition_clean = fct_relevel(condition_clean, "bj p12.5 -uv act", "bj p12.5 +uv act")), aes(condition_clean, sum_branches, fill = condition_clean, color = condition_clean)) + 
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
      y = 30500,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = actin_UV, comparisons = list(c("bj p12.5 -uv act", "bj p12.5 +uv act")), annotations = "ns", map_signif_level = TRUE, y_position = , color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/actUVSumBranch.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")

# boxplot for actin_UV mean_branches
summary_stats <- actin_UV %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p12.5 -uv act", "bj p12.5 +uv act")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_branches, na.rm = TRUE),
    sd = sd(mean_branches, na.rm = TRUE),
    ymax = max(mean_branches, na.rm = TRUE)
  )

ggplot(mutate(actin_UV, condition_clean = fct_relevel(condition_clean, "bj p12.5 -uv act", "bj p12.5 +uv act")), aes(condition_clean, mean_branches, fill = condition_clean, color = condition_clean)) + 
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
      y = 0.84,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = actin_UV, comparisons = list(c("bj p12.5 -uv act", "bj p12.5 +uv act")), annotations = "ns", map_signif_level = TRUE, y_position = 0.85, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/actUVMeanBranch.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")


# boxplot for actin_etop sum_length_2d
summary_stats <- actin_etop %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p8.5 act", "bj p9.5 etop vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(sum_length_2d, na.rm = TRUE),
    sd = sd(sum_length_2d, na.rm = TRUE),
    ymax = max(sum_length_2d, na.rm = TRUE)
  )

ggplot(mutate(actin_etop, condition_clean = fct_relevel(condition_clean, "bj p8.5 act", "bj p9.5 etop vim")), aes(condition_clean, sum_length_2d, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Condition", y = "Total 2D filament length per cell (μm)") +
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
      y = 78500,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = actin_etop, comparisons = list(c("bj p8.5 act","bj p9.5 etop vim")), annotations = "***", map_signif_level = TRUE, vjust = 0.7, y_position = 82000, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/actEtopsum2d.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")

# boxplot for actin_etop mean_length_2d
summary_stats <- actin_etop %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p8.5 act", "bj p9.5 etop vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_length_2d, na.rm = TRUE),
    sd = sd(mean_length_2d, na.rm = TRUE),
    ymax = max(mean_length_2d, na.rm = TRUE)
  )

ggplot(mutate(actin_etop, condition_clean = fct_relevel(condition_clean, "bj p8.5 act", "bj p9.5 etop vim")), aes(condition_clean, mean_length_2d, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Condition", y = "Mean 2D filament length per cell (μm)") +
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
      y = 1.14,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = actin_etop, comparisons = list(c("bj p8.5 act","bj p9.5 etop vim")), annotations = "**", map_signif_level = TRUE, vjust = 0.7, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/actEtopmean2d.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")

# boxplot for actin_etop sum_length_3d
summary_stats <- actin_etop %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p8.5 act", "bj p9.5 etop vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(sum_length_3d, na.rm = TRUE),
    sd = sd(sum_length_3d, na.rm = TRUE),
    ymax = max(sum_length_3d, na.rm = TRUE)
  )

ggplot(mutate(actin_etop, condition_clean = fct_relevel(condition_clean, "bj p8.5 act", "bj p9.5 etop vim")), aes(condition_clean, sum_length_3d, fill = condition_clean, color = condition_clean)) + 
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
      y = 420000,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = actin_etop, comparisons = list(c("bj p8.5 act","bj p9.5 etop vim")), annotations = "**", map_signif_level = TRUE, vjust = 0.7, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/actEtopSum3d.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")

# boxplot for actin_etop mean_length_3d
summary_stats <- actin_etop %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p8.5 act", "bj p9.5 etop vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_length_3d, na.rm = TRUE),
    sd = sd(mean_length_3d, na.rm = TRUE),
    ymax = max(mean_length_3d, na.rm = TRUE)
  )

ggplot(mutate(actin_etop, condition_clean = fct_relevel(condition_clean, "bj p8.5 act", "bj p9.5 etop vim")), aes(condition_clean, mean_length_3d, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Condition", y = "Mean 3D filament length per cell (μm)") +
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
      y = 9.4,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = actin_etop, comparisons = list(c("bj p8.5 act", "bj p9.5 etop vim")), annotations = "ns", map_signif_level = TRUE, y_position = 9.6, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/actEtopMean3d.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")


# boxplot for actin_etop mean_avg_angle
summary_stats <- actin_etop %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p8.5 act", "bj p9.5 etop vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_avg_angle, na.rm = TRUE),
    sd = sd(mean_avg_angle, na.rm = TRUE),
    ymax = max(mean_avg_angle, na.rm = TRUE)
  )

ggplot(mutate(actin_etop, condition_clean = fct_relevel(condition_clean, "bj p8.5 act", "bj p9.5 etop vim")), aes(condition_clean, mean_avg_angle, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Passage", y = "Mean average filament angle (deg)") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('p8.5 -etop', 'p9.5 +etop')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(1)) +
  geom_text(
    data = summary_stats,
    aes(
      x = condition_clean,
      y = 91.5,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = actin_etop, comparisons = list(c("bj p8.5 act", "bj p9.5 etop vim")), map_signif_level = TRUE, annotations = "***", vjust = 0.7, y_position = 92.5, color = 'black') 
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/actEtopMeanAvgAng.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")


# boxplot for actin_etop mean_deviation
summary_stats <- actin_etop %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p8.5 act", "bj p9.5 etop vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_deviation, na.rm = TRUE),
    sd = sd(mean_deviation, na.rm = TRUE),
    ymax = max(mean_deviation, na.rm = TRUE)
  )

ggplot(mutate(actin_etop, condition_clean = fct_relevel(condition_clean, "bj p8.5 act", "bj p9.5 etop vim")), aes(condition_clean, mean_deviation, fill = condition_clean, color = condition_clean)) + 
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
      y = 1.75,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = actin_etop, comparisons = list(c("bj p8.5 act", "bj p9.5 etop vim")), annotations = "ns", map_signif_level = TRUE, y_position = 1.8, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/actEtopmeanDev.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")


# boxplot for actin_etop mean_acc_curv
summary_stats <- actin_etop %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p8.5 act", "bj p9.5 etop vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_acc_curv, na.rm = TRUE),
    sd = sd(mean_acc_curv, na.rm = TRUE),
    ymax = max(mean_acc_curv, na.rm = TRUE)
  )

ggplot(mutate(actin_etop, condition_clean = fct_relevel(condition_clean, "bj p8.5 act", "bj p9.5 etop vim")), aes(condition_clean, mean_acc_curv, fill = condition_clean, color = condition_clean)) + 
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
      y = 4.25,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = actin_etop, comparisons = list(c("bj p8.5 act", "bj p9.5 etop vim")), annotations = "ns", map_signif_level = TRUE, y_position = 4.265, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/actEtopMeanAccCurv.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")

# boxplot for actin_etop mean_net_curv
summary_stats <- actin_etop %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p8.5 act", "bj p9.5 etop vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_net_curv, na.rm = TRUE),
    sd = sd(mean_net_curv, na.rm = TRUE),
    ymax = max(mean_net_curv, na.rm = TRUE)
  )

ggplot(mutate(actin_etop, condition_clean = fct_relevel(condition_clean, "bj p8.5 act", "bj p9.5 etop vim")), aes(condition_clean, mean_net_curv, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Condition", y = "Mean end-to-end filament curvature /cell") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('p8.5 -etop', 'p9.5 +etop')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(0.01)) +  geom_text(
    data = summary_stats,
    aes(
      x = condition_clean,
      y = 1.6,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = actin_etop, comparisons = list(c("bj p8.5 act", "bj p9.5 etop vim")), annotations = "**", map_signif_level = TRUE, vjust = 0.7, y_position = 1.62, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/actEtopMeanNetCurv.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")

# boxplot for actin_etop filament_count
summary_stats <- actin_etop %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p8.5 act", "bj p9.5 etop vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(filament_count, na.rm = TRUE),
    sd = sd(filament_count, na.rm = TRUE),
    ymax = max(filament_count, na.rm = TRUE)
  )

ggplot(mutate(actin_etop, condition_clean = fct_relevel(condition_clean, "bj p8.5 act", "bj p9.5 etop vim")), aes(condition_clean, filament_count, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Condition", y = "Number of filaments /cell") +
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
      y = 107000,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = actin_etop, comparisons = list(c("bj p8.5 act", "bj p9.5 etop vim")), annotations = "***", map_signif_level = TRUE, vjust = 0.7, y_position = 112000, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/actEtopFilCount.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")

# boxplot for actin_etop sum_branches
summary_stats <- actin_etop %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p8.5 act", "bj p9.5 etop vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(sum_branches, na.rm = TRUE),
    sd = sd(sum_branches, na.rm = TRUE),
    ymax = max(sum_branches, na.rm = TRUE)
  )

ggplot(mutate(actin_etop, condition_clean = fct_relevel(condition_clean, "bj p8.5 act", "bj p9.5 etop vim")), aes(condition_clean, sum_branches, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Condition", y = "Number of branches /cell") +
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
      y = 68000,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = actin_etop, comparisons = list(c("bj p8.5 act", "bj p9.5 etop vim")), annotations = "**", map_signif_level = TRUE, vjust = 0.7, y_position = 71000, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/actEtopSumBranch.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")

# boxplot for actin_etop mean_branches
summary_stats <- actin_etop %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p8.5 act", "bj p9.5 etop vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_branches, na.rm = TRUE),
    sd = sd(mean_branches, na.rm = TRUE),
    ymax = max(mean_branches, na.rm = TRUE)
  )

ggplot(mutate(actin_etop, condition_clean = fct_relevel(condition_clean, "bj p8.5 act", "bj p9.5 etop vim")), aes(condition_clean, mean_branches, fill = condition_clean, color = condition_clean)) + 
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
      y = 0.92,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = actin_etop, comparisons = list(c("bj p8.5 act", "bj p9.5 etop vim")), annotations = "***", map_signif_level = TRUE, vjust = 0.7, y_position = 0.95, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/actEtopMeanBranch.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")

# doplot for actin_etop sum_length_3d vs sum_branches
ggplot(data = actin_etop, aes(x = sum_length_3d, y = sum_branches, fill = condition_clean, color = condition_clean)) + 
  geom_point(size = 0.8) + 
  labs(x = "Total 3D filament length (μm)", y = "Total branches /cell") +
  theme_bw() +
  scale_color_manual(values = c("purple", "green")) +
  scale_fill_manual(values = c("purple", "green")) +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE), x = guide_axis(minor.ticks = TRUE)) +
  scale_y_continuous(minor_breaks = scales::breaks_width(5000)) +
  scale_x_continuous(minor_breaks = scales::breaks_width(100000)) 
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/vimPassSum2D.tiff", device = "tiff", dpi = 300, width = 177, height = 150, units = "mm")

# boxplot for actin_primary sum_length_2d
summary_stats <- actin_primary %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "primary fibro 18yo act", "primary fibro 29yo act", "primary fibro 32yo act", "primary fibro 59yo act")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(sum_length_2d, na.rm = TRUE),
    sd = sd(sum_length_2d, na.rm = TRUE),
    ymax = max(sum_length_2d, na.rm = TRUE)
  )

ggplot(mutate(actin_primary, condition_clean = fct_relevel(condition_clean, "primary fibro 18yo act", "primary fibro 29yo act", "primary fibro 32yo act", "primary fibro 59yo act")), aes(condition_clean, sum_length_2d, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Age", y = "Total 2D filament length per cell (μm)") +
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
      y = 65000,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = actin_primary, comparisons = list(c("primary fibro 18yo act", "primary fibro 59yo act")), annotations = "***", map_signif_level = TRUE, y_position = 79000, vjust = 0.7, color = 'black') +
  geom_signif(data = actin_primary, comparisons = list(c("primary fibro 29yo act", "primary fibro 59yo act")), annotations = "***", map_signif_level = TRUE, y_position = 76000, vjust = 0.7, color = 'black') +
  geom_signif(data = actin_primary, comparisons = list(c("primary fibro 32yo act", "primary fibro 59yo act")), annotations = "**", map_signif_level = TRUE, y_position = 73000, vjust = 0.7, color = 'black') 
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/actPrimSum2D.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")


# boxplot for actin_primary mean_length_2d
summary_stats <- actin_primary %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "primary fibro 18yo act", "primary fibro 29yo act", "primary fibro 32yo act", "primary fibro 59yo act")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_length_2d, na.rm = TRUE),
    sd = sd(mean_length_2d, na.rm = TRUE),
    ymax = max(mean_length_2d, na.rm = TRUE)
  )

ggplot(mutate(actin_primary, condition_clean = fct_relevel(condition_clean, "primary fibro 18yo act", "primary fibro 29yo act", "primary fibro 32yo act", "primary fibro 59yo act")), aes(condition_clean, mean_length_2d, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Age", y = "Mean 2D filament length per cell (μm)") +
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
      y = 1.58,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = actin_primary, comparisons = list(c("primary fibro 18yo act", "primary fibro 59yo act")), annotations = "ns", map_signif_level = TRUE, y_position = 1.63, color = 'black') 
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/actPrimMean2D.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")

# boxplot for actin_primary sum_length_3d
summary_stats <- actin_primary %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "primary fibro 18yo act", "primary fibro 29yo act", "primary fibro 32yo act", "primary fibro 59yo act")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(sum_length_3d, na.rm = TRUE),
    sd = sd(sum_length_3d, na.rm = TRUE),
    ymax = max(sum_length_3d, na.rm = TRUE)
  )

ggplot(mutate(actin_primary, condition_clean = fct_relevel(condition_clean, "primary fibro 18yo act", "primary fibro 29yo act", "primary fibro 32yo act", "primary fibro 59yo act")), aes(condition_clean, sum_length_3d, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Age", y = "Total 3D filament length per cell (μm)") +
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
      y = 480000,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = actin_primary, comparisons = list(c("primary fibro 18yo act", "primary fibro 59yo act")), annotations = "***", map_signif_level = TRUE, vjust = 0.7, y_position = 550000, color = 'black') +
  geom_signif(data = actin_primary, comparisons = list(c("primary fibro 29yo act", "primary fibro 59yo act")), annotations = "**", map_signif_level = TRUE, vjust = 0.7, y_position = 530000, color = 'black') +
  geom_signif(data = actin_primary, comparisons = list(c("primary fibro 32yo act", "primary fibro 59yo act")), annotations = "*", map_signif_level = TRUE, vjust = 0.7, y_position = 510000, color = 'black') 
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/actPrimSum3D.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")

# boxplot for actin_primary mean_length_3d
summary_stats <- actin_primary %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "primary fibro 18yo act", "primary fibro 29yo act", "primary fibro 32yo act", "primary fibro 59yo act")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_length_3d, na.rm = TRUE),
    sd = sd(mean_length_3d, na.rm = TRUE),
    ymax = max(mean_length_3d, na.rm = TRUE)
  )

ggplot(mutate(actin_primary, condition_clean = fct_relevel(condition_clean, "primary fibro 18yo act", "primary fibro 29yo act", "primary fibro 32yo act", "primary fibro 59yo act")), aes(condition_clean, mean_length_3d, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Age", y = "Mean 3D filament length (μm)") +
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
      y = 9.3,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = actin_primary, comparisons = list(c("primary fibro 18yo act", "primary fibro 59yo act")), annotations = "ns", map_signif_level = TRUE, y_position = 9.5, color = 'black') 
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/actPrimmean3d.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")

# boxplot for actin_primary mean_avg_angle
summary_stats <- actin_primary %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "primary fibro 18yo act", "primary fibro 29yo act", "primary fibro 32yo act", "primary fibro 59yo act")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_avg_angle, na.rm = TRUE),
    sd = sd(mean_avg_angle, na.rm = TRUE),
    ymax = max(mean_avg_angle, na.rm = TRUE)
  )

ggplot(mutate(actin_primary, condition_clean = fct_relevel(condition_clean, "primary fibro 18yo act", "primary fibro 29yo act", "primary fibro 32yo act", "primary fibro 59yo act")), aes(condition_clean, mean_avg_angle, fill = condition_clean, color = condition_clean)) + 
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
      y = 92,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = actin_primary, comparisons = list(c("primary fibro 18yo act", "primary fibro 59yo act")), annotations = "ns", map_signif_level = TRUE, y_position = 93.5, color = 'black') 
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/actPrimmeanAvgAng.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")


# boxplot for actin_primary mean_deviation
summary_stats <- actin_primary %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "primary fibro 18yo act", "primary fibro 29yo act", "primary fibro 32yo act", "primary fibro 59yo act")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_deviation, na.rm = TRUE),
    sd = sd(mean_deviation, na.rm = TRUE),
    ymax = max(mean_deviation, na.rm = TRUE)
  )

ggplot(mutate(actin_primary, condition_clean = fct_relevel(condition_clean, "primary fibro 18yo act", "primary fibro 29yo act", "primary fibro 32yo act", "primary fibro 59yo act")), aes(condition_clean, mean_deviation, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Age", y = "Mean filament deviation (μm)") +
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
      y = 1.63,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = actin_primary, comparisons = list(c("primary fibro 18yo act", "primary fibro 59yo act")), annotations = "ns", map_signif_level = TRUE, y_position = 1.67, color = 'black') 
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/actPrimmeanDev.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")

# boxplot for actin_primary mean_acc_curv
summary_stats <- actin_primary %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "primary fibro 18yo act", "primary fibro 29yo act", "primary fibro 32yo act", "primary fibro 59yo act")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_acc_curv, na.rm = TRUE),
    sd = sd(mean_acc_curv, na.rm = TRUE),
    ymax = max(mean_acc_curv, na.rm = TRUE)
  )

ggplot(mutate(actin_primary, condition_clean = fct_relevel(condition_clean, "primary fibro 18yo act", "primary fibro 29yo act", "primary fibro 32yo act", "primary fibro 59yo act")), aes(condition_clean, mean_acc_curv, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Age", y = "Mean accumulated filament curvature /cell") +
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
      y = 4.24,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = actin_primary, comparisons = list(c("primary fibro 18yo act", "primary fibro 59yo act")), annotations = "ns", map_signif_level = TRUE, y_position = 4.27, color = 'black') 
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/actPrimmeanAccCurv.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")


# boxplot for actin_primary mean_net_curv
summary_stats <- actin_primary %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "primary fibro 18yo act", "primary fibro 29yo act", "primary fibro 32yo act", "primary fibro 59yo act")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_net_curv, na.rm = TRUE),
    sd = sd(mean_net_curv, na.rm = TRUE),
    ymax = max(mean_net_curv, na.rm = TRUE)
  )

ggplot(mutate(actin_primary, condition_clean = fct_relevel(condition_clean, "primary fibro 18yo act", "primary fibro 29yo act", "primary fibro 32yo act", "primary fibro 59yo act")), aes(condition_clean, mean_net_curv, fill = condition_clean, color = condition_clean)) + 
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
      y = 1.58,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = actin_primary, comparisons = list(c("primary fibro 18yo act", "primary fibro 59yo act")), annotations = "ns", map_signif_level = TRUE, y_position = 1.6, color = 'black') 
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/actPrimmeanNetCurv.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")

# boxplot for actin_primary filament_count
summary_stats <- actin_primary %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "primary fibro 18yo act", "primary fibro 29yo act", "primary fibro 32yo act", "primary fibro 59yo act")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(filament_count, na.rm = TRUE),
    sd = sd(filament_count, na.rm = TRUE),
    ymax = max(filament_count, na.rm = TRUE)
  )

ggplot(mutate(actin_primary, condition_clean = fct_relevel(condition_clean, "primary fibro 18yo act", "primary fibro 29yo act", "primary fibro 32yo act", "primary fibro 59yo act")), aes(condition_clean, filament_count, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Age", y = "Filament count /cell") +
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
      y = 66000,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = actin_primary, comparisons = list(c("primary fibro 18yo act", "primary fibro 59yo act")), annotations = "***", map_signif_level = TRUE, vjust = 0.7, y_position = 75000, color = 'black') +
  geom_signif(data = actin_primary, comparisons = list(c("primary fibro 29yo act", "primary fibro 59yo act")), annotations = "**", map_signif_level = TRUE, vjust = 0.7, y_position = 72000, color = 'black') +
  geom_signif(data = actin_primary, comparisons = list(c("primary fibro 32yo act", "primary fibro 59yo act")), annotations = "**", map_signif_level = TRUE, vjust = 0.7, y_position = 69000, color = 'black') 
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/actPrimFilCount.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")

# boxplot for actin_primary sum_branches
summary_stats <- actin_primary %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "primary fibro 18yo act", "primary fibro 29yo act", "primary fibro 32yo act", "primary fibro 59yo act")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(sum_branches, na.rm = TRUE),
    sd = sd(sum_branches, na.rm = TRUE),
    ymax = max(sum_branches, na.rm = TRUE)
  )

ggplot(mutate(actin_primary, condition_clean = fct_relevel(condition_clean, "primary fibro 18yo act", "primary fibro 29yo act", "primary fibro 32yo act", "primary fibro 59yo act")), aes(condition_clean, sum_branches, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Age", y = "Total filament branches /cell") +
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
      y = 54000,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = actin_primary, comparisons = list(c("primary fibro 18yo act", "primary fibro 59yo act")), annotations = "***", map_signif_level = TRUE, vjust = 0.7, y_position = 63000, color = 'black') +
  geom_signif(data = actin_primary, comparisons = list(c("primary fibro 29yo act", "primary fibro 59yo act")), annotations = "**", map_signif_level = TRUE, vjust = 0.7, y_position = 60000, color = 'black') +
  geom_signif(data = actin_primary, comparisons = list(c("primary fibro 32yo act", "primary fibro 59yo act")), annotations = "**", map_signif_level = TRUE, vjust = 0.7, y_position = 57000, color = 'black') 
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/actPrimBranchCount.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")

# boxplot for actin_primary mean_branches
summary_stats <- actin_primary %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "primary fibro 18yo act", "primary fibro 29yo act", "primary fibro 32yo act", "primary fibro 59yo act")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_branches, na.rm = TRUE),
    sd = sd(mean_branches, na.rm = TRUE),
    ymax = max(mean_branches, na.rm = TRUE)
  )

ggplot(mutate(actin_primary, condition_clean = fct_relevel(condition_clean, "primary fibro 18yo act", "primary fibro 29yo act", "primary fibro 32yo act", "primary fibro 59yo act")), aes(condition_clean, mean_branches, fill = condition_clean, color = condition_clean)) + 
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
      y = 0.87,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = actin_primary, comparisons = list(c("primary fibro 32yo act", "primary fibro 59yo act")), annotations = "**", map_signif_level = TRUE, vjust = 0.7, y_position = 0.88, color = 'black') 
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/actPrimMeanBranch.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")

# doplot for actin_primary sum_length_3d vs sum_branches
ggplot(data = actin_primary, aes(x = sum_length_3d, y = sum_branches, fill = condition_clean, color = condition_clean)) + 
  geom_point(size = 0.8) + 
  labs(x = "Total 3D filament length (μm)", y = "Total branches /cell") +
  theme_bw() +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, y = guide_axis(minor.ticks = TRUE), x = guide_axis(minor.ticks = TRUE)) +
  scale_y_continuous(minor_breaks = scales::breaks_width(5000)) 
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/vimPassSum2D.tiff", device = "tiff", dpi = 300, width = 177, height = 150, units = "mm")
