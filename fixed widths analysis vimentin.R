library(tidyverse)
library(lme4)
library(car)
library(emmeans)
library(ggprism)
library(forcats)
library(dplyr)
library(ggsignif)
library(DHARMa)

cell_vimentin_fixed_widths <- read.csv("fixed_cell_vimentin_widths.csv")
vimentin_BJ <- filter(cell_vimentin_fixed_widths, condition_clean %in% c("bj p8.5 vim", "bj p9.5 etop vim", "bj p11.5 vim", "bj p12.5 -uv vim", "bj p12.5 +uv vim", "bj p26.5 vim", "bj p31.5 vim", "bj p36.5 vim", "bj p38.5 vim", "bj p40.5 vim"))
vimentin_primary <- cell_vimentin_fixed_widths[271:329,]
vimentin_UV <- filter(cell_vimentin_fixed_widths, condition_clean %in% c("bj p12.5 -uv vim", "bj p12.5 +uv vim"))
vimentin_etop <- filter(cell_vimentin_fixed_widths, condition_clean %in% c("bj p8.5 vim", "bj p9.5 etop vim"))
vimentin_pass <- filter(cell_vimentin_fixed_widths, condition_clean %in% c("bj p11.5 vim", "bj p26.5 vim", "bj p31.5 vim", "bj p36.5 vim", "bj p38.5 vim", "bj p40.5 vim"))


#stats for all actin BJ datasets
model <- lmer(mean_filament_width ~ condition_clean + (1|condition_clean:biorep) +(1|condition_clean:biorep:techrep), vimentin_BJ)
summary(model)
Anova(model)
em1 <- emmeans(model, ~condition_clean)
pairs(em1)


#stats for passage only 
model_pass <- lmer(mean_filament_width ~ condition_clean + (1|condition_clean:biorep:techrep), vimentin_pass)
summary(model_pass)
Anova(model_pass)
em_pass <- emmeans(model_pass, ~condition_clean)
pairs(em_pass)

simulationOutput <- simulateResiduals(fittedModel = model_pass, plot = T)

#stats for vimentin primary
model13 <- lmer(mean_filament_width ~ condition_clean + (1|condition_clean:biorep:techrep), vimentin_primary)
summary(model13)
Anova(model13)
em13 <- emmeans(model13, ~condition_clean)
pairs(em13)

simulationOutput <- simulateResiduals(fittedModel = model13, plot = T)

#stats for UV only
model_UV <- lmer(mean_filament_width ~ condition_clean + (1|condition_clean:biorep:techrep), vimentin_UV)
summary(model_UV)
Anova(model_UV)
em_UV <- emmeans(model_UV, ~condition_clean)
pairs(em_UV)

simulationOutput <- simulateResiduals(fittedModel = model_UV, plot = T)

#stats for etop only
model_etop <- lmer(mean_filament_width ~ condition_clean + (1|condition_clean:biorep) + (1|condition_clean:biorep:techrep), vimentin_etop)
summary(model_etop)
Anova(model_etop)
em_etop <- emmeans(model_etop, ~condition_clean)
pairs(em_etop)

simulationOutput <- simulateResiduals(fittedModel = model_etop, plot = T)

ggplot(vimentin_primary, aes(x=mean_filament_width)) + geom_density()

# boxplot for all vimentin 
ggplot(mutate(vimentin_BJ, condition_clean = fct_relevel(condition_clean, "bj p8.5 vim", "bj p9.5 etop vim", "bj p11.5 vim", "bj p12.5 -uv vim", "bj p12.5 +uv vim", "bj p26.5 vim", "bj p31.5 vim", "bj p36.5 vim", "bj p38.5 vim", "bj p40.5 vim")), aes(condition_clean, mean_filament_width, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Condition", y = "Total 2D filament length per cell (μm)") +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('p8.5', 'p9.5 +etop', 'p11.5', 'p12.5 -UV', 'p12.5 +UV', 'p26.5', 'p31.5', 'p36.5', 'p38.5', 'p40.5')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(5000))
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/vimBJsum2D.tiff", device = "tiff", dpi = 250, width = 177, height = 150, units = "mm")


# boxplot for vimentin_pass 
summary_stats <- vimentin_pass %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p11.5 vim","bj p26.5 vim","bj p31.5 vim",
                                       "bj p36.5 vim","bj p38.5 vim","bj p40.5 vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_filament_width, na.rm = TRUE),
    sd = sd(mean_filament_width, na.rm = TRUE),
    ymax = max(mean_filament_width, na.rm = TRUE)
  )

ggplot(mutate(vimentin_pass, condition_clean = fct_relevel(condition_clean, "bj p11.5 vim", "bj p26.5 vim", "bj p31.5 vim", "bj p36.5 vim", "bj p38.5 vim", "bj p40.5 vim")), aes(condition_clean, mean_filament_width, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Passage", y = "Cell average filament width (μm)") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('p11.5', 'p26.5', 'p31.5', 'p36.5', 'p38.5', 'p40.5')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(0.005)) +
  geom_text(
    data = summary_stats,
    aes(
      x = condition_clean,
      y = 0.27,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = vimentin_pass, comparisons = list(c("bj p11.5 vim","bj p40.5 vim")), map_signif_level = TRUE, annotations = "**", y_position = 0.273, vjust = 0.7, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/vimPasswidth.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")

# boxplot for vimentin_UV 
summary_stats <- vimentin_UV %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p12.5 -uv vim", "bj p12.5 +uv vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_filament_width, na.rm = TRUE),
    sd = sd(mean_filament_width, na.rm = TRUE),
    ymax = max(mean_filament_width, na.rm = TRUE)
  )

ggplot(mutate(vimentin_UV, condition_clean = fct_relevel(condition_clean, "bj p12.5 -uv vim", "bj p12.5 +uv vim")), aes(condition_clean, mean_filament_width, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Condition", y = "Cell average filament width (μm)") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('p12.5 -UV', 'p12.5 +UV')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(0.001)) +
  geom_text(
    data = summary_stats,
    aes(
      x = condition_clean,
      y = 0.244,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = vimentin_UV, comparisons = list(c("bj p12.5 -uv vim","bj p12.5 +uv vim")), annotations = "ns", map_signif_level = TRUE, y_position = 0.245, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/vimUVwidth.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")

# boxplot for vimentin_etop
summary_stats <- vimentin_etop %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p8.5 vim", "bj p9.5 etop vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_filament_width, na.rm = TRUE),
    sd = sd(mean_filament_width, na.rm = TRUE),
    ymax = max(mean_filament_width, na.rm = TRUE)
  )

ggplot(mutate(vimentin_etop, condition_clean = fct_relevel(condition_clean, "bj p8.5 vim", "bj p9.5 etop vim")), aes(condition_clean, mean_filament_width, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Condition", y = "Cell average filament width (μm)") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('p8.5 -etop', 'p9.5 +etop')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(0.002)) +
  geom_text(
    data = summary_stats,
    aes(
      x = condition_clean,
      y = 0.26,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = vimentin_etop, comparisons = list(c("bj p8.5 vim","bj p9.5 etop vim")), annotations = "ns", map_signif_level = TRUE, y_position = 0.262, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/vimEtopwidth.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")


# boxplot for vimentin_primary 
summary_stats <- vimentin_primary %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "primary fibro 18yo vim", "primary fibro 29yo vim", "primary fibro 32yo vim", "primary fibro 59yo vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_filament_width, na.rm = TRUE),
    sd = sd(mean_filament_width, na.rm = TRUE),
    ymax = max(mean_filament_width, na.rm = TRUE)
  )

ggplot(mutate(vimentin_primary, condition_clean = fct_relevel(condition_clean, "primary fibro 18yo vim", "primary fibro 29yo vim", "primary fibro 32yo vim", "primary fibro 59yo vim")), aes(condition_clean, mean_filament_width, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) +
  labs(x = "Age", y = "Cell average filament width (μm)") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('18 y/o', '29 y/o', '32 y/o', '59 y/o')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(0.002)) +
  geom_text(
    data = summary_stats,
    aes(
      x = condition_clean,
      y = 0.267,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = vimentin_primary, comparisons = list(c("primary fibro 18yo vim","primary fibro 59yo vim")), annotations = "ns", map_signif_level = TRUE, y_position = 0.269, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/vimPrimwidth.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")
