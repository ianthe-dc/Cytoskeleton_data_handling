library(tidyverse)
library(lme4)
library(car)
library(emmeans)
library(ggprism)
library(forcats)
library(dplyr)
library(ggsignif)
library(DHARMa)

cell_actin_fixed_widths <- read.csv("fixed_cell_actin_widths.csv")
actin_BJ <- filter(cell_actin_fixed_widths, condition_clean %in% c("bj p8.5 act", "bj p9.5 etop vim", "bj p11.5 act", "bj p12.5 -uv act", "bj p12.5 +uv act", "bj p26.5 act", "bj p31.5 act", "bj p36.5 act", "bj p38.5 act", "bj p40.5 act"))
actin_primary <- cell_actin_fixed_widths[258:309,]
actin_UV <- filter(cell_actin_fixed_widths, condition_clean %in% c("bj p12.5 -uv act", "bj p12.5 +uv act"))
actin_etop <- filter(cell_actin_fixed_widths, condition_clean %in% c("bj p8.5 act", "bj p9.5 etop vim"))
actin_pass <- filter(cell_actin_fixed_widths, condition_clean %in% c("bj p11.5 act", "bj p26.5 act", "bj p31.5 act", "bj p36.5 act", "bj p38.5 act", "bj p40.5 act"))

#stats for all actin BJ datasets
model <- lmer(mean_filament_width ~ condition_clean + (1|condition_clean:biorep) +(1|condition_clean:biorep:techrep), actin_BJ)
summary(model)
Anova(model)
em1 <- emmeans(model, ~condition_clean)
pairs(em1)


#stats for passage only 
model_pass <- lmer(mean_filament_width ~ condition_clean + (1|condition_clean:biorep:techrep), actin_pass)
summary(model_pass)
Anova(model_pass)
em_pass <- emmeans(model_pass, ~condition_clean)
pairs(em_pass)

simulationOutput <- simulateResiduals(fittedModel = model_pass, plot = T)

#stats for actin primary
model13 <- lm(mean_filament_width ~ condition_clean + techrep, actin_primary)
summary(model13)
Anova(model13)
em13 <- emmeans(model13, ~condition_clean)
pairs(em13)

simulationOutput <- simulateResiduals(fittedModel = model13, plot = T)

#stats for UV only
model_UV <- lm(mean_filament_width ~ condition_clean + techrep, actin_UV)
summary(model_UV)
Anova(model_UV) 
em_UV <- emmeans(model_UV, ~condition_clean)
pairs(em_UV)

simulationOutput <- simulateResiduals(fittedModel = model_UV, plot = T)

#stats for etop only
model_etop <- lmer(mean_filament_width ~ condition_clean + (1|condition_clean:biorep:techrep), actin_etop)
summary(model_etop)
Anova(model_etop)
em_etop <- emmeans(model_etop, ~condition_clean)
pairs(em_etop)

kruskal <- kruskal.test(mean_filament_width ~ condition_clean, data = actin_etop)
kruskal

simulationOutput <- simulateResiduals(fittedModel = model_etop, plot = T)

# boxplot for all actin 
ggplot(mutate(actin_BJ, condition_clean = fct_relevel(condition_clean, "bj p8.5 act", "bj p9.5 etop vim", "bj p11.5 act", "bj p12.5 -uv act", "bj p12.5 +uv act", "bj p26.5 act", "bj p31.5 act", "bj p36.5 act", "bj p38.5 act", "bj p40.5 act")), aes(condition_clean, mean_filament_width, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.5) +
  labs(x = "Condition", y = "Total 2D filament length per cell (μm)") +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('p8.5', 'p9.5 +etop', 'p11.5', 'p12.5 -UV', 'p12.5 +UV', 'p26.5', 'p31.5', 'p36.5', 'p38.5', 'p40.5')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(5000))
ggsave("vimBJsum2D.tiff", device = "tiff", dpi = 300, width = 177, height = 150, units = "mm")

ggplot(actin_etop, aes(x=mean_filament_width, color = condition_clean)) + geom_density()

# boxplot for actin_pass 
summary_stats <- actin_pass %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p11.5 act","bj p26.5 act","bj p31.5 act",
                                       "bj p36.5 act","bj p38.5 act","bj p40.5 act")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_filament_width, na.rm = TRUE),
    sd = sd(mean_filament_width, na.rm = TRUE),
    ymax = max(mean_filament_width, na.rm = TRUE)
  )

ggplot(mutate(actin_pass, condition_clean = fct_relevel(condition_clean, "bj p11.5 act", "bj p26.5 act", "bj p31.5 act", "bj p36.5 act", "bj p38.5 act", "bj p40.5 act")), aes(condition_clean, mean_filament_width, fill = condition_clean, color = condition_clean)) + 
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
      y = 0.283,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = actin_pass, comparisons = list(c("bj p11.5 act","bj p40.5 act")), map_signif_level = TRUE, annotations = "*", vjust = 0.7, y_position = 0.29, color = 'black') +
  geom_signif(data = actin_pass, comparisons = list(c("bj p11.5 act","bj p36.5 act")), annotations = "**", map_signif_level = TRUE, vjust = 0.7, y_position = 0.286, color = 'black') 
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/actPasswidth.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")


# boxplot for actin_UV 
summary_stats <- actin_UV %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p12.5 -uv act", "bj p12.5 +uv act")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_filament_width, na.rm = TRUE),
    sd = sd(mean_filament_width, na.rm = TRUE),
    ymax = max(mean_filament_width, na.rm = TRUE)
  )

ggplot(mutate(actin_UV, condition_clean = fct_relevel(condition_clean, "bj p12.5 -uv act", "bj p12.5 +uv act")), aes(condition_clean, mean_filament_width, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Condition", y = "Cell average filament width (μm)") +
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
      y = 0.268,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = actin_UV, comparisons = list(c("bj p12.5 -uv act", "bj p12.5 +uv act")), annotations = "ns", map_signif_level = TRUE, y_position = 0.27, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/actUVwidth.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")


# boxplot for actin_etop 
summary_stats <- actin_etop %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p8.5 act", "bj p9.5 etop vim")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_filament_width, na.rm = TRUE),
    sd = sd(mean_filament_width, na.rm = TRUE),
    ymax = max(mean_filament_width, na.rm = TRUE)
  )

ggplot(mutate(actin_etop, condition_clean = fct_relevel(condition_clean, "bj p8.5 act", "bj p9.5 etop vim")), aes(condition_clean, mean_filament_width, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Condition", y = "Cell average filament width (μm)") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('p8.5 -etop', 'p9.5 +etop')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(0.005)) +
  geom_text(
    data = summary_stats,
    aes(
      x = condition_clean,
      y = 0.287,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data =actin_etop, comparisons = list(c("bj p8.5 act","bj p9.5 etop vim")), annotations = "***", vjust = 0.7, y_position = 0.29, map_signif_level = TRUE, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/actEtopwidth.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")


# boxplot for actin_primary 
summary_stats <- actin_primary %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "primary fibro 18yo act", "primary fibro 29yo act", "primary fibro 32yo act", "primary fibro 59yo act")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_filament_width, na.rm = TRUE),
    sd = sd(mean_filament_width, na.rm = TRUE),
    ymax = max(mean_filament_width, na.rm = TRUE)
  )

ggplot(mutate(actin_primary, condition_clean = fct_relevel(condition_clean, "primary fibro 18yo act", "primary fibro 29yo act", "primary fibro 32yo act", "primary fibro 59yo act")), aes(condition_clean, mean_filament_width, fill = condition_clean, color = condition_clean)) + 
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
      y = 0.255,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data =actin_primary, comparisons = list(c("primary fibro 18yo act","primary fibro 59yo act")), annotations = "ns", map_signif_level = TRUE, y_position = 0.26, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/actPrimwidth.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")
