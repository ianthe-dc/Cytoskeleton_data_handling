library(tidyverse)
library(lme4)
library(car)
library(emmeans)
library(ggprism)
library(forcats)
library(dplyr)
library(ggsignif)
library(lmerTest)
library(DHARMa)

tubulin_df <- read.csv("fixed_cell_tubulin_widths.csv")
tubulin_BJ <- filter(tubulin_df, condition_clean %in% c("bj p11.5 tub", "bj p12.5 -uv tub", "bj p12.5 +uv tub", "bj p26.5 tub", "bj p31.5 tub", "bj p36.5 tub", "bj p38.5 tub", "bj p40.5 tub"))
tubulin_primary <- tubulin_df[167:201,]
tubulin_UV <- filter(tubulin_df, condition_clean %in% c("bj p12.5 -uv tub", "bj p12.5 +uv tub"))
tubulin_pass <- filter(tubulin_df, condition_clean %in% c("bj p11.5 tub", "bj p26.5 tub", "bj p31.5 tub", "bj p36.5 tub", "bj p38.5 tub", "bj p40.5 tub"))

#stats for all tubulin BJ datasets
model <- lmer(mean_filament_width ~ condition_clean + (1|condition_clean:biorep:techrep), tubulin_BJ)
summary(model)
Anova(model)
em1 <- emmeans(model, ~condition_clean)
pairs(em1)


kruskal <- kruskal.test(mean_filament_width ~ condition_clean, data = tubulin_pass)
summary(kruskal)
n <- pairwise.wilcox.test(tubulin_pass$mean_filament_width, tubulin_pass$condition_clean,
                     p.adjust.method = "BH")
n
#stats for passage only
model_pass <- lmer(mean_filament_width ~ condition_clean + (1|condition_clean:biorep:techrep), tubulin_pass)
summary(model_pass)
Anova(model_pass)
em_pass <- emmeans(model_pass, ~condition_clean)
pairs(em_pass)

simulationOutput <- simulateResiduals(fittedModel = model_pass, plot = T)

#stats for UV only
model_UV <- lmer(mean_filament_width ~ condition_clean + (1|condition_clean:biorep:techrep), tubulin_UV)
summary(model_UV)
Anova(model_UV)
em_UV <- emmeans(model_UV, ~condition_clean)
pairs(em_UV)

simulationOutput <- simulateResiduals(fittedModel = model_pass, plot = T)

#stats for tubulin_primary
model10_pass <- lm(mean_filament_width ~ condition_clean, tubulin_primary)
summary(model10_pass)
Anova(model10_pass)
em_pass10 <- emmeans(model10_pass, ~condition_clean)
pairs(em_pass10)

simulationOutput <- simulateResiduals(fittedModel = model_pass, plot = T)

welch_result <- oneway.test(mean_filament_width ~ condition_clean, data = tubulin_primary, var.equal = FALSE)
welch_result
pairwise.t.test(tubulin_primary$mean_filament_width, tubulin_primary$condition_clean, p.adj = "bonf", pool.sd = FALSE)

table(tubulin_primary$condition_clean, tubulin_primary$biorep)

ggplot(tubulin_pass, aes(x=mean_filament_width)) + geom_density()

# boxplot for all tubulin 
ggplot(mutate(tubulin_BJ, condition_clean = fct_relevel(condition_clean, "bj p11.5 tub", "bj p12.5 -uv tub", "bj p12.5 +uv tub", "bj p26.5 tub", "bj p31.5 tub", "bj p36.5 tub", "bj p38.5 tub", "bj p40.5 tub")), aes(condition_clean, mean_filament_width, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Condition", y = "Total 2D filament length per cell (μm)") +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('p11.5', 'p12.5 -UV', 'p12.5 +UV', 'p26.5', 'p31.5', 'p36.5', 'p38.5', 'p40.5')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(5000))
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/tubBJsum2D.tiff", device = "tiff", dpi = 300, width = 177, height = 150, units = "mm")


# boxplot for tubulin_pass 
summary_stats <- tubulin_pass %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p11.5 tub","bj p26.5 tub","bj p31.5 tub",
                                       "bj p36.5 tub","bj p38.5 tub","bj p40.5 tub")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_filament_width, na.rm = TRUE),
    sd = sd(mean_filament_width, na.rm = TRUE),
    ymax = max(mean_filament_width, na.rm = TRUE)
  )

ggplot(mutate(tubulin_pass, condition_clean = fct_relevel(condition_clean, "bj p11.5 tub", "bj p26.5 tub", "bj p31.5 tub", "bj p36.5 tub", "bj p38.5 tub", "bj p40.5 tub")), aes(condition_clean, mean_filament_width, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(x = "Passage", y = "Cell average filament width (μm)") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('p11.5', 'p26.5', 'p31.5', 'p36.5', 'p38.5', 'p40.5')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(0.002)) +
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
  geom_signif(data = tubulin_pass, comparisons = list(c("bj p31.5 tub","bj p36.5 tub")), map_signif_level = TRUE, annotations = "*", y_position = 0.247, vjust = 0.7, color = 'black') +
  geom_signif(data = tubulin_pass, comparisons = list(c("bj p36.5 tub","bj p38.5 tub")), map_signif_level = TRUE, annotations = "*", y_position = 0.246, vjust = 0.7, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/tubPasswidth.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")

# boxplot for tubulin_UV 
summary_stats <- tubulin_UV %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "bj p12.5 -uv tub", "bj p12.5 +uv tub")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_filament_width, na.rm = TRUE),
    sd = sd(mean_filament_width, na.rm = TRUE),
    ymax = max(mean_filament_width, na.rm = TRUE)
  )

ggplot(mutate(tubulin_UV, condition_clean = fct_relevel(condition_clean, "bj p12.5 -uv tub", "bj p12.5 +uv tub")), aes(condition_clean, mean_filament_width, fill = condition_clean, color = condition_clean)) + 
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
      y = 0.2235,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = tubulin_UV, comparisons = list(c("bj p12.5 -uv tub","bj p12.5 +uv tub")), annotations = "ns", map_signif_level = TRUE, y_position = 0.2239, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/tubUVwidth.tiff", device = "tiff", dpi = 300, width = 105, height = 120, units = "mm")

# boxplot for tubulin_primary 
summary_stats <- tubulin_primary %>%
  mutate(condition_clean = fct_relevel(condition_clean,
                                       "primary fibro 18yo tub", "primary fibro 29yo tub", "primary fibro 32yo tub", "primary fibro 59yo tub")) %>%
  group_by(condition_clean) %>%
  summarise(
    mean = mean(mean_filament_width, na.rm = TRUE),
    sd = sd(mean_filament_width, na.rm = TRUE),
    ymax = max(mean_filament_width, na.rm = TRUE)
  )

ggplot(mutate(tubulin_primary, condition_clean = fct_relevel(condition_clean, "primary fibro 18yo tub", "primary fibro 29yo tub", "primary fibro 32yo tub", "primary fibro 59yo tub")), aes(condition_clean, mean_filament_width, fill = condition_clean, color = condition_clean)) + 
  geom_boxplot() +
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
      y = 0.236,
      label = paste0("Mean=", round(mean, 2), "\nSD=", round(sd, 2))
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_signif(data = tubulin_primary, comparisons = list(c("primary fibro 18yo tub","primary fibro 59yo tub")), annotations = "**", map_signif_level = TRUE, y_position = 0.239, vjust = 0.7, color = 'black') +
  geom_signif(data = tubulin_primary, comparisons = list(c("primary fibro 18yo tub","primary fibro 32yo tub")), annotations = "*", map_signif_level = TRUE, y_position = 0.237, vjust = 0.7, color = 'black')
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/tubPrimwidth.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")
