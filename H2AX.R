library(tidyverse)
library(lme4)
library(car)
library(emmeans)
library(ggprism)
library(forcats)
library(dplyr)
library(ggsignif)
library(DHARMa)

h2ax_bj <- read.csv("all passages new.csv")
long_h2ax_bj <- stack(h2ax_bj)


# summary stats for h2ax staining intensity in BJ fibroblasts: UV, etoposide, and serially passaged 

summary_long_bj <- long_h2ax_bj %>%
  mutate(ind = fct_relevel(ind, "p11.5", "p8.5", "p9.5.etoposide", "p12.5.UV", "p12.5.UV.1", "p16.5", "p21.5", "p26.5", "p31.5", "p36.5", "p38.5", "p40.5")) %>%
  group_by(ind) %>%
  summarise(
    mean = mean(values, na.rm = TRUE),
    sd = sd(values, na.rm = TRUE)
  )

model <- lm(values ~ ind, long_h2ax_bj)
summary(model)
Anova(model)
em1 <- emmeans(model, ~ ind)
pairs(em1)

simulationOutput <- simulateResiduals(fittedModel = model, plot = T)

kruskal <- kruskal.test(values ~ ind, data = long_h2ax_bj)
kruskal
pairwise.wilcox.test(long_h2ax_bj$values, long_h2ax_bj$ind, p.adjust.method = "BH")


ggplot(mutate(long_h2ax_bj, ind = fct_relevel(ind, "p8.5", "p9.5.etoposide", "p11.5", "p12.5.UV", "p12.5.UV.1", "p16.5", "p21.5", "p26.5", "p31.5", "p36.5", "p38.5", "p40.5")), aes(ind, values, fill = ind, color = ind)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(y = "Average H2AX staining intensity") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('p8.5', 'p9.5 + etoposide', 'p11.5', 'p12.5 - UV', 'p12.5 +UV', 'p16.5', 'p21.5', 'p26.5', 'p31.5', 'p36.5', 'p38.5', 'p40.5')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(50)) 

ggsave("vimPassSum2D.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")


h2ax_prim_wf <- read.csv("Data 9.csv")
long_h2ax_prim_wf <- stack(h2ax_prim_wf)

# summary stats for h2ax staining intensity in primary fibroblasts 18yo vs 59yo widefield expt

summary_long_prim_wf <- long_h2ax_prim_wf %>%
  mutate(ind = fct_relevel(ind, "X18.yo", "X59.yo")) %>%
  group_by(ind) %>%
  summarise(
    mean = mean(values, na.rm = TRUE),
    sd = sd(values, na.rm = TRUE)
  )

model <- lm(values ~ ind, long_h2ax_prim_wf)
summary(model)
Anova(model)
em1 <- emmeans(model, ~ ind)
pairs(em1)

simulationOutput <- simulateResiduals(fittedModel = model, plot = T)

ggplot(mutate(long_h2ax_prim_wf, ind = fct_relevel(ind, "X18.yo", "X59.yo")), aes(ind, values, fill = ind, color = ind)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(y = "Average H2AX staining intensity") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('18 y/o', '59 y/o')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(1), limits = c(0, 25)) 


# summary stats for h2ax staining intensity in primary fibroblasts all

h2ax_prim <- read.csv("Data 8.csv")
long_h2ax_prim <- stack(h2ax_prim)

summary_long_prim <- long_h2ax_prim %>%
  mutate(ind = fct_relevel(ind, "X18.yo", "X29.yo", "X32.yo", "X59.yo")) %>%
  group_by(ind) %>%
  summarise(
    mean = mean(values, na.rm = TRUE),
    sd = sd(values, na.rm = TRUE)
  )

model <- lm(values ~ ind, long_h2ax_prim)
summary(model)
Anova(model)
em1 <- emmeans(model, ~ ind)
pairs(em1)

simulationOutput <- simulateResiduals(fittedModel = model, plot = T)

kruskal <- kruskal.test(values ~ ind, data = long_h2ax_prim)
kruskal
pairwise.wilcox.test(long_h2ax_prim$values, long_h2ax_prim$ind, p.adjust.method = "BH")

ggplot(mutate(long_h2ax_prim, ind = fct_relevel(ind, "X18.yo", "X29.yo", "X32.yo", "X59.yo")), aes(ind, values, fill = ind, color = ind)) + 
  geom_boxplot(outlier.size = 0.7) + 
  labs(y = "Average H2AX staining intensity") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(fill = FALSE, color = FALSE, y = guide_axis(minor.ticks = TRUE)) +
  scale_x_discrete(labels=c('18 y/o', '29 y/o', '32 y/o', '59 y/o')) +
  scale_y_continuous(minor_breaks = scales::breaks_width(2)) 


# bar plot of H2AX staining in primary cells

ggplot(mutate(summary_long_prim, ind = fct_relevel(ind, "X18.yo", "X29.yo", "X32.yo", "X59.yo")), aes(ind, mean, fill = ind, color = ind)) +
  geom_bar(stat = "identity", width = 0.8) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.1, colour="black", alpha=0.9, size=0.5) +
  labs(x = "Age", y = "Average H2AX staining intensity") +
  theme_bw() +
  theme(legend.position="none") +
  scale_color_prism(palette = "floral") +
  scale_fill_prism(palette = "floral") +
  guides(y = guide_axis(minor.ticks = TRUE), x = guide_axis(minor.ticks = TRUE)) 
