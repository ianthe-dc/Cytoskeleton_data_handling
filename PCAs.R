# PCA analyses

library(tidyverse)
library(lme4)
library(car)
library(emmeans)
library(ggprism)
library(forcats)
library(dplyr)
library(ggsignif)
library(DHARMa)
library(factoextra)


#ACTIN

actin_df <- read.csv("merged_output_actin.csv")

pc_columns_old <- c("mean_length_2d", "sum_length_2d", "mean_avg_angle", "mean_net_curv", "mean_acc_curv", "mean_length_3d", "sum_length_3d", "mean_deviation", "filament_count", "sum_branches", "mean_branches", "Cell.Volume..um.3.", "Skeleton.Density", "Branch.Point.Density...um.3.", "Average.Filament.Width..um.")

colnames(actin_df)[colnames(actin_df) == "mean_length_2d"] <- "Mean length (2D)"
colnames(actin_df)[colnames(actin_df) == "mean_length_3d"] <- "Mean length (3D)"
colnames(actin_df)[colnames(actin_df) == "sum_length_2d"] <- "Total length (2D)"
colnames(actin_df)[colnames(actin_df) == "sum_length_3d"] <- "Total length (3D)"
colnames(actin_df)[colnames(actin_df) == "mean_avg_angle"] <- "Average filament angle"
colnames(actin_df)[colnames(actin_df) == "mean_acc_curv"] <- "Total curvature"
colnames(actin_df)[colnames(actin_df) == "mean_net_curv"] <- "Net curvature"
colnames(actin_df)[colnames(actin_df) == "mean_deviation"] <- "Deviation"
colnames(actin_df)[colnames(actin_df) == "filament_count"] <- "Filament count"
colnames(actin_df)[colnames(actin_df) == "sum_branches"] <- "Branch count"
colnames(actin_df)[colnames(actin_df) == "mean_branches"] <- "Branch ratio"
colnames(actin_df)[colnames(actin_df) == "Cell.Volume..um.3."] <- "Cell volume (µm³)"
colnames(actin_df)[colnames(actin_df) == "Skeleton.Density"] <- "Network density"
colnames(actin_df)[colnames(actin_df) == "Branch.Point.Density...um.3."] <- "Branch density (/µm³)"
colnames(actin_df)[colnames(actin_df) == "Average.Filament.Width..um."] <- "Average filament width (μm)"

pc_columns <- c("Mean length (2D)", "Mean length (3D)", "Total length (2D)", "Total length (3D)", "Average filament angle", "Total curvature", "Net curvature", "Deviation", "Filament count", "Branch count", "Branch ratio", "Cell volume (µm³)", "Network density", "Branch density (/µm³)", "Average filament width (μm)")


actin_primary <- filter(actin_df, condition_clean %in% c("primary fibro 18yo act", "primary fibro 29yo act", "primary fibro 32yo act", "primary fibro 59yo act"))
actin_UV <- filter(actin_df, condition_clean %in% c("bj p12.5 -uv act", "bj p12.5 +uv act"))
actin_etop <- filter(actin_df, condition_clean %in% c("bj p8.5 act", "bj p9.5 etop vim"))
actin_pass <- filter(actin_df, condition_clean %in% c("bj p11.5 act", "bj p26.5 act", "bj p31.5 act", "bj p36.5 act", "bj p38.5 act", "bj p40.5 act"))

# PCA for serial passage data

pcApass <- princomp(na.omit(actin_pass[, pc_columns ]),cor=TRUE, scores=TRUE)
summary(pcApass)
biplot(pcApass)
pcApass$loadings

fviz_pca_ind(pcApass,
                habillage = factor(actin_pass$condition_clean, labels = c("p11.5", "p26.5", "p31.5", "p36.5", "p38.5", "p40.5")),
                label = "none",
                addEllipses = FALSE,
                pointshape = 19,      
                pointsize = 1.3) +
  labs(colour = "Passage")
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/PCAs/new/pcApassInd.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")

fviz_eig(pcApass)

fviz_pca_ind(pcApass,
             axes = c(3, 4),
             habillage = actin_pass$condition_clean,
             geom.ind = "point",
             pointshape = 19,
             pointsize = 2)

fviz_pca_var(pcApass,
             axes = c(3, 4),
             col.var = contrib[,4],
             gradient.cols = c("blue", "yellow", "red"),
             repel = TRUE) +
  labs(colour = "Contribution to PC3–PC4 (%)")


fviz_pca_var(pcApass,
             repel = TRUE) 
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/PCAs/new/pcApassVar.tiff", device = "tiff", dpi = 300, width = 202.5, height = 150, units = "mm")

var <- get_pca_var(pcApass)

contrib <- var$contrib

fviz_pca_var(pcApass,
             col.var = contrib[,1],   # PC1 contributions
             gradient.cols = c("blue", "yellow", "red"),
             repel = TRUE) +
  labs(colour = "Contribution to PC1 (%)")
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/PCAs/pcApassVar1.tiff", device = "tiff", dpi = 300, width = 202.5, height = 150, units = "mm")

fviz_pca_var(pcApass,
             col.var = contrib[,2],   # PC2 contributions
             gradient.cols = c("blue", "yellow", "red"),
             repel = TRUE) +
  labs(colour = "Contribution to PC2 (%)")
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/PCAs/pcApassVar2.tiff", device = "tiff", dpi = 300, width = 202.5, height = 150, units = "mm")


# PCA for UV data

pcAUV <- princomp(na.omit(actin_UV[, pc_columns ]),cor=TRUE, scores=TRUE)
summary(pcAUV)
biplot(pcAUV)
pcAUV$loadings

fviz_pca_ind(pcAUV,
             habillage = factor(actin_UV$condition_clean, labels = c("p12.5 -UV", "p12.5 +UV")),
             label = "none",
             addEllipses = FALSE,
             pointshape = 19,      
             pointsize = 1.3) +
  labs(colour = "Condition")
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/PCAs/new/pcAUVInd.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")


fviz_pca_var(pcAUV,
             repel = TRUE) 
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/PCAs/new/pcAUVVar.tiff", device = "tiff", dpi = 300, width = 202.5, height = 150, units = "mm")

fviz_pca_var(pcAUV,
             col.var = contrib[,1],   # PC1 contributions
             gradient.cols = c("blue", "yellow", "red"),
             repel = TRUE) +
  labs(colour = "Contribution to PC1 (%)")
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/PCAs/pcAUVVar1.tiff", device = "tiff", dpi = 300, width = 202.5, height = 150, units = "mm")

fviz_pca_var(pcAUV,
             col.var = contrib[,2],   # PC2 contributions
             gradient.cols = c("blue", "yellow", "red"),
             repel = TRUE) +
  labs(colour = "Contribution to PC2 (%)")
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/PCAs/pcAUVVar2.tiff", device = "tiff", dpi = 300, width = 202.5, height = 150, units = "mm")


# PCA for etoposide data

pcAetop <- princomp(na.omit(actin_etop[, pc_columns ]),cor=TRUE, scores=TRUE)
summary(pcAetop)
biplot(pcAetop)
pcAetop$loadings

fviz_pca_ind(pcAetop,
             habillage = factor(actin_etop$condition_clean, labels = c("p8.5 -etop", "p9.5 +etop")),
             label = "none",
             addEllipses = FALSE,
             pointshape = 19,      
             pointsize = 1.3) +
  labs(colour = "Condition")
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/PCAs/new/pcAetopInd.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")


fviz_pca_var(pcAetop,
             repel = TRUE) 
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/PCAs/new/pcAetopVar.tiff", device = "tiff", dpi = 300, width = 202.5, height = 150, units = "mm")

fviz_pca_var(pcAetop,
             col.var = contrib[,1],   # PC1 contributions
             gradient.cols = c("blue", "yellow", "red"),
             repel = TRUE) +
  labs(colour = "Contribution to PC1 (%)")
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/PCAs/pcAEtopVar1.tiff", device = "tiff", dpi = 300, width = 202.5, height = 150, units = "mm")

fviz_pca_var(pcAetop,
             col.var = contrib[,2],   # PC2 contributions
             gradient.cols = c("blue", "yellow", "red"),
             repel = TRUE) +
  labs(colour = "Contribution to PC2 (%)")
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/PCAs/pcAEtopVar2.tiff", device = "tiff", dpi = 300, width = 202.5, height = 150, units = "mm")


# PCA for primary fibroblast data

pcAprim <- princomp(na.omit(actin_primary[, pc_columns ]),cor=TRUE, scores=TRUE)
summary(pcAprim)
biplot(pcAprim)
pcAprim$loadings

fviz_pca_ind(pcAprim,
             habillage = factor(actin_primary$condition_clean, labels = c("18 y/o", "29 y/o", "32 y/o", "59 y/o")),
             label = "none",
             addEllipses = FALSE,
             pointshape = 19,      
             pointsize = 1.3) +
  labs(colour = "Age")
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/PCAs/new/pcAprimInd.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")


fviz_pca_var(pcAprim,
             repel = TRUE) 
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/PCAs/new/pcAprimVar.tiff", device = "tiff", dpi = 300, width = 202.5, height = 150, units = "mm")

fviz_pca_var(pcAprim,
             col.var = contrib[,1],   # PC1 contributions
             gradient.cols = c("blue", "yellow", "red"),
             repel = TRUE) +
  labs(colour = "Contribution to PC1 (%)")
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/PCAs/pcAprimVar1.tiff", device = "tiff", dpi = 300, width = 202.5, height = 150, units = "mm")

fviz_pca_var(pcAprim,
             col.var = contrib[,2],   # PC2 contributions
             gradient.cols = c("blue", "yellow", "red"),
             repel = TRUE) +
  labs(colour = "Contribution to PC2 (%)")
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/PCAs/pcAprimVar2.tiff", device = "tiff", dpi = 300, width = 202.5, height = 150, units = "mm")



#VIMENTIN

vimentin_df <- read.csv("merged_output_vimentin.csv")

pc_columns_old <- c("mean_length_2d", "sum_length_2d", "mean_avg_angle", "mean_net_curv", "mean_acc_curv", "mean_length_3d", "sum_length_3d", "mean_deviation", "filament_count", "sum_branches", "mean_branches", "Cell.Volume..um.3.", "Skeleton.Density", "Branch.Point.Density...um.3.", "Average.Filament.Width..um.")

colnames(vimentin_df)[colnames(vimentin_df) == "mean_length_2d"] <- "Mean length (2D)"
colnames(vimentin_df)[colnames(vimentin_df) == "mean_length_3d"] <- "Mean length (3D)"
colnames(vimentin_df)[colnames(vimentin_df) == "sum_length_2d"] <- "Total length (2D)"
colnames(vimentin_df)[colnames(vimentin_df) == "sum_length_3d"] <- "Total length (3D)"
colnames(vimentin_df)[colnames(vimentin_df) == "mean_avg_angle"] <- "Average filament angle"
colnames(vimentin_df)[colnames(vimentin_df) == "mean_acc_curv"] <- "Total curvature"
colnames(vimentin_df)[colnames(vimentin_df) == "mean_net_curv"] <- "Net curvature"
colnames(vimentin_df)[colnames(vimentin_df) == "mean_deviation"] <- "Deviation"
colnames(vimentin_df)[colnames(vimentin_df) == "filament_count"] <- "Filament count"
colnames(vimentin_df)[colnames(vimentin_df) == "sum_branches"] <- "Branch count"
colnames(vimentin_df)[colnames(vimentin_df) == "mean_branches"] <- "Branch ratio"
colnames(vimentin_df)[colnames(vimentin_df) == "Cell.Volume..um.3."] <- "Cell volume (µm³)"
colnames(vimentin_df)[colnames(vimentin_df) == "Skeleton.Density"] <- "Network density"
colnames(vimentin_df)[colnames(vimentin_df) == "Branch.Point.Density...um.3."] <- "Branch density (/µm³)"
colnames(vimentin_df)[colnames(vimentin_df) == "Average.Filament.Width..um."] <- "Average filament width (μm)"

pc_columns <- c("Mean length (2D)", "Mean length (3D)", "Total length (2D)", "Total length (3D)", "Average filament angle", "Total curvature", "Net curvature", "Deviation", "Filament count", "Branch count", "Branch ratio", "Cell volume (µm³)", "Network density", "Branch density (/µm³)", "Average filament width (μm)")


vimentin_primary <- filter(vimentin_df, condition_clean %in% c("primary fibro 18yo vim", "primary fibro 29yo vim", "primary fibro 32yo vim", "primary fibro 59yo vim"))
vimentin_UV <- filter(vimentin_df, condition_clean %in% c("bj p12.5 -uv vim", "bj p12.5 +uv vim"))
vimentin_etop <- filter(vimentin_df, condition_clean %in% c("bj p8.5 vim", "bj p9.5 etop vim"))
vimentin_pass <- filter(vimentin_df, condition_clean %in% c("bj p11.5 vim", "bj p26.5 vim", "bj p31.5 vim", "bj p36.5 vim", "bj p38.5 vim", "bj p40.5 vim"))

# PCA for serial passage data

pcVpass <- princomp(na.omit(vimentin_pass[, pc_columns ]),cor=TRUE, scores=TRUE)
summary(pcVpass)
biplot(pcVpass)
pcVpass$loadings

fviz_pca_ind(pcVpass,
             habillage = factor(vimentin_pass$condition_clean, labels = c("p11.5", "p26.5", "p31.5", "p36.5", "p38.5", "p40.5")),
             label = "none",
             addEllipses = FALSE,
             pointshape = 19,      
             pointsize = 1.3) +
  labs(colour = "Passage")
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/PCAs/new/pcVpassInd.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")


fviz_pca_var(pcVpass,
             repel = TRUE) 
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/PCAs/new/pcVpassVar.tiff", device = "tiff", dpi = 300, width = 202.5, height = 150, units = "mm")

fviz_pca_var(pcVpass,
             col.var = contrib[,1],   # PC1 contributions
             gradient.cols = c("blue", "yellow", "red"),
             repel = TRUE) +
  labs(colour = "Contribution to PC1 (%)")
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/PCAs/pcVpassVar1.tiff", device = "tiff", dpi = 300, width = 202.5, height = 150, units = "mm")

fviz_pca_var(pcVpass,
             col.var = contrib[,2],   # PC2 contributions
             gradient.cols = c("blue", "yellow", "red"),
             repel = TRUE) +
  labs(colour = "Contribution to PC2 (%)")
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/PCAs/pcVpassVar2.tiff", device = "tiff", dpi = 300, width = 202.5, height = 150, units = "mm")


# PCA for UV data

pcVUV <- princomp(na.omit(vimentin_UV[, pc_columns ]),cor=TRUE, scores=TRUE)
summary(pcVUV)
biplot(pcVUV)
pcVUV$loadings

fviz_pca_ind(pcVUV,
             habillage = factor(vimentin_UV$condition_clean, labels = c("p12.5 -UV", "p12.5 +UV")),
             label = "none",
             addEllipses = FALSE,
             pointshape = 19,      
             pointsize = 1.3) +
  labs(colour = "Condition")
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/PCAs/new/pcVUVInd.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")


fviz_pca_var(pcVUV,
             repel = TRUE) 
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/PCAs/new/pcVUVVar.tiff", device = "tiff", dpi = 300, width = 202.5, height = 150, units = "mm")

fviz_pca_var(pcVUV,
             col.var = contrib[,1],   # PC1 contributions
             gradient.cols = c("blue", "yellow", "red"),
             repel = TRUE) +
  labs(colour = "Contribution to PC1 (%)")
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/PCAs/pcVUVVar1.tiff", device = "tiff", dpi = 300, width = 202.5, height = 150, units = "mm")

fviz_pca_var(pcVUV,
             col.var = contrib[,2],   # PC2 contributions
             gradient.cols = c("blue", "yellow", "red"),
             repel = TRUE) +
  labs(colour = "Contribution to PC2 (%)")
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/PCAs/pcVUVVar2.tiff", device = "tiff", dpi = 300, width = 202.5, height = 150, units = "mm")


# PCA for etoposide data

pcVetop <- princomp(na.omit(vimentin_etop[, pc_columns ]),cor=TRUE, scores=TRUE)
summary(pcVetop)
biplot(pcVetop)
pcVetop$loadings

fviz_pca_ind(pcVetop,
             habillage = factor(vimentin_etop$condition_clean, labels = c("p8.5 -etop", "p9.5 +etop")),
             label = "none",
             addEllipses = FALSE,
             pointshape = 19,      
             pointsize = 1.3) +
  labs(colour = "Condition")
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/PCAs/new/pcVetopInd.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")


fviz_pca_var(pcVetop,
             repel = TRUE) 
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/PCAs/new/pcVetopVar.tiff", device = "tiff", dpi = 300, width = 202.5, height = 150, units = "mm")

fviz_pca_var(pcVetop,
             col.var = contrib[,1],   # PC1 contributions
             gradient.cols = c("blue", "yellow", "red"),
             repel = TRUE) +
  labs(colour = "Contribution to PC1 (%)")
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/PCAs/pcVEtopVar1.tiff", device = "tiff", dpi = 300, width = 202.5, height = 150, units = "mm")

fviz_pca_var(pcVetop,
             col.var = contrib[,2],   # PC2 contributions
             gradient.cols = c("blue", "yellow", "red"),
             repel = TRUE) +
  labs(colour = "Contribution to PC2 (%)")
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/PCAs/pcVEtopVar2.tiff", device = "tiff", dpi = 300, width = 202.5, height = 150, units = "mm")


# PCA for primary fibroblast data

pcVprim <- princomp(na.omit(vimentin_primary[, pc_columns ]),cor=TRUE, scores=TRUE)
summary(pcVprim)
biplot(pcVprim)
pcVprim$loadings

fviz_pca_ind(pcVprim,
             habillage = factor(vimentin_primary$condition_clean, labels = c("18 y/o", "29 y/o", "32 y/o", "59 y/o")),
             label = "none",
             addEllipses = FALSE,
             pointshape = 19,      
             pointsize = 1.3) +
  labs(colour = "Age")
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/PCAs/new/pcVprimInd.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")


fviz_pca_var(pcVprim,
             repel = TRUE) 
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/PCAs/new/pcVprimVar.tiff", device = "tiff", dpi = 300, width = 202.5, height = 150, units = "mm")

fviz_pca_var(pcVprim,
             col.var = contrib[,1],   # PC1 contributions
             gradient.cols = c("blue", "yellow", "red"),
             repel = TRUE) +
  labs(colour = "Contribution to PC1 (%)")
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/PCAs/pcVprimVar1.tiff", device = "tiff", dpi = 300, width = 202.5, height = 150, units = "mm")

fviz_pca_var(pcVprim,
             col.var = contrib[,2],   # PC2 contributions
             gradient.cols = c("blue", "yellow", "red"),
             repel = TRUE) +
  labs(colour = "Contribution to PC2 (%)")
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/PCAs/pcVprimVar2.tiff", device = "tiff", dpi = 300, width = 202.5, height = 150, units = "mm")



#TUBULIN

tubulin_df <- read.csv("merged_output_tubulin.csv")

pc_columns_old <- c("mean_length_2d", "sum_length_2d", "mean_avg_angle", "mean_net_curv", "mean_acc_curv", "mean_length_3d", "sum_length_3d", "mean_deviation", "filament_count", "sum_branches", "mean_branches", "Cell.Volume..um.3.", "Skeleton.Density", "Branch.Point.Density...um.3.", "Average.Filament.Width..um.")

colnames(tubulin_df)[colnames(tubulin_df) == "mean_length_2d"] <- "Mean length (2D)"
colnames(tubulin_df)[colnames(tubulin_df) == "mean_length_3d"] <- "Mean length (3D)"
colnames(tubulin_df)[colnames(tubulin_df) == "sum_length_2d"] <- "Total length (2D)"
colnames(tubulin_df)[colnames(tubulin_df) == "sum_length_3d"] <- "Total length (3D)"
colnames(tubulin_df)[colnames(tubulin_df) == "mean_avg_angle"] <- "Average filament angle"
colnames(tubulin_df)[colnames(tubulin_df) == "mean_acc_curv"] <- "Total curvature"
colnames(tubulin_df)[colnames(tubulin_df) == "mean_net_curv"] <- "Net curvature"
colnames(tubulin_df)[colnames(tubulin_df) == "mean_deviation"] <- "Deviation"
colnames(tubulin_df)[colnames(tubulin_df) == "filament_count"] <- "Filament count"
colnames(tubulin_df)[colnames(tubulin_df) == "sum_branches"] <- "Branch count"
colnames(tubulin_df)[colnames(tubulin_df) == "mean_branches"] <- "Branch ratio"
colnames(tubulin_df)[colnames(tubulin_df) == "Cell.Volume..um.3."] <- "Cell volume (µm³)"
colnames(tubulin_df)[colnames(tubulin_df) == "Skeleton.Density"] <- "Network density"
colnames(tubulin_df)[colnames(tubulin_df) == "Branch.Point.Density...um.3."] <- "Branch density (/µm³)"
colnames(tubulin_df)[colnames(tubulin_df) == "Average.Filament.Width..um."] <- "Average filament width (μm)"

pc_columns <- c("Mean length (2D)", "Mean length (3D)", "Total length (2D)", "Total length (3D)", "Average filament angle", "Total curvature", "Net curvature", "Deviation", "Filament count", "Branch count", "Branch ratio", "Cell volume (µm³)", "Network density", "Branch density (/µm³)", "Average filament width (μm)")


tubulin_primary <- filter(tubulin_df, condition_clean %in% c("primary fibro 18yo tub", "primary fibro 29yo tub", "primary fibro 32yo tub", "primary fibro 59yo tub"))
tubulin_UV <- filter(tubulin_df, condition_clean %in% c("bj p12.5 -uv tub", "bj p12.5 +uv tub"))
tubulin_pass <- filter(tubulin_df, condition_clean %in% c("bj p11.5 tub", "bj p26.5 tub", "bj p31.5 tub", "bj p36.5 tub", "bj p38.5 tub", "bj p40.5 tub"))

# PCA for serial passage data

pcTpass <- princomp(na.omit(tubulin_pass[, pc_columns ]),cor=TRUE, scores=TRUE)
summary(pcTpass)
biplot(pcTpass)
pcTpass$loadings

fviz_pca_ind(pcTpass,
             habillage = factor(tubulin_pass$condition_clean, labels = c("p11.5", "p26.5", "p31.5", "p36.5", "p38.5", "p40.5")),
             label = "none",
             addEllipses = FALSE,
             pointshape = 19,      
             pointsize = 1.3) +
  labs(colour = "Passage")
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/PCAs/new/pcTpassInd.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")


fviz_pca_var(pcTpass,
             repel = TRUE) 
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/PCAs/new/pcTpassVar.tiff", device = "tiff", dpi = 300, width = 202.5, height = 150, units = "mm")


fviz_pca_var(pcTpass,
             col.var = contrib[,1],   # PC1 contributions
             gradient.cols = c("blue", "yellow", "red"),
             repel = TRUE) +
  labs(colour = "Contribution to PC1 (%)")
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/PCAs/pcTpassVar1.tiff", device = "tiff", dpi = 300, width = 202.5, height = 150, units = "mm")

fviz_pca_var(pcTpass,
             col.var = contrib[,2],   # PC2 contributions
             gradient.cols = c("blue", "yellow", "red"),
             repel = TRUE) +
  labs(colour = "Contribution to PC2 (%)")
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/PCAs/pcTpassVar2.tiff", device = "tiff", dpi = 300, width = 202.5, height = 150, units = "mm")


# PCA for UV data

pcTUV <- princomp(na.omit(tubulin_UV[, pc_columns ]),cor=TRUE, scores=TRUE)
summary(pcTUV)
biplot(pcTUV)
pcTUV$loadings

fviz_pca_ind(pcTUV,
             habillage = factor(tubulin_UV$condition_clean, labels = c("p12.5 -UV", "p12.5 +UV")),
             label = "none",
             addEllipses = FALSE,
             pointshape = 19,      
             pointsize = 1.3) +
  labs(colour = "Condition")
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/PCAs/new/pcTUVInd.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")


fviz_pca_var(pcTUV,
             repel = TRUE) 
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/PCAs/new/pcTUVVar.tiff", device = "tiff", dpi = 300, width = 202.5, height = 150, units = "mm")

fviz_pca_var(pcTUV,
             col.var = contrib[,1],   # PC1 contributions
             gradient.cols = c("blue", "yellow", "red"),
             repel = TRUE) +
  labs(colour = "Contribution to PC1 (%)")
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/PCAs/pcTUVVar1.tiff", device = "tiff", dpi = 300, width = 202.5, height = 150, units = "mm")

fviz_pca_var(pcTUV,
             col.var = contrib[,2],   # PC2 contributions
             gradient.cols = c("blue", "yellow", "red"),
             repel = TRUE) +
  labs(colour = "Contribution to PC2 (%)")
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/PCAs/pcTUVVar2.tiff", device = "tiff", dpi = 300, width = 202.5, height = 150, units = "mm")


# PCA for primary fibroblast data

pcTprim <- princomp(na.omit(tubulin_primary[, pc_columns ]),cor=TRUE, scores=TRUE)
summary(pcTprim)
biplot(pcTprim)
pcTprim$loadings

fviz_pca_ind(pcTprim,
             habillage = factor(tubulin_primary$condition_clean, labels = c("18 y/o", "29 y/o", "32 y/o", "59 y/o")),
             label = "none",
             addEllipses = FALSE,
             pointshape = 19,      
             pointsize = 1.3) +
  labs(colour = "Age")
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/PCAs/new/pcTprimInd.tiff", device = "tiff", dpi = 300, width = 135, height = 100, units = "mm")


fviz_pca_var(pcTprim,
             repel = TRUE) 
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/PCAs/new/pcTprimVar.tiff", device = "tiff", dpi = 300, width = 202.5, height = 150, units = "mm")


fviz_pca_var(pcTprim,
             col.var = contrib[,1],   # PC1 contributions
             gradient.cols = c("blue", "yellow", "red"),
             repel = TRUE) +
  labs(colour = "Contribution to PC1 (%)")
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/PCAs/pcTprimVar1.tiff", device = "tiff", dpi = 300, width = 202.5, height = 150, units = "mm")

fviz_pca_var(pcTprim,
             col.var = contrib[,2],   # PC2 contributions
             gradient.cols = c("blue", "yellow", "red"),
             repel = TRUE) +
  labs(colour = "Contribution to PC2 (%)")
ggsave("/Volumes/T7/Fixed data quantitative analysis/final graphs/PCAs/new/pcTprimVar2.tiff", device = "tiff", dpi = 300, width = 202.5, height = 150, units = "mm")


