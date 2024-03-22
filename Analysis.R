
install.packages("sesame", "sesameData")
library(tidyverse)
library(tidyr)
library(table1)
library(sesame)
library(sesameData)
library(data.table)
library(GenomicRanges)
library(performance)
library(dplyr)
library(flextable)


##read in phenotypic, dried bloodspot card, immune cell-type, and aging data

d <- read.csv("data_raw/covariates_outcomes.csv")
e <- read.csv("data_raw/DNAm_real_age.csv")
f <- read.csv("data_raw/Pred_EPIC.csv")
g <- read.csv("data_raw/idats/SampleSheet_DBSCard_2023.csv")
h <- read.csv("data_raw/d_pace.csv")
i <- read.csv("data_raw/clock_preds_rev.csv")

d <- left_join(d, g, by=c("studyid" = "Study_ID"))
d <- left_join(d,e, by = c("Sample_Name" = "SampleID"))
d <- left_join(d, h, by =c("Filename" = "id"))
d <- left_join(d,f, by= c("Sample_Name" = "Description"))
d <- left_join(d,i, by= c("Filename" = "id"))

list(d)


rm(e,f,g,h,i)

par(mfrow = c(1, 5))

# Plot histograms
hist(d$YingAdaptAge, main = "YingAdaptAge Histogram")
hist(d$YingCausAge, main = "YingCausAge Histogram")
hist(d$YingDamAge, main = "YingDamAge Histogram")
hist(d$grim, main = "GrimAge Histogram")
hist(d$dunedin, main = "Dunedin Pace Histogram")


##set factors
d$group <- as.factor(d$group)
d$exposed <- ifelse(d$group=="Non exposed", 0, 1)
d$bio_sex <- as.factor(d$bio_sex)
d$incomecategory <- as.factor(d$incomecategory)
d$education <- as.factor(d$education)
d$ptsd <- as.factor(d$ptsd)
d$group_factor<- factor(d$group, levels = c("Non exposed", "Exposed to genocide", "Exposed to genocide and rape"))
# Rename levels of the group_factor variable
d$group_factor <- factor(d$group_factor, 
                         levels = c("Non exposed", "Exposed to genocide", "Exposed to genocide and rape"),
                         labels = c("Control", "Singly Exposed", "Doubly Exposed"))



## Make immune cell principle components and calculate MLR and NLR

d$mlr <- d$Mono/(d$Bmem+d$Bnv+d$CD4mem+d$CD4nv+d$CD8mem+d$CD8nv+d$NK+d$Treg) ## monocyte-lymphocyte ratio
d$nlr <- d$Neu/(d$Bmem+d$Bnv+d$CD4mem+d$CD4nv+d$CD8mem+d$CD8nv+d$NK+d$Treg) ## neutrophil-lymphocyte ratio

subset_vars <- d[, c("Bas", "Bmem", "Bnv", "CD4mem", "CD4nv", "CD8nv", "Eos", "Mono", "Neu", "NK", "Treg")]


pca <- prcomp(subset_vars)


# Principal components
pcs <- pca$x

# Proportion of variance explained
var_explained <- pca$sdev^2 / sum(pca$sdev^2)

summary(pca)
pca$rotation


d <- cbind(d, pcs)


##calculate and save residuals from predicted vs. chronological age

m1 <- lm(d$mAge_Hannum ~ d$Age)
d$residuals_Hannum <- residuals(m1)
m2 <- lm(d$mAge_Hovath ~ d$Age)
d$residuals_Horvath <- residuals(m2)
m3 <- lm(d$PhenoAge ~ d$Age)
d$residuals_PhenoAge <- residuals(m3)
m4 <- lm(d$YingCausAge ~ d$Age)
d$residuals_YingCausAge <- residuals(m4)
m5 <- lm(d$YingAdaptAge ~ d$Age)
d$residuals_YingAdaptAge <- residuals(m5)
m6 <- lm(d$YingDamAge ~ d$Age)
d$residuals_YingDamAge <- residuals(m6)
m8 <- lm(d$dunedin ~ d$Age)
d$residuals_dunedin <- residuals(m8)

## take a look at them ##
plot(d$residuals_Hannum)
plot(d$residuals_Horvath)
plot(d$residuals_PhenoAge)
plot(d$residuals_YingCausAge)
plot(d$residuals_YingAdaptAge)
plot(d$residuals_YingDamAge)
plot(d$residuals_dunedin)



write.csv(d,"data_raw/d.csv")

d <- read.csv("data_raw/d.csv")

view(d$group_factor)


########### Descriptive Plots ###########

# Create Table 1 Descriptive Statistics
library(table1)

t1 <- table1(~ factor(bio_sex) + factor(education) + factor(incomecategory) +factor(ptsd) +ace_total | group_factor, data=d)

ft1 <- t1flex(t1) %>% 
  save_as_docx(path="figures_tables/table1.docx")

t2 <- table1(~ realAge+ mAge_Hovath + mAge_Hannum + PhenoAge + dunedin + grim + YingDamAge + YingAdaptAge | group_factor, data=d)

table2 <- t1flex(t2) %>% 
  save_as_docx(path="figures_tables/clock_age_table.docx")


#Immune cell pca

res.pca <- prcomp(subset_vars, scale = F)
fviz_eig(res.pca)

fviz_pca_var(res.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             title = "Immune Cell Principle Components Analysis",
             repel = TRUE     # Avoid text overlapping
)




fviz_pca_ind(res.pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

#########Descriptive plots ##################

##ACEs by exposure group ###################

shapiro.test(d$ace_total)## quick test to see if they are normally distributed between groups, they're not.

ggplot(d, aes(x = ace_total, fill= "pink")) +
  geom_density(binwidth = 1, position = "identity", alpha = 0.7) +
  labs(x = "ace_total", y = "Count") +
  ggtitle("Histogram of ACEs ") +
  theme_minimal() +
  theme(legend.position = "none",
        panel.border = element_blank())

ggplot(d, aes(x = ace_total, fill = group_factor)) +
  geom_density(binwidth = 1, position = "identity", alpha = 0.7) +
  labs(x = "ace_total", y = "Count") +
  ggtitle("Histogram of ACEs by Group") +
  facet_wrap(~ group_factor, ncol = 1) +
  theme_minimal() +
  theme(legend.position = "none",
        panel.border = element_blank())

### Violin plots of immune outcomes ########

library(viridisLite)

fill_cols <- viridis(3)
pc1plot <- ggplot(d, aes(x = group_factor, y = PC1, fill = group_factor)) +
  geom_violin(trim = FALSE, color = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.5) +
  scale_fill_viridis_d() +
  scale_color_viridis_d() +
  labs(x = "", y = "PC1") +
  ggtitle("Violin Plot of PC1 by Group") +
  theme_minimal(18) +
  theme(legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels at 45 degrees

pc1plot



mlrplot <- ggplot(d, aes(x = group_factor, y = mlr, fill = group_factor)) +
  geom_violin(trim = FALSE, color = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.5) +
  scale_fill_viridis_d() +
  scale_color_viridis_d() +
  labs(x = "", y = "Monocyte-Lymphocyte Ratio") +
  ggtitle("Violin Plot of Monocyte-Lymphocyte Ratio by Group") +
  theme_minimal(18) +
  theme(legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels at 45 degrees

mlrplot

nlrplot <- ggplot(d, aes(x = group_factor, y = nlr, fill = group_factor)) +
  geom_violin(trim = FALSE, color = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.5) +
  scale_fill_viridis_d() +
  scale_color_viridis_d() +
  labs(x = "", y = "Neutrophil-Lymphocyte Ratio") +
  ggtitle("Violin Plot of Neutrophil-Lymphocyte Ratio by Group") +
  theme_minimal(18) +
  theme(legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels at 45 degrees

nlrplot

## Violin Plots of  Epigenetic Aging by Group ########
# 
# ##phenoplot <- ggplot(d, aes(x = group_factor, y = residuals_PhenoAge, fill = group_factor)) +
#   geom_violin(trim = FALSE, color = NA, alpha = 0.7) +
#   geom_jitter(width = 0.2, size = 2, alpha = 0.5) +
#   scale_fill_viridis_d() +
#   scale_color_viridis_d() +
#   labs(x = "", y = "Age Residual") +
#   ggtitle("Violin Plot of PhenoAge Residual by Group") +
#   theme_minimal(18) +
#   ##geom_text(aes(label = studyid), position = position_jitter(width = 0.2, height = 0), vjust = -0.7, size = 3)+ ## de-comment if you want to see studyids labelled
#   theme(legend.title = element_blank(),
#         axis.text.x = element_text(angle = 45, hjust = 1))
# print(phenoplot)
# 
# horvathplot <- ggplot(d, aes(x = group_factor, y = residuals_Horvath, fill = group_factor)) +
#   geom_violin(trim = FALSE, color = NA, alpha = 0.7) +
#   geom_jitter(width = 0.2, size = 2, alpha = 0.5) +
#   scale_fill_viridis_d() +
#   scale_color_viridis_d() +
#   labs(x = "", y = "Age Residual") +
#   ggtitle("Violin Plot of Horvath Residual by Group") +
#   theme_minimal(18) +
#   ##geom_text(aes(label = studyid), position = position_jitter(width = 0.2, height = 0), vjust = -0.7, size = 3)+## de-comment if you want to see studyids labelled
#   theme(legend.title = element_blank(),
#         axis.text.x = element_text(angle = 45, hjust = 1))
# print(horvathplot)
# 
# 
# hannumplot <- ggplot(d, aes(x = group, y = residuals_Hannum, fill = group)) +
#   geom_violin(trim = FALSE, color = NA, alpha = 0.7) +
#   geom_jitter(width = 0.2, size = 2, alpha = 0.5) +
#   scale_fill_viridis_d() +
#   scale_color_viridis_d() +
#   labs(x = "", y = "Age Residual") +
#   ggtitle("Violin Plot of Hannum Residual by Group") +
#   theme_minimal(18) +
#  ## geom_text(aes(label = studyid), position = position_jitter(width = 0.2, height = 0), vjust = -0.7, size = 3)+## de-comment if you want to see studyids labelled
#   theme(legend.title = element_blank(),
#         axis.text.x = element_text(angle = 45, hjust = 1))
# 
# print(hannumplot)
# 
# dunedinplot <- ggplot(d, aes(x = group, y = dunedin, fill = group)) +
#   geom_violin(trim = FALSE, color = NA, alpha = 0.7) +
#   geom_jitter(width = 0.2, size = 2, alpha = 0.5) +
#   scale_fill_viridis_d() +
#   scale_color_viridis_d() +
#   labs(x = "", y = "Dunedin Pace of Aging") +
#   ggtitle("Violin Plot of Dunedin Pace of Aging by Group") +
#   theme_minimal(18) +
#   ## geom_text(aes(label = studyid), position = position_jitter(width = 0.2, height = 0), vjust = -0.7, size = 3)+## de-comment if you want to see studyids labelled
#   theme(legend.title = element_blank(),
#         axis.text.x = element_text(angle = 45, hjust = 1))
# 
# print(dunedinplot)
# 
# 
# causageplot <- ggplot(d, aes(x = group, y = YingCausAge, fill = group)) +
#   geom_violin(trim = FALSE, color = NA, alpha = 0.7) +
#   geom_jitter(width = 0.2, size = 2, alpha = 0.5) +
#   scale_fill_viridis_d() +
#   scale_color_viridis_d() +
#   labs(x = "", y = "Ying Causal Predicted Age") +
#   ggtitle("Violin Plot of Causal Aging Group") +
#   theme_minimal(18) +
#   ## geom_text(aes(label = studyid), position = position_jitter(width = 0.2, height = 0), vjust = -0.7, size = 3)+## de-comment if you want to see studyids labelled
#   theme(legend.title = element_blank(),
#         axis.text.x = element_text(angle = 45, hjust = 1))
# 
# print(causageplot)
# 
# adaptageplot <- ggplot(d, aes(x = group, y = YingAdaptAge, fill = group)) +
#   geom_violin(trim = FALSE, color = NA, alpha = 0.7) +
#   geom_jitter(width = 0.2, size = 2, alpha = 0.5) +
#   scale_fill_viridis_d() +
#   scale_color_viridis_d() +
#   labs(x = "", y = "Ying Adaptive Predicted Age") +
#   ggtitle("Violin Plot of Adaptive Aging Group") +
#   theme_minimal(18) +
#   ## geom_text(aes(label = studyid), position = position_jitter(width = 0.2, height = 0), vjust = -0.7, size = 3)+## de-comment if you want to see studyids labelled
#   theme(legend.title = element_blank(),
#         axis.text.x = element_text(angle = 45, hjust = 1))
# 
# print(adaptageplot)
# 
# 
# damageplot <- ggplot(d, aes(x = group, y = YingDamAge, fill = group)) +
#   geom_violin(trim = FALSE, color = NA, alpha = 0.7) +
#   geom_jitter(width = 0.2, size = 2, alpha = 0.5) +
#   scale_fill_viridis_d() +
#   scale_color_viridis_d() +
#   labs(x = "", y = "Ying Damage Predicted Age") +
#   ggtitle("Violin Plot of Damage Aging Group") +
#   theme_minimal(18) +
#   ## geom_text(aes(label = studyid), position = position_jitter(width = 0.2, height = 0), vjust = -0.7, size = 3)+## de-comment if you want to see studyids labelled
#   theme(legend.title = element_blank(),
#         axis.text.x = element_text(angle = 45, hjust = 1))
# 
# print(damageplot)
# 
# 
# grimplot <- ggplot(d, aes(x = group, y = grim, fill = group)) +
#   geom_violin(trim = FALSE, color = NA, alpha = 0.7) +
#   geom_jitter(width = 0.2, size = 2, alpha = 0.5) +
#   scale_fill_viridis_d() +
#   scale_color_viridis_d() +
#   labs(x = "", y = "GrimAge Predicted Age") +
#   ggtitle("Violin Plot of GrimAge Group") +
#   theme_minimal(18) +
#   ## geom_text(aes(label = studyid), position = position_jitter(width = 0.2, height = 0), vjust = -0.7, size = 3)+## de-comment if you want to see studyids labelled
#   theme(legend.title = element_blank(),
#         axis.text.x = element_text(angle = 45, hjust = 1))
# 
# print(grimplot)##
# 
# 
# 


## Combine plot for all clocks by group ##########
d_clock_long <- d %>% 
  pivot_longer(cols = c("residuals_PhenoAge", "residuals_Hannum", "residuals_Horvath"), names_to = "clock", values_to = "age_residual") %>% 
  mutate(clock = sapply( str_split(clock, "_"), "[", 2 ))

all_clock_plot <- ggplot(d_clock_long, aes(x = group_factor, y = age_residual, fill = group_factor)) +
  facet_wrap(~clock) +
  geom_violin(trim = FALSE, color = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.5) +
  scale_fill_viridis_d() +
  scale_color_viridis_d() +
  labs(x = "", y = "Age Residual") +
  theme_minimal(18) +
  ## geom_text(aes(label = studyid), position = position_jitter(width = 0.2, height = 0), vjust = -0.7, size = 3)+ ## de-comment if you want to see studyids labelled
  theme(legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

all_clock_plot


d_clock_long2 <- d %>% 
  pivot_longer(cols = c( "residuals_YingDamAge", "residuals_YingAdaptAge"), names_to = "clock", values_to = "age_residual") %>% 
  mutate(clock = sapply( str_split(clock, "_"), "[", 2 ))

all_clock_plot2 <- ggplot(d_clock_long2, aes(x = group_factor, y = age_residual, fill = group_factor)) +
  facet_wrap(~clock) +
  geom_violin(trim = FALSE, color = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.5) +
  scale_fill_viridis_d() +
  scale_color_viridis_d() +
  labs(x = "", y = "Age Residual") +
  theme_minimal(18) +
  ## geom_text(aes(label = studyid), position = position_jitter(width = 0.2, height = 0), vjust = -0.7, size = 3)+ ## de-comment if you want to see studyids labelled
  theme(legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

all_clock_plot2

d_clock_long3 <- d %>% 
  pivot_longer(cols = c("residuals_dunedin",), names_to = "clock", values_to = "age_residual") %>% 
  mutate(clock = sapply( str_split(clock, "_"), "[", 2 ))

all_clock_plot3 <- ggplot(d_clock_long3, aes(x = group_factor, y = age_residual, fill = group_factor)) +
  facet_wrap(~clock) +
  geom_violin(trim = FALSE, color = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.5) +
  scale_fill_viridis_d() +
  scale_color_viridis_d() +
  labs(x = "", y = "Age Residual") +
  theme_minimal(18) +
  ## geom_text(aes(label = studyid), position = position_jitter(width = 0.2, height = 0), vjust = -0.7, size = 3)+ ## de-comment if you want to see studyids labelled
  theme(legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

all_clock_plot3

d_clock_long4 <- d %>% 
  pivot_longer(cols = c("grim",), names_to = "clock", values_to = "age_residual") 

all_clock_plot4 <- ggplot(d_clock_long4, aes(x = group_factor, y = age_residual, fill = group_factor)) +
  facet_wrap(~clock) +
  geom_violin(trim = FALSE, color = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.5) +
  scale_fill_viridis_d() +
  scale_color_viridis_d() +
  labs(x = "", y = "Age Residual") +
  theme_minimal(18) +
  ## geom_text(aes(label = studyid), position = position_jitter(width = 0.2, height = 0), vjust = -0.7, size = 3)+ ## de-comment if you want to see studyids labelled
  theme(legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

all_clock_plot4


d_clock_long5 <- d %>% 
  pivot_longer(cols = c("DNAmTL",), names_to = "clock", values_to = "DNAm_telomere_length") 

all_clock_plot5 <- ggplot(d_clock_long5, aes(x = group_factor, y = DNAm_telomere_length, fill = group_factor)) +
  facet_wrap(~clock) +
  geom_violin(trim = FALSE, color = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.5) +
  scale_fill_viridis_d() +
  scale_color_viridis_d() +
  labs(x = "", y = "DNAm Telomere Length") +
  theme_minimal(18) +
  ## geom_text(aes(label = studyid), position = position_jitter(width = 0.2, height = 0), vjust = -0.7, size = 3)+ ## de-comment if you want to see studyids labelled
  theme(legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

all_clock_plot5

### Epigenetic age prediction overlapping density plot ######


mean_realAge <- mean(d$realAge)

mean(d$realAge)
sd(d$realAge)

clock_density <- ggplot(d, aes(x = mAge_Hovath, fill = "Horvath Clock")) +
  geom_density(alpha = 0.5) +
  geom_density(aes(x = mAge_Hannum, fill = "Hannum Clock"), alpha = 0.5) +
  geom_density(aes(x = PhenoAge, fill = "PhenoAge Clock"), alpha = 0.5) +
  geom_density(aes(x = YingAdaptAge, fill = "Adaptive Age Clock"), alpha = 0.5) +
  geom_density(aes(x = YingDamAge, fill = "Damage Age Clock"), alpha = 0.5) +
  geom_vline(xintercept = mean_realAge, linetype = "solid", color = "black", linewidth = 1) + 
  annotate("text", x = 30, y = 0.2, label = "Chronological Age (mean = 24.12, sd = .10)", vjust = -1, color = "black") + 
  labs(title = "Predicted Epigenetic Ages of Study Participants",
       x = "Predicted Age",
       y = "Density",
       fill = "Clock Type") +
  scale_fill_manual(values = c("Horvath Clock" = "blue", "Hannum Clock" = "green",  "PhenoAge Clock" = "red", "Causal Age Clock" = "pink", "Adaptive Age Clock" = "cyan", "Damage Age Clock" = "purple" )) +
  theme_minimal()

clock_density


#### Pairs plot for clock correlations with each other 
library(GGally)
clocks <- data.frame(Hannum = d$mAge_Hannum, Horvath = d$mAge_Hovath, PhenoAge = d$PhenoAge, AdaptiveAge = d$YingAdaptAge, DamageAge = d$YingDamAge, GrimAge =d$grim, DunedinPace=d$dunedin, DNAmTL = d$DNAmTL )  
ggpairs(clocks, title="Correlations between Predicted Ages of Epigenetic Clocks") 

round(cor(clocks), 3)





## Immune Linear models: MLR, NLR, and PC1 #################################
# 
# mlr1 <- lm(mlr ~ group_factor + bio_sex, data = d)
# mlr2 <- lm(mlr ~ group_factor+ ace_total + bio_sex, data = d)
# 
# summary(mlr1)
# summary(mlr2)
# 
# ## remove outliers
# 
# 
# remove_outliers <- function(d, z_threshold = 3) {
#   d %>%
#     mutate(z_score_mlr = abs(scale(mlr))) %>%
#     filter(z_score_mlr <= z_threshold) %>%
#     select(-z_score_mlr)
# }
# 
# # Remove outliers from the data using the function
# data_no_outliers <- remove_outliers(d)
# 
# mlr1.1 <- lm(mlr ~ group_factor + bio_sex, data = data_no_outliers)
# mlr2.1 <- lm(mlr ~ group_factor+ ace_total + bio_sex, data = data_no_outliers)
# 
# summary(mlr1.1)
# summary(mlr2.1)
# 
# 
# nlr1 <- lm(nlr ~ group_factor + bio_sex, data = d)
# nlr2 <- lm(nlr ~ group_factor + ace_total + bio_sex, data = d)
# 
# summary(nlr1)
# summary(nlr2)
# 
# remove_outliers <- function(d, z_threshold = 3) {
#   d %>%
#     mutate(z_score_nlr = abs(scale(nlr))) %>%
#     filter(z_score_nlr <= z_threshold) %>%
#     select(-z_score_nlr)
# }
# 
# # Remove outliers from the data using the function
# data_no_outliers <- remove_outliers(d)
# 
# nlr1.1 <- lm(nlr ~ group_factor + bio_sex, data = data_no_outliers)
# nlr2.1 <- lm(nlr ~ group_factor+ ace_total + bio_sex, data = data_no_outliers)
# 
# summary(nlr1.1)
# summary(nlr2.1)
# 
# pc1 <- lm(PC1 ~ group_factor +  bio_sex, data = d)
# pc2 <- lm(PC1 ~ group_factor + ace_total+ bio_sex, data = d)
# 
# summary(pc1)
# summary(pc2)
# 
# remove_outliers <- function(d, z_threshold = 3) {
#   d %>%
#     mutate(z_score_PC1 = abs(scale(PC1))) %>%
#     filter(z_score_PC1 <= z_threshold) %>%
#     select(-z_score_PC1)
# }
# 
# # Remove outliers from the data using the function
# data_no_outliers <- remove_outliers(d)
# 
# pc1.1 <- lm(PC1 ~ group_factor + bio_sex, data = data_no_outliers)
# pc2.1 <- lm(PC1 ~ group_factor+ ace_total + bio_sex, data = data_no_outliers)
# 
# summary(pc1.1)
# summary(pc2.1)
# 
# 
# ## Save model outputs as tables in Word###########
# 
# ftmlr1 <- as_flextable(mlr1)
# ftmlr2 <- as_flextable(mlr2)
# ftnlr2 <- as_flextable(nlr2)
# ftnlr1 <- as_flextable(nlr1)
# ftpc1 <- as_flextable(pc1)  
# ftpc2 <- as_flextable(pc2)
# 
# save_as_docx(
#   `Monocyte-Lymphocyte Ratio Model 1` = ftmlr1, `Monocyte-Lymphocyte Ratio Model 2` = ftmlr2, `Neutrophil-Lymphocyte Ratio Model 1` = ftnlr1,`Neutrophil-Lymphocyte Ratio Model 2` = ftnlr2,`Immune PC1 Model 1` = ftpc1, `Immune PC1 Model 2` = ftpc2, path ="figures_tables/ImmuneModelv2.docx")
# 
# 
# 


## Epigenetic aging by group linear models#######################

Horvath1 <- lm(residuals_Horvath ~ group_factor + bio_sex +PC1, data = d)
Horvath2 <- lm(residuals_Horvath ~ ace_total + group_factor + bio_sex+ PC1, data = d)


Hannum1 <- lm(residuals_Hannum ~ group_factor + bio_sex+ PC1, data = d)
Hannum2 <- lm(residuals_Hannum ~ ace_total + group_factor + bio_sex+PC1, data = d)


Pheno1 <- lm(residuals_PhenoAge ~ group_factor + bio_sex+PC1, data = d)
Pheno2 <- lm(residuals_PhenoAge ~ ace_total + group_factor + bio_sex+PC1, data = d)


YingAdaptAge1 <- lm(residuals_YingAdaptAge~ group_factor + bio_sex+PC1, data = d)
YingAdaptAge2 <- lm(residuals_YingAdaptAge ~ ace_total + group_factor + bio_sex+PC1, data = d)

YingDamAge1 <- lm(residuals_YingDamAge~ group_factor + bio_sex+PC1, data = d)
YingDamAge2 <- lm(residuals_YingDamAge ~ ace_total + group_factor + bio_sex+PC1, data = d)

Dunedin1 <- lm(residuals_dunedin~ group_factor + bio_sex+PC1, data = d)
Dunedin2 <- lm(residuals_dunedin ~ ace_total + group_factor + bio_sex+PC1, data = d)

Grim1 <- lm(grim~ group_factor + bio_sex+PC1, data = d)
Grim2 <- lm(grim ~ ace_total + group_factor + bio_sex+PC1, data = d)

DNAmTL1 <- lm(DNAmTL~ group_factor + bio_sex+PC1, data = d)
DNAmTL2 <- lm(DNAmTL ~ ace_total + group_factor + bio_sex+PC1, data = d)


# Print the summary of the model
summary(Horvath1)
summary(Horvath2)
summary(Hannum1)
summary(Hannum2)
summary(Pheno1)
summary(Pheno2)
summary(YingAdaptAge1)
summary(YingAdaptAge2)
summary(YingDamAge1)
summary(YingDamAge2)
summary(Dunedin1)
summary(Dunedin2)
summary(Grim1)
summary(Grim2)
summary(DNAmTL1)
summary(DNAmTL2)




 # Take a look at model performance (not great)
library(performance)

check_model(Hannum1)
check_model(Hannum2)
check_model(Horvath1)
check_model(Horvath2)
check_model(Pheno1)
check_model(Pheno2)
check_model(YingAdaptAge1)
check_model(YingAdaptAge2)
check_model(YingDamAge1)
check_model(YingDamAge2)

## save as Word tables ###

  ft1 <- as_flextable(Horvath1)
  ft2 <- as_flextable(Horvath2)
  ft3 <- as_flextable(Hannum1)
  ft4 <- as_flextable(Hannum2)
  ft5 <- as_flextable(Pheno1)  
  ft6 <- as_flextable(Pheno2)
  ft
  
  save_as_docx(
    `Horvath Clock Model 1` = ft1, `Horvath Clock Model 2` = ft2,`Hannum Clock Model 1` = ft3, `Hannum Clock Model 2` = ft4,`PhenoAge Clock Model 1` = ft5, `PhenoAge Clock Model 2` = ft6,
    path ="figures_tables/ClockModelv3.docx")

  
  ##sensitivity analysis for Hannum outliers
  
  remove_outliers <- function(d, z_threshold = 3) {
    d %>%
      mutate(z_score_residuals_Hannum = abs(scale(residuals_Hannum))) %>%
      filter(z_score_residuals_Hannum <= z_threshold) %>%
      select(-z_score_residuals_Hannum)
  }
  
  
  # Remove outliers from the data using the function
  data_no_outliers <- remove_outliers(d)
  
  Hannum1.1 <- lm(residuals_Hannum ~ group_factor + bio_sex+ PC1, data = data_no_outliers)
  Hannum2.1 <- lm(residuals_Hannum ~ ace_total + group_factor + bio_sex+PC1, data = data_no_outliers)
  
  summary(Hannum1.1)
  summary(Hannum2.1)

  
  write.csv(d, file ="d.csv")
  