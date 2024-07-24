
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


d <- d %>%
  mutate(factor_education = case_when(
    education %in% c("adv diploma", "bachelors") ~ "university level",
    education %in% c("primary", "some primary") ~ "primary level",
    education %in% c("secondary", "some secondary", "vocation") ~ "secondary level",
    TRUE ~ as.character(education)
  ))

d$ptsd <- as.factor(d$ptsd)
d$group_factor<- factor(d$group, levels = c("Non exposed", "Exposed to genocide", "Exposed to genocide and rape"))
# Rename levels of the group_factor variable
d$group_factor <- factor(d$group_factor, 
                         levels = c("Exposed to genocide and rape","Exposed to genocide", "Non exposed"),
                         labels = c("Double Exposed", "Single Exposed", "Control"))
d$GrimAgeAccel=d$grim
d$DunedinPACE = d$dunedin
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

