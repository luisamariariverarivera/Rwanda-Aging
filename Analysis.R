
library(GGally)
library(tidyverse)
library(tidyr)
library(table1)
library(performance)
library(dplyr)
library(flextable)
library(broom)
library (viridis)
library(factoextra)

d <- read.csv("data_raw/d.csv")

##reorder factor so control is the reference
d$group_factor <- factor(d$group, 
                         levels = c("Non exposed","Exposed to genocide","Exposed to genocide and rape"),
                         labels = c("Control", "Single Exposed","Double Exposed"))


# Create Table 1 Descriptive Statistics
library(table1)

t1 <- table1(~ factor(bio_sex) + factor(factor_education) + factor(incomecategory) +ace_total | group_factor, data=d)

ft1 <- t1flex(t1) %>% 
  save_as_docx(path="figures_tables/table1.docx")

t2 <- table1(~ realAge+ mAge_Hovath + mAge_Hannum + PhenoAge + dunedin + GrimAgeAccel + YingDamAge + YingAdaptAge | group_factor, data=d)

table2 <- t1flex(t2) %>% 
  save_as_docx(path="figures_tables/clock_age_table.docx")


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


## Save as table

ft1 <- as_flextable(Horvath1)
ft2 <- as_flextable(Horvath2)
ft3 <- as_flextable(Hannum1)
ft4 <- as_flextable(Hannum2)
ft5 <- as_flextable(Pheno1)
ft6 <- as_flextable(Pheno2)

save_as_docx(
  `Horvath Model 1` = ft1, `Horvath Model 2` = ft2,`Hannum Model 1` = ft3, `Hannum Model 2` = ft4,`PhenoAge Model 1` = ft5, `PhenoAge Model 2` = ft6,
  path ="figures_tables/First_gen_clock_models.docx")

  ft1 <- as_flextable(YingDamAge1)
  ft2 <- as_flextable(YingDamAge2)
  ft3 <- as_flextable(YingAdaptAge1)
  ft4 <- as_flextable(YingAdaptAge2)
  ft5 <- as_flextable(Grim1)
  ft6 <- as_flextable(Grim1)
  ft7 <- as_flextable(Dunedin1)  
  ft8 <- as_flextable(Dunedin2)
  
  
  save_as_docx(
    `YingDamAge  Model 1` = ft1, `YingDamAge Model 2` = ft2,`YingAdaptAge Model 1` = ft3, `YingAdaptAge Model 2` = ft4,`GrimAgeAccel Model 1` = ft5, `GrimAgeAccel Model 2` = ft6,
    `DunedinPACE Model 1` = ft7, `DunedinPACE Model 2` = ft8,
    path ="figures_tables/second_gen_models.docx")

  
  ##Calculate standardized mean difference and plot across clocks and exposure groups
  
  ##Unadjusted
  
tidyHorvath1 <- tidy(Horvath1)
tidy_Horvath1single_exposed <- tidyHorvath1 %>%
    filter(term == "group_factorSingle Exposed")
tidyHannum1 <- tidy(Hannum1)
tidy_Hannum1single_exposed <- tidyHannum1 %>%
  filter(term == "group_factorSingle Exposed")
tidyPheno1 <- tidy(Pheno1)
tidy_Pheno1single_exposed <- tidyPheno1 %>%
  filter(term == "group_factorSingle Exposed")
tidyGrim1 <- tidy(Grim1)
tidy_Grim1single_exposed <- tidyGrim1 %>%
  filter(term == "group_factorSingle Exposed")
tidyDunedin1 <- tidy(Dunedin1)
tidy_Dunedin1single_exposed <- tidyDunedin1 %>%
  filter(term == "group_factorSingle Exposed")
tidyDamage1 <- tidy(YingDamAge1)
tidy_Damage1single_exposed <- tidyDamage1 %>%
  filter(term == "group_factorSingle Exposed")
tidyAdaptage1 <- tidy(YingAdaptAge1)
tidy_Adaptage1single_exposed <- tidyAdaptage1 %>%
  filter(term == "group_factorSingle Exposed")
tidy_Horvath1double_exposed <- tidyHorvath1 %>%
  filter(term == "group_factorDouble Exposed")
tidy_Hannum1double_exposed <- tidyHannum1 %>%
  filter(term == "group_factorDouble Exposed")
tidy_Pheno1double_exposed <- tidyPheno1 %>%
  filter(term == "group_factorDouble Exposed")
tidy_Grim1double_exposed <- tidyGrim1 %>%
  filter(term == "group_factorDouble Exposed")
tidy_Dunedin1double_exposed <- tidyDunedin1 %>%
  filter(term == "group_factorDouble Exposed")
tidy_Damage1double_exposed <- tidyDamage1 %>%
  filter(term == "group_factorDouble Exposed")
tidy_Adaptage1double_exposed <- tidyAdaptage1 %>%
  filter(term == "group_factorDouble Exposed")
  
combined_tidy <- bind_rows(
  mutate(tidy_Horvath1single_exposed, model = "Horvath1"),
  mutate(tidy_Horvath1double_exposed, model = "Horvath1"),
  mutate(tidy_Hannum1single_exposed, model = "Hannum1"),
  mutate(tidy_Hannum1double_exposed, model = "Hannum1"),
  mutate(tidy_Pheno1single_exposed, model = "Pheno1"),
  mutate(tidy_Pheno1double_exposed, model = "Pheno1"),
  mutate(tidy_Grim1single_exposed, model = "Grim1"),
  mutate(tidy_Grim1double_exposed, model = "Grim1"),
  mutate(tidy_Dunedin1single_exposed, model = "Dunedin1"),
  mutate(tidy_Dunedin1double_exposed, model = "Dunedin1"),
  mutate(tidy_Damage1single_exposed, model = "YingDamAge1"),
  mutate(tidy_Damage1double_exposed, model = "YingDamAge1"),
  mutate(tidy_Adaptage1single_exposed, model = "YingAdaptAge1"),
  mutate(tidy_Adaptage1double_exposed, model = "YingAdaptAge1")
)


d_sd <- data.frame(
  sd = c(
    sd(d$residuals_dunedin), sd(d$residuals_Hannum), sd(d$residuals_Horvath), sd(d$residuals_PhenoAge), sd(d$residuals_YingAdaptAge), sd(d$residuals_YingDamAge), sd(d$grim)
  ),
  
  model = c(
    "Dunedin1", "Hannum1", "Horvath1", "Pheno1", "YingAdaptAge1", "YingDamAge1", "Grim1"
  )
)

d_combined <- left_join(combined_tidy, d_sd)

d_combined$model_2 <- substr(d_combined$model, 1, nchar(d_combined$model)-1)
d_combined$model_2 <- factor(d_combined$model_2, levels = c("Horvath", "Hannum", "Pheno", "Grim", "Dunedin", "YingAdaptAge", "YingDamAge"))

d_combined$term2 <- substr(d_combined$term, nchar("group_factor")+1, nchar(d_combined$term))

# make effect size dataset
d_combined_export1 <- d_combined %>% 
  group_by(model_2, term2) %>% 
  mutate(
    std_mean_diff = estimate/sd,
    lower_95_CI = (estimate - std.error*1.96)/sd,
    upper_95_CI = (estimate + std.error*1.96)/sd
  )


p_noadj <- ggplot(d_combined_export1, aes(x = estimate/sd, y = model_2, color = model_2)) +
  facet_wrap(~fct_rev(term2)) +
  geom_point(position = position_dodge(width = 0.3), size = 3) +  
  geom_errorbarh(aes(xmin = (estimate - std.error*1.96)/sd, xmax = (estimate + std.error*1.96)/sd),
                 position = position_dodge(width = 0.3), height = 0) +  
  geom_vline(xintercept = 0, linetype = "dashed", color = "indianred") +
  labs(title = "Model 1",
       x = "Standardized Mean Difference",
       y = "") +
  scale_color_viridis_d(option = "viridis", begin = 0, end = 1) +
  theme_minimal(base_size = 14) +
  theme(plot.margin = margin(10, 50, 10, 10), legend.title = element_blank(), legend.position = "none", panel.spacing = unit(2, "lines")) +
  scale_x_continuous(breaks = c(-1, 0, 1)) 

p_noadj

# export the estimates
write_csv(d_combined_export1 %>% select(model_2, term2, std_mean_diff, lower_95_CI, upper_95_CI), file = "effect_sizes_unadjusted.csv")


table(d_combined$estimate/d_combined$sd)

## Redo this all but now for the models adjusted for ACEs

tidyHorvath2 <- tidy(Horvath2)
tidy_Horvath2single_exposed <- tidyHorvath2 %>%
  filter(term == "group_factorSingle Exposed")
tidyHannum2 <- tidy(Hannum2)
tidy_Hannum2single_exposed <- tidyHannum2 %>%
  filter(term == "group_factorSingle Exposed")
tidyPheno2 <- tidy(Pheno2)
tidy_Pheno2single_exposed <- tidyPheno2 %>%
  filter(term == "group_factorSingle Exposed")
tidyGrim2 <- tidy(Grim2)
tidy_Grim2single_exposed <- tidyGrim2 %>%
  filter(term == "group_factorSingle Exposed")
tidyDunedin2 <- tidy(Dunedin2)
tidy_Dunedin2single_exposed <- tidyDunedin2 %>%
  filter(term == "group_factorSingle Exposed")
tidyDamage2 <- tidy(YingDamAge2)
tidy_Damage2single_exposed <- tidyDamage2 %>%
  filter(term == "group_factorSingle Exposed")
tidyAdaptage2 <- tidy(YingAdaptAge2)
tidy_Adaptage2single_exposed <- tidyAdaptage2 %>%
  filter(term == "group_factorSingle Exposed")
tidy_Horvath2double_exposed <- tidyHorvath2 %>%
  filter(term == "group_factorDouble Exposed")
tidy_Hannum2double_exposed <- tidyHannum2 %>%
  filter(term == "group_factorDouble Exposed")
tidy_Pheno2double_exposed <- tidyPheno2 %>%
  filter(term == "group_factorDouble Exposed")
tidy_Grim2double_exposed <- tidyGrim2 %>%
  filter(term == "group_factorDouble Exposed")
tidy_Dunedin2double_exposed <- tidyDunedin2 %>%
  filter(term == "group_factorDouble Exposed")
tidy_Damage2double_exposed <- tidyDamage2 %>%
  filter(term == "group_factorDouble Exposed")
tidy_Adaptage2double_exposed <- tidyAdaptage2 %>%
  filter(term == "group_factorDouble Exposed")

combined_tidy2 <- bind_rows(
  mutate(tidy_Horvath2single_exposed, model = "Horvath2"),
  mutate(tidy_Horvath2double_exposed, model = "Horvath2"),
  mutate(tidy_Hannum2single_exposed, model = "Hannum2"),
  mutate(tidy_Hannum2double_exposed, model = "Hannum2"),
  mutate(tidy_Pheno2single_exposed, model = "Pheno2"),
  mutate(tidy_Pheno2double_exposed, model = "Pheno2"),
  mutate(tidy_Grim2single_exposed, model = "Grim2"),
  mutate(tidy_Grim2double_exposed, model = "Grim2"),
  mutate(tidy_Dunedin2single_exposed, model = "Dunedin2"),
  mutate(tidy_Dunedin2double_exposed, model = "Dunedin2"),
  mutate(tidy_Damage2single_exposed, model = "YingDamAge2"),
  mutate(tidy_Damage2double_exposed, model = "YingDamAge2"),
  mutate(tidy_Adaptage2single_exposed, model = "YingAdaptAge2"),
  mutate(tidy_Adaptage2double_exposed, model = "YingAdaptAge2")
)


d_sd <- data.frame(
  sd = c(
    sd(d$residuals_dunedin), sd(d$residuals_Hannum), sd(d$residuals_Horvath), sd(d$residuals_PhenoAge), sd(d$residuals_YingAdaptAge), sd(d$residuals_YingDamAge), sd(d$grim)
  ),
  
  model = c(
    "Dunedin2", "Hannum2", "Horvath2", "Pheno2", "YingAdaptAge2", "YingDamAge2", "Grim2"
  )
)

d_combined2 <- left_join(combined_tidy2, d_sd)


d_combined2$model_2 <- substr(d_combined2$model, 1, nchar(d_combined2$model)-1)
d_combined2$model_2 <- factor(d_combined2$model_2, levels = c("Horvath", "Hannum", "Pheno", "Grim", "Dunedin", "YingAdaptAge", "YingDamAge"))

d_combined2$term2 <- substr(d_combined$term, nchar("group_factor")+1, nchar(d_combined$term))

# make effect size dataset
d_combined_export2 <- d_combined2 %>% 
  group_by(model_2, term2) %>% 
  mutate(
    std_mean_diff = estimate/sd,
    lower_95_CI = (estimate - std.error*1.96)/sd,
    upper_95_CI = (estimate + std.error*1.96)/sd
    )

p_adj <- ggplot(d_combined_export2, aes(x = std_mean_diff, y = model_2, color = model_2)) +
  facet_wrap(~fct_rev(term2)) +
  geom_point(position = position_dodge(width = 0.3), size = 3) +  
  geom_errorbarh(aes(xmin = lower_95_CI, xmax = upper_95_CI),
                 position = position_dodge(width = 0.3), height = 0) +  # 
  geom_vline(xintercept = 0, linetype = "dashed", color = "indianred") +
  labs(title = "Model 2",
       x = "Standardized Mean Difference",
       y = "") +
  scale_color_viridis_d(option = "viridis", begin = 0, end = 1) +
  theme_minimal(base_size = 14) +
  theme(plot.margin = margin(10, 50, 10, 10), legend.title = element_blank(), legend.position = "none", panel.spacing = unit(2, "lines", ))  

# export the estimates
write_csv(d_combined_export2 %>% select(model_2, term2, std_mean_diff, lower_95_CI, upper_95_CI), file = "effect_sizes_adjusted.csv")



# compose the two plots together
library(patchwork)
p_noadj + p_adj

############Supplementary plots###############
  
  ## Immune cell pca
  
  res.pca <- prcomp(subset_vars, scale = F)
  fviz_eig(res.pca)
  
  fviz_pca_var(res.pca,
               col.var = "contrib", # Color by contributions to the PC
               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
               title = "Immune Cell Principle Components Analysis",
               repel = TRUE     # Avoid text overlapping
  )
  
  

  
  fviz_pca_ind(res.pca,
               col.ind = "cos2", 
               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
               repel = TRUE     
  )
  
  
  ##ACEs by exposure group ###################
  
  shapiro.test(d$ace_total)## quick test to see if they are normally distributed between groups, they're not.
  
  ## Descriptive distribution plot for all clocks by group ##########
  d_clock_long <- d %>% 
    pivot_longer(cols = c("residuals_PhenoAge", "residuals_Hannum", "residuals_Horvath", "residuals_YingDamAge", "residuals_YingAdaptAge", "grim", residuals_dunedin), names_to = "clock", values_to = "age_residual") %>% 
    mutate(clock = ifelse(clock %in% c("grim"), clock, sapply( str_split(clock, "_"), "[", 2 )))
  

  
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
          axis.text.x = element_blank())
  
  all_clock_plot
  
  
  d_clock_long3 <- d %>% 
    pivot_longer(cols = c("DunedinPACE",), names_to = "clock", values_to = "age_residual") 
  
  #Make a separate one for Dunedin because its range is much smaller
  dunedin_plot <- ggplot(d_clock_long3, aes(x = group_factor, y = age_residual, fill = group_factor)) +
    facet_wrap(~clock) +
    geom_violin(trim = FALSE, color = NA, alpha = 0.7) +
    geom_jitter(width = 0.2, size = 2, alpha = 0.5) +
    scale_fill_viridis_d() +
    scale_color_viridis_d() +
    labs(x = "", y = "Pace of aging") +
    theme_minimal(18) +
    ## geom_text(aes(label = studyid), position = position_jitter(width = 0.2, height = 0), vjust = -0.7, size = 3)+ ## de-comment if you want to see studyids labelled
    theme(legend.title = element_blank(),
          axis.text.x = element_blank())
  
dunedin_plot
  
  ### Epigenetic age prediction overlapping density plot ######
  
  mean_realAge <- mean(d$realAge)
  
  mean(d$realAge)
  sd(d$realAge)
  
  clock_density <- ggplot(d, aes(x = residuals_Horvath, fill = "Horvath")) +
    geom_density(alpha = 0.5) +
    geom_density(aes(x = residuals_Hannum, fill = "Hannum"), alpha = 0.5) +
    geom_density(aes(x = residuals_PhenoAge, fill = "PhenoAge"), alpha = 0.5) +
    geom_density(aes(x = residuals_YingAdaptAge, fill = "YingAdaptAge"), alpha = 0.5) +
    geom_density(aes(x = residuals_YingDamAge, fill = "YingDamAge"), alpha = 0.5) +
    geom_density(aes(x = residuals_dunedin, fill = "Dunedin Pace"), alpha = 0.5) +
    annotate("text", x = 30, y = 0.2, label = "Chronological Age (mean = 24.12, sd = .10)", vjust = -1, color = "black") + 
    labs(title = "Age residuals",
         x = "Predicted Age",
         y = "Density",
         fill = "Clock Type") +
    scale_fill_manual(values = c("Horvath" = "blue", "Hannum" = "green",  "PhenoAge" = "red", "YingAdaptAge" = "cyan", "YingDamAge" = "purple" )) +
    theme_minimal()
  
  clock_density
  
  
  #### Pairs plot for clock correlations with each other 
  library(GGally)
  clocks <- data.frame(Hannum = d$mAge_Hannum, Horvath = d$mAge_Hovath, PhenoAge = d$PhenoAge, YingAdaptAge = d$YingAdaptAge, YingDamAge = d$YingDamAge, GrimAgeAccel =d$grim, DunedinPACE=d$dunedin)  
  ggpairs(clocks, title="Correlations between  Epigenetic Clocks") 
  
  round(cor(clocks), 3)
  
## scatterplot
  
  
  correlation <- cor.test(d$PhenoAge, d$Age)
  cor_value <- correlation$estimate
  p_value <- correlation$p.value
  formatted_p_value <- ifelse(p_value < 0.001, "< 0.001", format(round(p_value, 3)))
  

  ggplot(d, aes(x = Age, y = PhenoAge)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE, color = "blue") +
    labs(title = "",
         x = "Chronological Age",
         y = "PhenoAge",
         caption = paste("Correlation:", round(cor_value, 2), "| p-value:", formatted_p_value)) +
    theme_minimal()
  
  
  correlation <- cor.test(d$PhenoAge, d$Age)
  cor_value <- correlation$estimate
  p_value <- correlation$p.value
  formatted_p_value <- ifelse(p_value < 0.001, "< 0.001", format(round(p_value, 3)))
  
  
 phenocor <-  ggplot(d, aes(x = Age, y = PhenoAge)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE, color = "blue") +
    labs(title = "",
         x = "Chronological Age",
         y = "PhenoAge",
         caption = paste("Correlation:", round(cor_value, 2), "| p-value:", formatted_p_value)) +
    theme_minimal()
  
 
 
 correlation <- cor.test(d$mAge_Hannum, d$Age)
 cor_value <- correlation$estimate
 p_value <- correlation$p.value
 formatted_p_value <- ifelse(p_value < 0.001, "< 0.001", format(round(p_value, 3)))
 
hannumcor <-  ggplot(d, aes(x = Age, y = mAge_Hannum)) +
   geom_point() +
   geom_smooth(method = "lm", se = FALSE, color = "red") +
   labs(title = "",
        x = "Chronological Age",
        y = "Hannum Age",
        caption = paste("Correlation:", round(cor_value, 2), "| p-value:", formatted_p_value)) +
   theme_minimal()


 
correlation <- cor.test(d$mAge_Hovath, d$Age)
cor_value <- correlation$estimate
p_value <- correlation$p.value
formatted_p_value <- ifelse(p_value < 0.001, "< 0.001", format(round(p_value, 3)))

horvathcor <-  ggplot(d, aes(x = Age, y = mAge_Hovath)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "darkgreen") +
  labs(title = "",
       x = "Chronological Age",
       y = "Horvath Age",
       caption = paste("Correlation:", round(cor_value, 2), "| p-value:", formatted_p_value)) +
  theme_minimal()

horvathcor

library(patchwork)
firstcor <- horvathcor+hannumcor+phenocor

print(firstcor)



correlation <- cor.test(d$YingAdaptAge, d$Age)
cor_value <- correlation$estimate
p_value <- correlation$p.value
formatted_p_value <- ifelse(p_value < 0.001, "< 0.001", format(round(p_value, 3)))


adaptcor <- ggplot(d, aes(x = Age, y = YingAdaptAge)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "lightgreen") +
  labs(title = "",
       x = "Chronological Age",
       y = "YingAdaptAge",
       caption = paste("Correlation:", round(cor_value, 2), "| p-value:", formatted_p_value)) +
  theme_minimal()


correlation <- cor.test(d$YingDamAge, d$Age)
cor_value <- correlation$estimate
p_value <- correlation$p.value
formatted_p_value <- ifelse(p_value < 0.001, "< 0.001", format(round(p_value, 3)))


damagecor <-  ggplot(d, aes(x = Age, y = YingDamAge)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "lightblue") +
  labs(title = "",
       x = "Chronological Age",
       y = "YingDamAge",
       caption = paste("Correlation:", round(cor_value, 2), "| p-value:", formatted_p_value)) +
  theme_minimal()


library(patchwork)
secondcor <- damagecor+adaptcor

print(secondcor)

allcor <- firstcor+secondcor
print(allcor)



#Immune cell pca visualization

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


  