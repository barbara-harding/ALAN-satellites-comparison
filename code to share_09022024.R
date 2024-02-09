### LOADING DATA and packages
#################################################################################
#Clean space
rm(list=ls())

#load needed packages
library(haven)
library(dplyr)
library(gam)
library(ggplot2)
library(blandr)
library(gmodels)

## pull in dataset
data <- read_dta("//fs01.isglobal.lan/mcc-spain-ext/04_Analysis/Datos ALAN/03_Data/MCC_ALAN_analysis_requested_variables_20230307.dta") 
data$night_ever <- as.factor(data$night_ever)
summary(data$night_ever)
data$night_ever_recoded <- ifelse(data$night_ever %in% c(6, 0), 0,
                                  ifelse(data$night_ever == 1,1,'NA'))

data$night_ever_recoded <- as.factor(data$night_ever_recoded)
summary(data$night_ever_recoded)
View(data)

##create subset of dataset with COMPLETE CRC case status
data$controls <- as.factor(data$casoc)
data<- data[!is.na(data$casoc),]
View(data)

#we also want to remove those missing data on one or more satellite types
summary(data$impvlggr) 
data <- data[!is.na(data$impvlggr),]
summary(data$impvlggr)
summary(data$msigr)
data <- data[!is.na(data$msigr),]
summary(data$msigr)
summary(data$dmsp11cal_avg_vis)
data <- data[!is.na(data$dmsp11cal_avg_vis),]
summary(data$dmsp11cal_avg_vis)
summary(data$viirs_vcmcfg_avg) 
data <- data[!is.na(data$viirs_vcmcfg_avg),]
summary(data$viirs_vcmcfg_avg)
View(data)
summary(data$controls)
#N=2314
#################################################################################

### SD increase
#################################################################################
# standarise all exposure variables
#imp from CRC paper (same years but better image)
summary(data$impvlggr)
imp_sd <- sd(data$impvlggr)
imp_mean <- mean(data$impvlggr)
data$visual_n <- (data$impvlggr - imp_mean)/imp_sd
summary(data$visual_n)

#mel from CRC paper (same years but better image)
summary(data$msigr)
mel_sd <- sd(data$msigr)
mel_mean <- mean(data$msigr)
data$mel_n <- (data$msigr - mel_mean)/mel_sd
summary(data$mel_n)

summary(data$dmsp11cal_avg_vis)
dmsp11cal_avg_vis_sd <- sd(data$dmsp11cal_avg_vis)
dmsp11cal_avg_vis_mean <- mean(data$dmsp11cal_avg_vis)
data$dmsp11 <- (data$dmsp11cal_avg_vis - dmsp11cal_avg_vis_mean)/dmsp11cal_avg_vis_sd
summary(data$dmsp11)

summary(data$viirs_vcmcfg_avg)
viirs_vcmcfg_avg_sd <- sd(data$viirs_vcmcfg_avg)
viirs_vcmcfg_avg_mean <- mean(data$viirs_vcmcfg_avg)
data$viirs <- (data$viirs_vcmcfg_avg - viirs_vcmcfg_avg_mean)/viirs_vcmcfg_avg_sd
summary(data$viirs)

#################################################################################

### CORRELATION
#################################################################################
variables <- data[, c("visual_n", 
  "mel_n", "dmsp11", "viirs"
)]

variables <- data.matrix(variables) 
cormat <- cor(variables, use="complete.obs")
cormat <- round(cormat,3)

head(cormat)
library(reshape2)
melted_cormat <- melt(cormat)
head(melted_cormat)
library(ggplot2)
ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()

# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

upper_tri <- get_upper_tri(cormat)
upper_tri

# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)
# Heatmap
ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()
# In the figure above:negative correlations are in blue color and positive correlations in red. The function scale_fill_gradient2 is used with the argument limit = c(-1,1) as correlation coefficients range from -1 to 1.
# coord_fixed() : this function ensures that one unit on the x-axis is the same length as one unit on the y-axis.

# reorder correlation matrix
reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}
# Reorder the correlation matrix
cormat <- reorder_cormat(cormat)
upper_tri <- get_upper_tri(cormat)
# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)
# Create a ggheatmap
ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()
# Print the heatmap
print(ggheatmap)

#Add correlation coefficients on the heatmap
ggheatmap + 
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))

##scatter plot##
ggplot(data, aes(x = visual_n, y = viirs)) +
  geom_point() +
  labs(title = "Scatter Plot", x = "iss visual", y = "viirs visual")
ggplot(data, aes(x = visual_n, y = dmsp11)) +
  geom_point() +
  labs(title = "Scatter Plot", x = "iss visual", y = "dmsp visual")

#################################################################################

### DESCRIPTIVE
#################################################################################
library(data.table)
library(tableone)
vardesc<- c("impvlggr", "msigr","visual", "melat", "dmsp11cal_avg_vis", "viirs_vcmcfg_avg")

tab <-CreateTableOne(vardesc, data=data,includeNA =T)
try<-print(tab, nonnormal = c("impvlggr", "msigr","visual", "melat", "dmsp11cal_avg_vis", "viirs_vcmcfg_avg"),  
           noSpaces = TRUE, varLabels=T, includeNA=T)

View(try)
write.csv(try,"//fs01.isglobal.lan/mcc-spain-ext/04_Analysis/Datos ALAN/04_Analysis/Results/Descriptive/Exposures_des_all_08112023.csv")

### tertiles among controls 
data$casoc <- as.factor(data$casoc)
data$casoc_d <- as.factor(data$casoc)
data$casoc_d[is.na(data$casoc_d)] = 0
controls <- data[!data$casoc_d=="1", ]

#CRC paper
data$visual_n_tert <- cut(data$impvlggr, breaks=c(quantile(controls$impvlggr, probs = seq(0, 1, 1/3))), 
                         labels=c("Low","Medium", "High"), include.lowest=TRUE) 
summary(data$visual_n_tert)

#CRC paper
data$mel_n_tert <- cut(data$msigr, breaks=c(quantile(controls$msigr, probs = seq(0, 1, 1/3))), 
                         labels=c("Low","Medium", "High"), include.lowest=TRUE) 
data$mel_n_tert[is.na(data$mel_n_tert)] = "Low"
summary(data$mel_n_tert)

#DMSP11
data$dmsp_11_tert <- cut(data$dmsp11cal_avg_vis, breaks=c(quantile(controls$dmsp11cal_avg_vis, probs = seq(0, 1, 1/3))), 
                             labels=c("Low","Medium", "High"), include.lowest=TRUE)
data$dmsp_11_tert[is.na(data$dmsp_11_tert)] = "Low"
summary(data$dmsp_11_tert)

#VIIRS
data$viirs_tert <- cut(data$viirs_vcmcfg_avg, breaks=c(quantile(controls$viirs_vcmcfg_avg, probs = seq(0, 1, 1/3))), 
                           labels=c("Low","Medium", "High"), include.lowest=TRUE) 
data$viirs_tert[is.na(data$viirs_tert)] = "Low"
summary(data$viirs_tert)


#prepare covariates
library(car)
summary(data$education_basic_final)
data$education_basic_final <- as.factor(data$education_basic_final)
data$education_basic_final <- Recode(data$education_basic_final, '1= "Less than primary"; 2 = "Primary school";
3 = "Secondary school"; 4= "University"', as.factor=TRUE)
summary(data$education_basic_final)

data$score_SE <- as.factor(data$score_SE)
data$score_SE <- Recode(data$score_SE, '1= "Low"; 2 = "Medium";
3 = "High"', as.factor=TRUE)
summary(data$score_SE)

data$fumador_ever <- as.factor(data$fumador_ever)
data$fumador_ever <- Recode(data$fumador_ever, '0= "No"; 1 = "Yes"', as.factor=TRUE)
summary(data$fumador_ever)

data$a3_sexe <- as.factor(data$a3_sexe)
data$a3_sexe <- Recode(data$a3_sexe, '1= "Male"; 2 = "Female"', as.factor=TRUE)
summary(data$a3_sexe)

data$g14_problemas_sue_o <- as.factor(data$g14_problemas_sue_o)
data$g14_problemas_sue_o <- Recode(data$g14_problemas_sue_o, '1= "Yes"; 0 = "No"', as.factor=TRUE)
summary(data$g14_problemas_sue_o)

data$msf_c3 <- as.factor(data$msf_c3)
summary(data$msf_c3)
data$msf_c3 <- Recode(data$msf_c3, '1= "Morning"; 2 = "Neither";
3 = "Evening"', as.factor=TRUE)
summary(data$msf_c3)

data$FH_cancer_Colon <- factor(data$FH_cancer_Colon) #fix
summary(data$FH_cancer_Colon)
data$FH_cancer_Colon[data$FH_cancer_Colon=="2" | data$FH_cancer_Colon=="3" | 
                        data$FH_cancer_Colon=="0" | is.na(data$FH_cancer_Colon)] = "0"
data$FH_cancer_Colon <- factor(data$FH_cancer_Colon) #fix
summary(data$FH_cancer_Colon)

data$night_ever_recoded
summary(data$night_ever_recoded) 

##summarize how those without info on night shift differ from those who work night shift or who do not work night shift
data$night_stratify <- ifelse(is.na(data$night_stratify), 3, data$night_stratify)
summary(data$night_stratify)
data$night_stratify <- factor(ifelse(data$night_stratify == 2, 1, as.numeric(data$night_stratify)))
summary <- table(data$night_stratify)
print(summary)

#genearte table with descriptives based on whether complete info exists on night shift work status or not
vardesc<- c("casoc", "a4_edat", "bmi_complementario", "education_basic_final", "a3_sexe",
            "fumador_ever", "t_ethanolP", "t_sA_energy", "score_SE", "score_WCRF_simple", "g14_problemas_sue_o",
            "g13_tiempo_dormir", "msf_c3", "uvi", "night_ever_recoded")

tab <-CreateTableOne(vardesc, strata = c("night_stratify"), data=data,includeNA =T)
try<-print(tab, noSpaces = TRUE)

View(try)
write.csv(try,"//fs01.isglobal.lan/mcc-spain-ext/04_Analysis/Datos ALAN/04_Analysis/Results/Descriptive/Descriptive_bynightstatus_08022024.csv")

##generate table 1 with descriptives by tertiles of blue light exposure
vardesc<- c("casoc", "a4_edat", "bmi_complementario", "education_basic_final", "a3_sexe",
            "fumador_ever", "t_ethanolP", "t_sA_energy", "score_SE", "score_WCRF_simple", "g14_problemas_sue_o",
            "g13_tiempo_dormir", "msf_c3", "uvi", "score_WCRF_simple", "night_ever_recoded")

tab <-CreateTableOne(vardesc, strata = c("mel_n_tert"), data=data,includeNA =T)
try<-print(tab, noSpaces = TRUE)

View(try)
write.csv(try,"//fs01.isglobal.lan/mcc-spain-ext/04_Analysis/Datos ALAN/04_Analysis/Results/Descriptive/Descriptive_by_MSI_withnights24112023.csv")

#################################################################################

### SCATTER
#################################################################################
library(tidyverse)
par(mfrow=c(3,3)) 
hist(data$impvlggr) 
hist(data$msigr) 
hist(data$visual) 
hist(data$melat) 
hist(data$dmsp11cal_avg_vis) 
hist(data$viirs_vcmcfg_avg) 
#################################################################################

### BA PLOT
#################################################################################
run_ba <- function(data, methoda, methodb) {
  
  # Get values 
  a <- data[[methoda]]
  b <- data[[methodb]]
  mean_diff <- mean(a - b)
  sd_diff <- sd(a - b)
  average <- (a+b)/2
  differences <- a-b
  #gam_model <- gam(differences ~ s(average, k = 3), data = data)
  
  # Create plot
  plot <- ggplot(data = data, aes(x = (a + b) / 2, y = a - b)) +
    geom_point(size = 2) +
    geom_hline(yintercept = mean_diff, linetype = "dashed") +
    geom_hline(yintercept = mean_diff + 1.96 * sd_diff, linetype = "dashed", color = "red") +
    geom_hline(yintercept = mean_diff - 1.96 * sd_diff, linetype = "dashed", color = "red") +
    geom_text(aes(x = max((a + b) / 2), y = mean_diff + 1.96 * sd_diff + 0.1, label = sprintf("%.2f", mean_diff + 1.96 * sd_diff), hjust = 0)) +
    geom_text(aes(x = max((a + b) / 2), y = mean_diff - 1.96 * sd_diff - 0.1, label = sprintf("%.2f", mean_diff - 1.96 * sd_diff), hjust = 0)) +
    geom_smooth(data = data, aes(x = average, y = differences), method = "gam", formula = y ~ s(x, k=3), color = "green", se = TRUE) +
    labs(title = "Bland-Altman Plot with Non-linear Trend",
         x = "Averages",
         y = "Differences") +
    ggtitle("Bland-Altman Plot with Non-linear Trend (N=2,314)") +
    ylab("Difference Between Measurements") +
    xlab("Average Measurement") +
    theme_bw()
    theme(text = element_text(size = 10))
  
  # Display summary
  blandr.display.and.draw(a, b, plotter = "ggplot", method1name = "a",
                          method2name = "b")
  
  return(plot)
  
}

##using gam plots to explore non-linear trends in the BA plots
#title: Bland-Altman Plot with Non-linear Trend for bias including CI using GAM
# DMSP 11 (calibrated) vs melat new
run_ba(data = data, methoda = "dmsp11", methodb = "visual_n")
# VIIRS 13 vs melat new
run_ba(data = data, methoda = "viirs", methodb = "visual_n")
# VISUAL paper PC&BC vs visual paper CRC
##run_ba(data = data, methoda = "visual_n", methodb = "visual_o")
# MELAT paper PC&BC vs melat paper CRC
##run_ba(data = data, methoda = "mel_n", methodb = "mel_o")
#################################################################################

### GAMS
#################################################################################
## CRC
par(mfrow=c(2,2)) 

model0 <- gam(casoc~ s(visual_n) + a4_edat + area + education_basic_final + a3_sexe +
                bmi_complementario + score_SE + FH_cancer_Prostate + fumador_ever + uvi , data=data, family=binomial("logit"))
summary(model0)
plot.Gam(model0, se = TRUE, terms = "s(visual_n)", xlab= "visual_n", 
         ylab= "Colorectal cancer risk") #se stands for standard error Bands

model0 <- gam(casoc~ s(mel_n) + a4_edat + area + education_basic_final + a3_sexe+
                bmi_complementario + score_SE + FH_cancer_Prostate + fumador_ever + uvi , data=data, family=binomial("logit"))
summary(model0)
plot.Gam(model0, se = TRUE, terms = "s(mel_n)", xlab= "mel_n", 
         ylab= "Colorectal cancer risk") #se stands for standard error Bands

model0 <- gam(casoc~ s(visual_o) + a4_edat + area + education_basic_final + a3_sexe + 
                bmi_complementario + score_SE + FH_cancer_Prostate + fumador_ever + uvi , data=data, family=binomial("logit"))
summary(model0)
plot.Gam(model0, se = TRUE, terms = "s(visual_o)", xlab= "visual_o", 
         ylab= "Colorectal cancer risk") #se stands for standard error Bands

model0 <- gam(casoc~ s(mel_o) + a4_edat + area + education_basic_final + a3_sexe +
                bmi_complementario + score_SE + FH_cancer_Prostate + fumador_ever + uvi , data=data, family=binomial("logit"))
summary(model0)
plot.Gam(model0, se = TRUE, terms = "s(mel_o)", xlab= "mel_o", 
         ylab= "Colorectal cancer risk") #se stands for standard error Bands

model0 <- gam(casoc~ s(dmsp11) + a4_edat + area + education_basic_final + a3_sexe+
                bmi_complementario + score_SE + FH_cancer_Prostate + fumador_ever + uvi , data=data, family=binomial("logit"))
summary(model0)
plot.Gam(model0, se = TRUE, terms = "s(dmsp11)", xlab= "dmsp11", 
         ylab= "Colorectal cancer risk") #se stands for standard error Bands

model0 <- gam(casoc~ s(viirs) + a4_edat + area + education_basic_final + a3_sexe+
                bmi_complementario + score_SE + FH_cancer_Prostate + fumador_ever + uvi , data=data, family=binomial("logit"))
summary(model0)
plot.Gam(model0, se = TRUE, terms = "s(viirs)", xlab= "viirs", 
         ylab= "Colorectal cancer risk, OR (95% CI)") #se stands for standard error Bands

#################################################################################

### GLM
#################################################################################
run_model_continuous <- function (outcome, exposure, data, covars, output_dir) {
  
  adjust <- c(exposure, covars)
  f <- as.formula(paste0(outcome, "~", paste(adjust, collapse=" + ")))
  
  result <- glm(f,family=binomial("logit"),data=data)
  
  # Get statistics 
  OR <- exp(as.data.frame(result$coefficients)[c(exposure),])
  CI <- t(exp(as.data.frame(confint(result)[c(exposure),])))
  pval <- summary(result)$coefficients[exposure,"Pr(>|z|)" ]
  caso <- table(data[,outcome])[2]
  control <- table(data[,outcome])[1]
 
  # Get final df
  
  complete_result <- cbind(caso, control, OR, CI, pval)
  rownames(complete_result) <- exposure
  return(as.data.frame(complete_result))
  
}
run_model_categorical<- function (outcome, exposure, data, covars, output_dir) {

  print(paste0('OUR EXPOSURE IS: ', exposure))
  print(paste0('OUR COVARS ARE: ', covars))
  
  adjust <- c(exposure, covars)
  f <- as.formula(paste0(outcome, "~", paste(adjust, collapse=" + ")))
  result <- glm(f,family=binomial("logit"),data=data)
  
  # Get statistics 
  OR <- exp(as.data.frame(result$coefficients)[2,])
  CI <- t(exp(as.data.frame(confint(result)[2,])))
  pval <- summary(result)$coefficients[2,"Pr(>|z|)" ]
  ORhigh <- exp(as.data.frame(result$coefficients)[3,])
  CIhigh <- t(exp(as.data.frame(confint(result)[3,])))
  pvalhigh <- summary(result)$coefficients[3,"Pr(>|z|)" ]
  tbl <- table(data[,c(outcome,exposure)])
  #caso_r <- tbl[2,1] # I will need to fix this
  #control_r <- tbl[1,1] # I will need to fix this
  caso <- tbl[2,2]
  control <- tbl[1,2]
  caso_high <- tbl[2,3]
  ctrl_high  <- tbl[1,3]

  # Get final df
  complete_result_rows <- cbind(caso, control, OR, CI, pval)
  rownames(complete_result_rows) <- names(result$coeff[2])
  complete_result_rows_2 <- cbind(caso_high, ctrl_high, ORhigh, CIhigh, pvalhigh)
  rownames(complete_result_rows_2) <- names(result$coeff[3])
  complete_result <- rbind(complete_result_rows, complete_result_rows_2)
  return(as.data.frame(complete_result))
  
}

# Continuous exposure variables 
exposures_cont <- c("visual_n", "mel_n", 
                    "dmsp11", "viirs")
# Categorical exposure variables 
exposures_cat <- c("visual_n_tert", "mel_n_tert", 
                   "dmsp_11_tert", "viirs_tert")

# COLORECTAL
# datasets

#no information on chronotype
data$casoc <- as.factor(data$casoc)
colorect <- data[!is.na(data$casoc),]
colorect.2 <- colorect[!is.na(colorect$a4_edat) & !is.na(colorect$area) & !is.na(colorect$education_basic_final) & !is.na(colorect$a3_sexe) & !is.na(colorect$FH_cancer_Colon), ] 
colorect.3 <- colorect.2[!is.na(colorect.2$bmi_complementario) & !is.na(colorect.2$fumador_ever) & !is.na(colorect.2$uvi) & !is.na(colorect.2$score_WCRF_simple), ] 
colorect.4 <- colorect.3[!is.na(colorect.3$night_ever_recoded), ] 

summary(colorect$night_ever_recoded)
summary(colorect$score_SE)
#covariates
covars.1 <- c("a3_sexe", "a4_edat", "area",  "education_basic_final", "FH_cancer_Colon" )
covars.2 <- c(covars.1, "bmi_complementario", "fumador_ever", "uvi", "score_WCRF_simple")
covars.3 <- c(covars.2, "night_ever_recoded")
covars.4 <- c(covars.3, "ndvi_300", "pm25", "no2")

#models with continuous exposures
list_result.1 <- lapply(outcome = "casoc", exposures_cont, run_model_continuous, data=colorect.2, covars = covars.1)
list_result.2 <- lapply(outcome = "casoc", exposures_cont, run_model_continuous, data=colorect.3, covars = covars.2)
list_result.3 <- lapply(outcome = "casoc", exposures_cont, run_model_continuous, data=colorect.4, covars = covars.3)
list_result.4 <- lapply(outcome = "casoc", exposures_cont, run_model_continuous, data=colorect.4, covars = covars.4)


df_c1 <- bind_rows(list_result.1)
df_c2 <- bind_rows(list_result.2)
df_c3 <- bind_rows(list_result.3)
df_c4 <- bind_rows(list_result.4)
df_c <- cbind(df_c1,df_c2,df_c3, df_c4)
df_c <- round(df_c,2)

#models with categorical exposues
list_result.1 <- lapply(outcome = "casoc", exposures_cat, run_model_categorical, data=colorect.2, covars = covars.1)
list_result.2 <- lapply(outcome = "casoc", exposures_cat, run_model_categorical, data=colorect.3, covars = covars.2)
list_result.3 <- lapply(outcome = "casoc", exposures_cat, run_model_categorical, data=colorect.4, covars = covars.3)
list_result.4 <- lapply(outcome = "casoc", exposures_cat, run_model_categorical, data=colorect.4, covars = covars.4)

df_c1 <- bind_rows(list_result.1)
df_c2 <- bind_rows(list_result.2)
df_c3 <- bind_rows(list_result.3)
df_c4 <- bind_rows(list_result.4)
df_c_cat <- cbind(df_c1,df_c2,df_c3, df_c4)
df_c_cat <- round(df_c_cat,2)

####
df <- rbind(df_c, df_c_cat)
View(df)
write.csv(df,"//fs01.isglobal.lan/mcc-spain-ext/04_Analysis/Datos ALAN/04_Analysis/Results/GLM/Results_glm_withnights_24112023.csv")

summary(colorect$controls)
summary(colorect.2$controls)
summary(colorect.3$controls)
summary(colorect.4$controls)

