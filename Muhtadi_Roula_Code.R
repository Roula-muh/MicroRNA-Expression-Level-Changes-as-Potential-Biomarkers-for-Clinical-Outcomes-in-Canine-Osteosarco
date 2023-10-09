#install.packages("readxl")
#install.packages("tidyverse")
#install.packages("dplyr")
#install.packages("gplots")
#install.packages("survival")
#install.packages("glmnet")
#install.packages("ggplot2")
#install.packages("survminer")

#if (!require("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")
#BiocManager::install("ddCt")

library(readxl)
library("tidyverse")
library(dplyr)
library(gplots)
library(survival)
library(glmnet)
library(ggplot2)
library(survminer)


# Data pre-processing  ----------------------------------------------------


#Import data into a dataframe
Osteosarcoma_Data <- read_excel("Osteosarcoma_Data.xlsx", sheet = "Original Population", col_names = FALSE)

#Transpose the dataframe and convert matrix produced to a dataframe
transposed_df <- t(Osteosarcoma_Data)
DogData <- as.data.frame(transposed_df)

#Merge first 2 columns to have unique headers
DogData <- DogData %>%
  unite(DogName, V1, V2, sep = "_")

#Make first column into rows header
rownames(DogData) <- DogData[ , 1]
DogData <- DogData[, -1]

#Make first row into columns header
colnames(DogData) <- DogData[1, ]
DogData <- DogData[-1, ]

#Transpose data and convert to a dataframe for downstream analysis
transposedDogData <- t(DogData)
transposedDogData <- as.data.frame(transposedDogData)

#remove control miRNAs and miRNAs with 0s 
transposedDogData <- transposedDogData[-c(54:57), ]
transposedDogData <- transposedDogData[rowSums(transposedDogData == 0) == 0, ]

#Lable post data as group 1 and pre data as group 2
group_labels <- ifelse(grepl("^Post_", names(transposedDogData)), 1, ifelse(grepl("^pre_", names(transposedDogData)), 0, 2))
transposedDogData <- rbind(transposedDogData, Group = group_labels)

#convert dataframe into numeric for downstream analysis
transposedDogData <- mutate_all(transposedDogData, as.numeric)


#Import survival data 
Survival_data <- read_excel("Osteosarcoma_Data.xlsx", sheet = "Survival", col_names = FALSE)

Covariates <- read_excel("Osteosarcoma_Data.xlsx", sheet = "Dog Info OVC", col_names = FALSE)

#Clean Survival dataframe
#Make first column into rows header
Survival_data <- as.data.frame(Survival_data)
rownames(Survival_data) <- Survival_data[ ,1]
Survival_data <- Survival_data[, -1]


#Make first row into columns header
colnames(Survival_data) <- Survival_data[1, ]
Survival_data <- Survival_data[-1, ]


#Clean covariates dataframe
#Make first column into rows header
Covariates <- as.data.frame(Covariates)
rownames(Covariates) <- Covariates[ ,1]
Covariates <- Covariates[, -1]

#Make first row into columns header
colnames(Covariates) <- Covariates[1, ]
Covariates <- Covariates[-1, ]

#Clean up
rm(group_labels,transposed_df, Osteosarcoma_Data)

# NormFinder --------------------------------------------------------------

#Norm Finder on all dataset

#Download the script as a function
source("r.NormOldStab5.txt")

#Export DogData into a text file since NormFinder can only work on text files
write.table(transposedDogData, file = "DogData.txt", sep = "\t", quote = FALSE, row.names = TRUE)

#Apply NormFinder function on data
Result=Normfinder("DogData.txt")

#Explore result from applying NormFinder to data and find the most stable genes to be used as reference genes for normalization
Result$Ordered
Result$UnOrdered
Result$PairOfGenes


# ddCT  -------------------------------------------------------------------

#Pre-processing 
#Remove Group row because no longer required in analysis
transposedDogData <- transposedDogData[-45, ]

#extract reference genes into their own data frame from the pre-amputation dataset
reference_genes <-  transposedDogData[c("cfa.miR.140", "hsa.miR.22.3p","hsa.miR.126.5p","cfa.miR.23b"),]


# Calculate the row-wise average and add the average row to the dataframe
reference_genes_average <- colMeans(reference_genes)
reference_genes <- rbind(reference_genes, Average= reference_genes_average)
reference_genes_average <- reference_genes[-c(1:4), ]


# Create an empty data frame to store the dCt values
dCt_values <- data.frame(matrix(ncol = ncol(transposedDogData), nrow = nrow(transposedDogData)))
colnames(dCt_values) <- colnames(transposedDogData)
rownames(dCt_values) <- rownames(transposedDogData)


# Calculate dCt values for each sample and miRNA
for (sample in colnames(transposedDogData)) {
  dCt_values[, sample] <- transposedDogData[, sample] - reference_genes_average[, sample]
}

# Extract the dog names from the column names
dog_names <- unique(sub("(Post_|Pre_)", "", colnames(dCt_values)))

# Calculate delta delta Ct (ddCt) by subtracting the pre-amputation dCt values from the post-amputation dCt values
ddCt_data <- data.frame(row.names = row.names(dCt_values))
for (dog in dog_names) {
  post_col <- paste0("Post_", dog)
  pre_col <- paste0("Pre_", dog)
  ddCt_data[[paste(dog)]] <- dCt_values[[post_col]] - dCt_values[[pre_col]]
}

rm(dCt_values, reference_genes, reference_genes_average, DogData, transposedDogData, dog, dog_names, post_col, pre_col, sample, Normfinder, Result)


# Visualizing Data --------------------------------------------------------

#Pre-process
survival_days <- as.data.frame(Survival_data[ , "Overall", drop = FALSE])

#To order based on ascending survival days

#Switch values to numeric
survival_days$Overall <- as.numeric(survival_days$Overall)

#Order in ascending order based on survival days 
survival_days <- as.data.frame(t(survival_days[order(survival_days$Overall), , drop = FALSE]))

#Extract column names
column_order <- colnames(survival_days)

#Order based on the survival_days dataframe 
sorted_ddCt_data <- ddCt_data[, column_order]


#Create the heatmap for visualization
colnames(sorted_ddCt_data) <- NULL

heatmap.2(
  as.matrix(sorted_ddCt_data),
  Rowv = TRUE,
  Colv = FALSE,
  col = colorRampPalette(c("blue", "white", "red"))(100),
  scale = "none",
  key = FALSE,  # Turn off the default legend
  trace = "none",
  dendrogram = "none",
  xlab = "Dogs",
  cexRow = 0.5,
  )

legend(
  x = "topright",  # Position the legend at the top-right corner
  legend = c("Low Expression", "High Expression"),
  fill = c("blue", "red"),
  title = "Expression Levels",
  cex = 0.8,
  inset = c(0, -0.2),  # Adjust the legend position (move it up and to the left by increasing the second value)
  xpd = TRUE  # Plot the legend outside the plotting area
)

#Clean up
rm(sorted_ddCt_data, survival_days, column_order)


# Prepare Data for analysis -----------------------------------------------


#Clean ddCT data 
ddCt_data <- t(ddCt_data)
ddCt_data <- as.data.frame(ddCt_data)
ddCt_data <- mutate_all(ddCt_data, as.numeric)


# Pre-processing before Main analysis ---------------------------------------------------------

#Extract covariates that will be used in downstream analysis
Covariates <- data.frame(Covariates$Age, Covariates$Weight)
colnames(Covariates) <- c("Age", "Weight")
Covariates <- mutate_all(Covariates, as.numeric)


#Overall survival data
OverallSurv <- Survival_data[, c("Overall","Censor (OS)"), drop = FALSE ]
colnames(OverallSurv) <- c("Time", "Censor")
OverallSurv <- mutate_all(OverallSurv, as.numeric)


#DFI survival data 
DFISurv <- Survival_data[, c("DFI","Censor (DFI)"), drop = FALSE ]
colnames(DFISurv) <- c("Time", "Censor")
DFISurv <- mutate_all(DFISurv, as.numeric)


#Create a dataset that includes age and weight as covariates in addition to miRNA expression levels
With_covariates <- cbind(ddCt_data, Covariates)



# MAIN ANALYSIS  ----------------------------------------------------------



# Overall Survival Data Exploration  -------------------------------------------------------

#Count how many events total, how many events are censored and how many aren't censored
length(OverallSurv$Time)
event_counts <- sum(OverallSurv$Censor == 1, na.rm = TRUE)
censored_counts <- sum(OverallSurv$Censor == 0, na.rm = TRUE)

#Percentage of events censored is 25% of data (loss to follow-up)
sum(OverallSurv$Censor == 0)/length(OverallSurv$Time)*100


# Create a bar plot for events vs censored 
barplot(c(event_counts, censored_counts), 
        names.arg = c("Event", "Censored"), 
        xlab = "Status", 
        ylab = "Count", 
        main = "Distribution of Censoring and Events")



# DFI Survival Data Exploration -------------------------------------------


#Count how many events total, how many events are censored and how many aren't censored
length(DFISurv$Time)
event_counts <- sum(DFISurv$Censor == 1, na.rm = TRUE)
censored_counts <- sum(DFISurv$Censor == 0, na.rm = TRUE)

#Percentage of events censored is 25% of data (loss to follow-up)
sum(DFISurv$Censor == 0)/length(DFISurv$Time)*100

# Create a bar plot for events vs censored 
barplot(c(event_counts, censored_counts), 
        names.arg = c("Event", "Censored"), 
        xlab = "Status", 
        ylab = "Count", 
        main = "Distribution of Censoring and Events")

rm(event_counts, censored_counts)


# Survival Model Pre-processing  --------------------------------------------------------

#Create training and test sets for downstream use

dogs_to_leave_out <- c("Smokey", "Buffy", "Chubbs", "Ellie", "Letta")

test_indices <- which(rownames(ddCt_data) %in% dogs_to_leave_out)


# Separate the test set
test_with_covariates <- With_covariates[test_indices, ]
test_OS <- OverallSurv[test_indices,]
test_without_covariates <- ddCt_data[test_indices,]
test_DFI <- DFISurv[test_indices,]

# Separate the training set
training_with_covariates <- With_covariates[-test_indices, ]
training_OS <- OverallSurv[-test_indices,]
training_without_covariates <- ddCt_data[-test_indices,]
training_DFI <- DFISurv[-test_indices,]

#Create survival objects
test_OS_surv <- Surv(time = test_OS$Time, event = test_OS$Censor)

test_DFI_surv <- Surv(time = test_DFI$Time, event = test_DFI$Censor)

training_OS_surv <- Surv(time = training_OS$Time, event = training_OS$Censor)

training_DFI_surv <- Surv(time = training_DFI$Time, event = training_DFI$Censor)


# Model 1: Overall survival with covariates -------------------------------

#cross-validate
set.seed(123)

#Fit the cross validated model 
cvfit1 <- cv.glmnet(as.matrix(training_with_covariates), training_OS_surv, family = "cox", alpha = 1, nfolds = 4)


cvfit1$lambda.min
cvfit1$glmnet.fit

#Plot the model to visualize
plot(cvfit1)


#Fit the model using ideal lambda and find coefficients then exponentiation to interpret using hazard ratio 
final_model1 <- glmnet(as.matrix(training_with_covariates), training_OS_surv, family = "cox", alpha = 1, lambda = 0.31890)


#Plot the general survival curve
fit1 <- survival::survfit(final_model1, s = 0.31890, x = as.matrix(training_with_covariates), y = training_OS_surv)


plot(fit1, xlab = "Days of Survival", ylab = "Probability of Survival", lwd = 2)

# Identify censored events
censored_indices <- which(training_OS$Censor == 0)

# Add ticks on censored events
points(fit1$time[censored_indices], fit1$surv[censored_indices], pch = "|", col = "red")

#Add median survival line
median_survival_time <- median(training_OS$Time)
abline(v = median_survival_time, col = "blue", lty = 2)

#Label the ticks
legend("topright", legend = c("Loss to Follow Up", "Median Survival Time"), pch = c("|", NA), col = c("red", "blue"), lty = c(NA, 2), cex = 1, box.lwd = 0, bty = "n")


#Find coefficients to get significant predictors
model1 <- coef(cvfit1, s = 0.31890)
ifelse(model1 != 0, exp(model1), 0) 




# Model 2:  Overall survival without covariates -----------------------------------------------------------------

#Fit the cross validated model 
cvfit2 <- cv.glmnet(as.matrix(training_without_covariates), training_OS_surv, family = "cox", alpha = 1, nfolds = 4)

cvfit2$lambda.min
cvfit2$glmnet.fit

#Plot the model to visualize
plot(cvfit2)

#Fit the model using ideal lambda and find coefficients then exponentiation to interpret using hazard ratio 
final_model2 <- glmnet(as.matrix(training_without_covariates), training_OS_surv, family = "cox", alpha = 1, lambda = 0.31890)

#Plot the general survival curve
fit2 <- survival::survfit(final_model2, s = 0.31890, x = as.matrix(training_without_covariates), y = training_OS_surv)


#Find coefficients to get significant predictors
model2 <- coef(cvfit2, s = 0.31890)
ifelse(model2 != 0, exp(model2), 0) 



#Model 3: DFI with covariates  -----------------------------------------------------------

#Fit the cross validated model 
cvfit3 <- cv.glmnet(as.matrix(training_with_covariates), training_DFI_surv, family = "cox", alpha = 1, nfolds = 4)

cvfit3$lambda.min
cvfit3$glmnet.fit

#Plot 
plot(cvfit3)

#Get predictors using coefficients 
coef(cvfit3, s = 0.33540)

#Fit the model using ideal lambda and find coefficients then exponentiation to interpret using hazard ratio 
final_model3 <- glmnet(as.matrix(training_with_covariates), training_DFI_surv, family = "cox", alpha = 1, lambda = 0.33540)

#Plot the general survival curve
fit3 <- survival::survfit(final_model3, s = 0.33540, x = as.matrix(training_with_covariates), y = training_DFI_surv)


plot(fit3, xlab = "Days of Survival", ylab = "Probability of Survival", lwd = 2)

# Identify censored events
censored_indices <- which(training_DFI$Censor == 0)

# Add ticks on censored events
points(fit3$time[censored_indices], fit3$surv[censored_indices], pch = "|", col = "red")

#Add median survival line
median_survival_time <- median(training_DFI$Time)
abline(v = median_survival_time, col = "blue", lty = 2)

#Label the ticks
legend("topright", legend = c("Loss to Follow Up", "Median Survival Time"), pch = c("|", NA), col = c("red", "blue"), lty = c(NA, 2), cex = 1, box.lwd = 0, bty = "n")

#Exponentiate the non-zero coefficients for interpretation
model3 <- coef(cvfit3, s = 0.33540)
ifelse(model3 != 0, exp(model3), 0) 


# #Model 4: DFI without covariates ----------------------------------------

#Fit the cross validated model 
cvfit4 <- cv.glmnet(as.matrix(training_without_covariates), training_DFI_surv, family = "cox", alpha = 1, nfolds = 4)

cvfit4$lambda.min
cvfit4$glmnet.fit

#Visualize lambda
plot(cvfit4)

#Fit the model using ideal lambda and find coefficients then exponentiation to interpret using hazard ratio 
final_model4 <- glmnet(as.matrix(training_without_covariates), training_DFI_surv, family = "cox", alpha = 1, lambda = 0.33540)

#Get predictors using coefficients 
coef(cvfit4, s = 0.33540)


#Exponentiate the non-zero coefficients for interpretation
model4 <- coef(cvfit4, s = 0.33540)
ifelse(model4 != 0, exp(model4), 0) 


#remove dataframes not required
#rm(OverallSurv, Survival_data, With_covariates, dogs_to_leave_out, test_indices, Covariates, ddCt_data, DFISurv)

#rm(censored_indices, median_survival_time, order_dogs, p_value, r_squared, test_indices, cox_model3, cox_model4, fit, fit1, fit2, fit3, fit4, model1, model2, model3, model4, plot_HR_DFI, plot_HR_OS, predicted_model1, predicted_model2, prediction_model1, prediction_model1_test, prediction_model2, prediction_model2_test, predictor_model1, sorted_predicted_model1, sorted_predicted_model2, sorted_survival_days, sorted_survival_days_DFI, Covariates, final_model1, final_model2, final_model3, final_model3, Survival_data, survival_days_DFI, survival_days_OS, survival_days_test, survival_days_test_DFI, survival_days_training, survival_days_training_DFI, test_DFI, test_OS, test_with_covariates, test_without_covariates, training_DFI, training_OS, With_covariates, test_DFI_surv, test_OS_surv, DFISurv, final_model4, training_without_covariates, training_with_covariates, training_DFI_surv, training_OS_surv, dogs_to_leave_out)



# Prediction --------------------------------------------------------------

prediction_model1_test <- predict(final_model1, newx = as.matrix(test_with_covariates), s = 0.31890, type = "response")

prediction_model2_test <- predict(final_model3, newx = as.matrix(test_with_covariates), s = 0.33540, type = "response")


# Refitting the model using significant predictors ------------------------

#Use univariate unpenalized cox regression to refit the model 
survival_object <- Surv(time = DFISurv$Time, event = DFISurv$Censor)

cox_data <- data.frame(survival_object, miRNA_1 = ddCt_data$hsa.miR.20a.5p, miRNA_2 = ddCt_data$hsa.miR.451a, miRNA_3 = ddCt_data$hsa.miR.93.5p)


#Pvalue of miR-451a is 0.07887196 (not significant)
cox_model_miRNA2 <- coxph(survival_object ~ miRNA_2, data = cox_data)
(summary(cox_model_miRNA2))$coefficients[, "Pr(>|z|)"]

#Pvalue of miR-93-5p is 0.09278443 (not significant)
cox_model_miRNA3 <- coxph(survival_object ~ miRNA_3, data = cox_data)
(summary(cox_model_miRNA3))$coefficients[, "Pr(>|z|)"]
#Fit models and get p-values then plot 

#Pvalue of miR-20a-5p is 0.027* 
cox_model_miRNA1 <- coxph(survival_object ~ miRNA_1, data = cox_data)
(summary(cox_model_miRNA1))$coefficients[, "Pr(>|z|)"]



# Create the plot to compare overall survival curve to hsa.miR.20a --------


fit_specific <- survfit(cox_model_miRNA1)


fit_overall <- survfit(survival_object ~ 1, data = ddCt_data)

plot_general_to_miRNa <- data.frame(time = fit_specific$time,
                                    survival_specific = fit_specific$surv,
                                    survival_overall = fit_overall$surv)

ggplot(plot_general_to_miRNa, aes(x = time)) +
  geom_step(aes(y = survival_specific, color = "Specific miRNA"), size = 1) +
  geom_step(aes(y = survival_overall, color = "Overall cases"), linetype = "dashed", size = 1) +
  labs(title = "Kaplan-Meier Curves for Specific miRNA and Overall Cases",
       y = "Survival Probability", x = "Time",
       color = "Group") +
  theme_minimal() +
  theme(legend.position = "bottom")


# Create a plot to compare higher and lower expression levels  ------------

# Extract the specific miRNA expression
miRNA_expression <- ddCt_data$hsa.miR.20a.5p

# Calculate the median expression
median_expression <- median(miRNA_expression)

# Create a binary variable indicating high or low expression based on median
expression_group <- ifelse(miRNA_expression >= median_expression, "High", "Low")

data <- data.frame(survival_object, expression_group)

# Extract miRNA data for high expression group
miRNA_expression_high <- as.data.frame(miRNA_expression[data$expression_group == "High"])


#extract data_high datapoints to create a survival object 
data_high_survival <- DFISurv[c(2, 3, 4, 5, 7, 9, 10, 13, 15, 16, 20, 21, 22, 25), ]

data_high_survival <- Surv(time = data_high_survival$Time, event = data_high_survival$Censor)

# Fit Kaplan-Meier curve for high expression
fit_high <- survfit(data_high_survival ~ 1, data = miRNA_expression_high)

# Filter data for low expression group
# Extract miRNA data for low expression group

miRNA_expression_low <- as.data.frame(miRNA_expression[data$expression_group == "Low"])

#extract data_high datapoints to create a survival object 
data_low_survival <- DFISurv[c(1, 6, 8, 11, 12, 14, 17, 18, 19, 23, 24, 26, 27), ]

data_low_survival <- Surv(time = data_low_survival$Time, event = data_low_survival$Censor)

# Fit Kaplan-Meier curve for high expression
fit_low <- survfit(data_low_survival ~ 1, data = miRNA_expression_low)


# Create the plot data
plot_data <- data.frame(
  time = c(fit_specific$time, fit_high$time, fit_low$time),
  survival_specific = c(fit_specific$surv, rep(NA, length(fit_high$time)), rep(NA, length(fit_low$time))),
  survival_high = c(rep(NA, length(fit_specific$time)), fit_high$surv, rep(NA, length(fit_low$time))),
  survival_low = c(rep(NA, length(fit_specific$time)), rep(NA, length(fit_high$time)), fit_low$surv)
)

# Create the plot
ggplot(plot_data, aes(x = time)) +
  geom_step(aes(y = survival_specific, color = "General Curve"), size = 1) +
  geom_step(aes(y = survival_high, color = "Higher expression"), linetype = "dashed", size = 1) +
  geom_step(aes(y = survival_low, color = "Lower expression"), linetype = "dashed", size = 1) +
  labs(y = "Survival Probability", x = "Time",
       color = "Group") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  scale_color_manual(values = c("General Curve" = "black",
                                "Higher expression" = "magenta",
                                "Lower expression" = "blue"))


