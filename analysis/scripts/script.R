library(tidyverse)
library(reshape2)
library(fmsb)

data <- read.csv("data/arginine_data.csv")

adData <- data %>% filter(Name == "AD") %>% select ("Name", "Arginine", "Citrulline", "Ornithine", "Proline")
caaData <- data %>% filter(Name == "CAA") %>% select ("Arginine", "Citrulline", "Ornithine", "Proline")
controlData <- data %>% filter(Name == "CONTROL") %>% select ("Name", "Arginine", "Citrulline", "Ornithine", "Proline")

data <- rbind(adData, controlData)
data$Name <- factor(data$Name)

# AD means and standard deviations
arginineAD <- adData%>% select ("Arginine") %>% unlist
citrullineAD <- adData %>% select ("Citrulline") %>% unlist
ornithineAD <- adData %>% select ("Ornithine") %>% unlist
prolineAD <- adData %>% select ("Proline") %>% unlist

argMeanAD <- mean(arginineAD)
argSDAD <- sd(arginineAD)

citrullineMeanAD <- mean(citrullineAD)
citrullineSDAD <- sd(citrullineAD)

ornithineMeanAD <- mean(ornithineAD)
ornithineSDAD <- sd(ornithineAD)

prolineMeanAD <- mean(prolineAD)
prolineSDAD <- sd(prolineAD)

# CAA means and standard deviations
arginineCAA <- caaData%>% select ("Arginine") %>% unlist
citrullineCAA <- caaData %>% select ("Citrulline") %>% unlist
ornithineCAA <- caaData %>% select ("Ornithine") %>% unlist
prolineCAA <- caaData %>% select ("Proline") %>% unlist

argMeanCAA <- mean(arginineCAA)
argSDCAA <- sd(arginineCAA)

citrullineMeanCAA <- mean(citrullineCAA)
citrullineSDCAA <- sd(citrullineCAA)

ornithineMeanCAA <- mean(ornithineCAA)
ornithineSDCAA <- sd(ornithineCAA)

prolineMeanCAA <- mean(prolineCAA)
prolineSDCAA <- sd(prolineCAA)

# Control means and standard deviations
arginineControl <- controlData%>% select ("Arginine") %>% unlist
citrullineControl <- controlData %>% select ("Citrulline") %>% unlist
ornithineControl <- controlData %>% select ("Ornithine") %>% unlist
prolineControl <- controlData %>% select ("Proline") %>% unlist

argMeanControl <- mean(arginineControl)
argSDControl <- sd(arginineControl)

citrullineMeanControl <- mean(citrullineControl)
citrullineSDControl <- sd(citrullineControl)

ornithineMeanControl <- mean(ornithineControl)
ornithineSDControl <- sd(ornithineControl)

prolineMeanControl <- mean(prolineControl)
prolineSDControl <- sd(prolineControl)

# Combining data into table
adMeanSD <- data.frame(c("L-arginine", "Citrulline", "Ornithine", "Proline"), 
                       c(argMeanAD, citrullineMeanAD, ornithineMeanAD, prolineMeanAD),
                       c(argMeanControl, citrullineMeanControl, ornithineMeanControl, prolineMeanControl),
                       c(argSDAD, citrullineSDAD, ornithineSDAD, prolineSDAD),
                       c(argSDControl, citrullineSDControl, ornithineSDControl, prolineSDControl))
colnames(adMeanSD) <- c("Metabolite", "AD Mean","Control Mean","AD Standard Deviation", "Control Standard Deviation")

# Test normality of data https://stats.stackexchange.com/questions/101274/how-to-interpret-a-qq-plot
qqnorm(arginineAD, main = "Q-Q Plot for L-arginine")
qqline(arginineAD) # Approximately normal

qqnorm(citrullineAD, main = "Q-Q Plot for Citrulline")
qqline(citrullineAD) # Light tailed

qnorm(ornithineAD, main = "Q-Q Plot for Ornithine")
qqline(ornithineAD) # NaN are produced, not sure if this 

qnorm(prolineAD) #NaN produced
qqline(prolineAD)

#f test on variances 
var.test(arginineAD, arginineControl, alternative = "two.sided")
var.test(citrullineAD, citrullineControl, alternative = "two.sided")
var.test(ornithineAD, ornithineControl, alternative = "two.sided")
var.test(prolineAD, prolineControl, alternative = "two.sided")

# 2-sample T-test
pArg <- t.test(arginineAD, arginineControl)
pOrn <- t.test(ornithineAD, ornithineControl, var.equal = TRUE)
pCit <- t.test(citrullineAD, citrullineControl, var.equal = TRUE)
pProEqual <- t.test(prolineAD, prolineControl, var.equal = TRUE)
pProUnequal <- t.test(prolineAD, prolineControl)

# Heatmaps
cormat <- round(cor(adData),2) #Correlation matrix

get_lower_tri <- function(cormat) {
  cormat[upper.tri(cormat)] <- NA
  return (cormat)
}

get_upper_tri <- function(cormat) {
  cormat[lower.tri(cormat)] <- NA
  return (cormat)
}

upper_tri <- get_upper_tri(cormat)
melted_cormat <- melt(upper_tri, na.rm = TRUE)
ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed() +
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
  ggtitle("Correlation Heatmap of Metabolites")
  

# Density plots
density_data_ad <- data %>% select("Name", "Arginine")
density_data_ad[,"Name"] <- as.factor(density_data_ad[,"Name"])
density_data_ad[,"Condition"] <- density_data_ad[,"Name"]
ggplot(density_data_ad, aes(x = Arginine, color = Condition)) + geom_density(alpha = 0.4) +
  xlab(NULL) + ylab(NULL)

density_data_orn <- data %>% select("Name", "Ornithine")
density_data_orn[,"Name"] <- as.factor(density_data_orn[,"Name"])
density_data_orn[,"Condition"] <- density_data_orn[,"Name"]
ggplot(density_data_orn, aes(x = Ornithine, color = Condition)) + geom_density(alpha = 0.4)+
  labs(x = NULL, y = NULL)

density_data_cit <- data %>% select("Name", "Citrulline")
density_data_cit[,"Name"] <- as.factor(density_data_cit[,"Name"])
density_data_cit[,"Condition"] <- density_data_cit[,"Name"]
ggplot(density_data_cit, aes(x = Citrulline, color = Condition)) + geom_density(alpha = 0.4)+
  labs(x = NULL, y = NULL)

density_data_pro <- data %>% select("Name", "Proline")
density_data_pro[,"Name"] <- as.factor(density_data_pro[,"Name"])
density_data_pro[,"Condition"] <- density_data_pro[,"Name"]
ggplot(density_data_pro, aes(x = Proline, color = Name)) + geom_density(alpha = 0.4)+
  labs(x = NULL, y = NULL)

# Boxplot
boxplot(Arginine ~ Name, data, main = NULL, 
        xlab = "Condition", ylab = "Relative Concentration",
        col = c("red", "lightblue"))

boxplot(Ornithine ~ Name, data, main = "Ornithine Changes in AD and Control", 
        xlab = "Condition", ylab = "Relative Concentration",
        col = c("red", "lightblue"))

boxplot(Citrulline ~ Name, data, main = "Citrulline Changes in AD and Control", 
        xlab = "Condition", ylab = "Relative Concentration",
        col = c("red", "lightblue"))

boxplot(Proline ~ Name, data, main = "Proline Changes in AD and Control", 
        xlab = "Condition", ylab = "Relative Concentration",
        col = c("red", "lightblue"))

# Violin plots: combining density and boxplots
ggplot(data, aes(x = Name, y = Arginine, fill = Name)) + geom_violin(draw_quantiles = TRUE) + 
  geom_boxplot(width = 0.1) +
  labs(title = "Arginine in AD and Control")

ggplot(data, aes(x = Name, y = Ornithine, fill = Name)) + geom_violin(draw_quantiles = TRUE) + 
  geom_boxplot(width = 0.1) +
  labs(title = "Ornithine in AD and Control")

ggplot(data, aes(x = Name, y = Citrulline, fill = Name)) + geom_violin(draw_quantiles = TRUE) + 
  geom_boxplot(width = 0.1) +
  labs(title = "Citrulline in AD and Control")

ggplot(data, aes(x = Name, y = Proline, fill = Name)) + geom_violin(draw_quantiles = TRUE) + 
  geom_boxplot(width = 0.1) +
  labs(title = "Proline in AD and Control")

# Visualizing covariates
newData <- read.csv("covariate-arginine.csv")
newData <- subset(newData, select = c(-X, -Subject.ID) )
newData$Female <- as.factor(newData$Female)

newDataAD <- newData %>% filter(Name == "AD")
newDataControl <- newData %>% filter(Name == "CONTROL")

ggplot(newDataAD, aes(x = Arginine, y = Memory, color = Female)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  labs(title = "L-arginine vs Memory Scores") # No significant relation

memory_model <- lm(Memory ~ Arginine, data = newData)
summary(memory_model)
  
ggplot(newDataAD, aes(x = Arginine, y = Age, color = Female)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  labs(title = "L-arginine vs Age") # No significant relation

age_model <- lm(Age ~ Arginine, data = newData)
summary(age_model)

ggplot(newDataAD, aes(x = Arginine, y = Exec, color = Female)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  labs(title = "L-arginine vs Executive Scores") # Some clusters

exec_model <- lm (Exec ~ Arginine, data = newData)
summary(exec_model)

ggplot(newDataAD, aes(x = Arginine, y = Speed, color = Female)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  labs(title = "L-arginine vs Speed Scores") # No significant relation

speed_model <- lm(Speed ~ Arginine, data = newData)
summary(speed_model)

ggplot(newDataAD, aes(x = Female, y = Arginine, fill = Female)) +
  geom_boxplot() + 
  labs(title = "L-arginine vs Gender")

adMemoryMean <- mean(newDataAD %>% select ("Memory") %>% unlist)
adSpeedMean <- mean(newDataAD %>% select ("Speed") %>% unlist)
adExecutiveMean <- mean(newDataAD %>% select ("Exec") %>% unlist)

adMemorySD <- sd(newDataAD %>% select ("Memory") %>% unlist)
adSpeedSD <- sd(newDataAD %>% select ("Speed") %>% unlist)
adExecutiveSD <- sd(newDataAD %>% select ("Exec") %>% unlist)

controlMemoryMean <- mean(newDataControl %>% select("Memory") %>% unlist)
controlSpeedMean <- mean(newDataControl %>% select("Speed") %>% unlist)
controlExecutiveMean <- mean(newDataControl %>% select("Exec") %>% unlist)

controlMemorySD <- sd(newDataControl %>% select("Memory") %>% unlist)
controlSpeedSD <- sd(newDataControl %>% select("Speed") %>% unlist)
controlExecutiveSD <- sd(newDataControl %>% select("Exec") %>% unlist)

# Testing normality
shapiro.test(adData$Arginine) #p: 0.2343
shapiro.test(adData$Citrulline) #p: 0.04248
shapiro.test(adData$Ornithine) # p: 0.007392
shapiro.test(adData$Proline) #p: 0.005537
