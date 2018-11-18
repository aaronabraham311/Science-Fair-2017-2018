library(tidyverse)
library(caret)
library(gridExtra)

# Reading data
df <- read.csv("data/brain_volume_mlData.csv")
brain_vol_data <- df[,c("Subject.ID", "lh_MeanThickness_thickness.mm", "rh_MeanThickness_thickness.mm", "BrainSegVolNotVent.mm3", "estimated.total.intracranial.volume.mm3")]

# Normalizing brain_volume:
brain_vol_dataNorm <- as.data.frame(scale(brain_vol_data))

# Finding best cluster value
bss <- numeric()
wss <- numeric()

for(i in 1:10){
  
  # For each k, calculate betweenss and tot.withinss
  bss[i] <- kmeans(brain_vol_dataNorm, centers=i)$betweenss
  wss[i] <- kmeans(brain_vol_dataNorm, centers=i)$tot.withinss
}

# Between sum of squares (maximize)
p3 <- qplot(1:10, bss, geom=c("point", "line"), 
            xlab="Number of clusters", ylab="Between-cluster sum of squares") +
  scale_x_continuous(breaks=seq(0, 10, 1))

# Within sum of squares (minimize)
p4 <- qplot(1:10, wss, geom=c("point", "line"),
            xlab="Number of clusters", ylab="Total within-cluster sum of squares") +
  scale_x_continuous(breaks=seq(0, 10, 1))

grid.arrange(p3, p4, ncol=2)

# k = 5 looks like best cluster value
brain_vol_cluster <- kmeans(brain_vol_dataNorm, centers = 5) #Read in

# Visualizing clusters
ggplot(brain_vol_dataNorm, aes(BrainSegVolNotVent.mm3, estimated.total.intracranial.volume.mm3, 
                               color = factor(brain_vol_cluster$cluster))) + geom_point() + ggtitle("Clusters") + labs(x = "BrainSegVolNotVent", y = "Estimated Total Intracranial Volume", colour = "Clusters")

df$Name <- factor(df$Name)
kClusters <- factor(brain_vol_cluster$cluster)
summarizedResults <- table(df$Name, kClusters)


ggplot(df, aes(kClusters, BrainSegVolNotVent.mm3, fill = kClusters)) + 
  geom_bar(stat = "identity") + ggtitle("Total Brain Volume Per Cluster") + labs(x = "Cluster", y = "Total Brain Volume", colour = "Cluster")
ggplot(df, aes(kClusters, lh_MeanThickness_thickness.mm, fill = kClusters)) + 
  geom_bar(stat = "identity") + ggtitle("Left Hemisphere Mean Cortical Thickness Per Cluster") + labs(x = "Cluster", y = "Left Hemisphere Mean Cortical Thickness", color = "Cluster")
ggplot(df, aes(kClusters, rh_MeanThickness_thickness.mm, fill = kClusters)) + 
  geom_bar(stat = "identity") + ggtitle("Right Hemisphere Mean Cortical Thickness Per Cluster") + labs(x = "Cluster", y = "Right Hemisphere Mean Cortical Thickness", color = "Cluster")
ggplot(df, aes(kClusters, estimated.total.intracranial.volume.mm3, fill = kClusters)) + 
  geom_bar(stat = "identity") + ggtitle("Total Skull Volume Per Cluster") + labs(x = "Cluster", y = "Total Skull Volume", color = "Cluster")
