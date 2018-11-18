library(tidyverse)

nd <- read.csv("new-data.csv")

adData <- nd %>% filter(group_recruit == "AD") %>% select ("Arg.Average", "Cit.Average", "Pro.Average", "Orn.Average", "memory_z", "exec_z", "speed_z")
controlData <- nd %>% filter(group_recruit == "Control") %>% select ("Arg.Average", "Cit.Average", "Pro.Average", "Orn.Average", "memory_z", "exec_z", "speed_z")

# Note that there are much more control than AD

# AD Mean and SD
arginineAD <- adData%>% select ("Arg.Average") %>% unlist
citrullineAD <- adData %>% select ("Cit.Average") %>% unlist
ornithineAD <- adData %>% select ("Orn.Average") %>% unlist
prolineAD <- adData %>% select ("Pro.Average") %>% unlist

argMeanAD <- mean(arginineAD)
argSDAD <- sd(arginineAD)

citrullineMeanAD <- mean(citrullineAD)
citrullineSDAD <- sd(citrullineAD)

ornithineMeanAD <- mean(ornithineAD)
ornithineSDAD <- sd(ornithineAD)

prolineMeanAD <- mean(prolineAD)
prolineSDAD <- sd(prolineAD)

#Control Mean and SD
arginineControl <- controlData%>% select ("Arg.Average") %>% unlist
citrullineControl <- controlData %>% select ("Cit.Average") %>% unlist
ornithineControl <- controlData %>% select ("Orn.Average") %>% unlist
prolineControl <- controlData %>% select ("Pro.Average") %>% unlist

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

# f-tests to find variances
var.test(arginineAD, arginineControl, alternative = "two.sided") # p = 0.1434 with CI (0.1907, 1.299)
var.test(citrullineAD, citrullineControl, alternative = "two.sided") # p = 0.509 with CI (0.2861158, 1.9486019)
var.test(ornithineAD, ornithineControl, alternative = "two.sided") # p = 0.9305 with CI (0.3783057, 2.5764641)
var.test(prolineAD, prolineControl, alternative = "two.sided") # p = 0.2237 with CI (0.7030441 4.7881059)

# t-tests based on variances
pArg <- t.test(arginineAD, arginineControl, var.equal = TRUE) #p = 0.00228 with CI (-0.38195215, -0.09014024)
pOrn <- t.test(ornithineAD, ornithineControl, var.equal = TRUE)# p = 0.9304 with CI(-0.3649159,  0.3980137)
pCit <- t.test(citrullineAD, citrullineControl, var.equal = TRUE) # p = 0.2375 with CI (-0.29820734, 0.07627799)
pPro <- t.test(prolineAD, prolineControl, var.equal = TRUE) # p = 0.5151 with CI ( -0.300165, 0.153127)
