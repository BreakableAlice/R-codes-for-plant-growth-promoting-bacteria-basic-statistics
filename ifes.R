# Phosphate solubilizing bacteria greenhouse assay ----

# The bio assay consisted in five bacterial strain treatments, one consortium 
# containing all bacteria, one positive control containing the bio fertilizer 
# Tango and one negative control containing just bare culture medium. Following
# 15 days separated inoculations, plants were taken until tomato production 
# stages. This code accomplishes ANOVA for all the variants reviewed.
# Created by Luckie-Duque C. Alicia 2023-01-04 

## Settings ----
#### Loading packages ----
rm(list = ls())
library(readr)
library(agricolae)
library(ggplot2)
library(multcompView)
library(dplyr)
library(RColorBrewer)

#### Create data frames for appending ----

sumTab <- data.frame()
compTab <- data.frame()

#### File reading from directory containing all downloaded csv files from project repository. ----

#setwd("C:/Users/lime_/Documents/work/files") One can specify the directory here
#so the code file is saved on a higher directory independent from the one 
#containing all csv files. In this rmd file is specified upstream because 
#it requires off-chunk direction. 
file_list <- list.files(pattern = "*.csv")

file_list


## CODE ----

#### Main loop to iterate through files. ----
j=1
for (j in 1:length(file_list)) {
  filename = file_list[j]
  df <- read.csv(filename)
  
  # "Sampling" is to get id from each name file.
  sampling <- gsub("_clean|.csv", "", filename)
  print(sampling)
  
  filename
  df
  
  cat("Currently working with: ", filename, "\n")
  
  # Change variables to correct data type.
  df$bed <- as.factor(df$bed)
  
  #### Set elements for model. ----
  Y = df[,grep("bed|treatment|replicate|condition|X", colnames(df), invert = T)]
  x = df$treatment
  b = df$bed
  
  # Set data frame for output display.
  anovaTab = data.frame(matrix(0,nrow=ncol(Y), ncol=3))
  rownames(anovaTab) = colnames(Y)
  colnames(anovaTab) = c("anovaA", "meanShapiro", "Sample")
  
  #### Secondary loop to iterate through each variable (or column) in each file. ----
  i=1
  for(i in 1:ncol(Y)){
    tag = colnames(Y)[i]
    plot_tag = tag
    y = Y[,i]
    # Anova with previous set model.
    anovaA = aov(y~x+b, data = df)
    value = summary(anovaA)[[1]][1,5]
    # If to only run Shapiro and HSD on significant variables.
    if (value <= 0.05) {
      cat(tag, "", "was significant" ,"\n")
      anovaTab[i,"anovaA"] = summary(anovaA)[[1]][1,5]
      anovaTab[i, "meanShapiro"] = mean(tapply(y,x,function(x)shapiro.test(x)$p.value))  
      anovaTab[i, "Sample"] = sampling
      comparaciones = HSD.test(anovaA, "x", group=TRUE)
    
      # Group retrieving to include on plot.
      groupsf = cbind(comparaciones$groups, comparaciones$means[rownames(comparaciones$groups),-1])
      groupsf$tratamientos = rownames(groupsf)
      
      cat("\n")
      print(groupsf[1:2])
      
      # Automatize title with time, variable, and P-value.
      titling = paste("For", sampling, tag, "was significant", "P =", value)
      
      # Plot settings. ----
      mycolors = brewer.pal(9, "Set1")
      figure <- ggplot() + 
        geom_boxplot(df, mapping = aes(x = treatment, y = y, col = treatment, fill = treatment), alpha = 0.6) +
        scale_color_manual(values = mycolors) +
        scale_fill_manual(values = mycolors) +
        geom_text(data = groupsf, aes(label = groups, x = tratamientos, y = y)) +
        labs(x = "Treatment", y = plot_tag, title = titling) +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
      
      print(figure)
        
      } else {
      # Prints a "was NOT significant" warning.   
      cat(tag, "", "was NOT significant" ,"\n")
    }
    
  }

  # Prints summary tables including singificant variables and scores. 
  well = anovaTab[anovaTab$anovaA <= 0.05 & anovaTab$anovaA != 0,]
  sumTab = rbind(well, sumTab) 

  cat("\n")
  
  
  cat("All done: ", filename, "\n")
  cat("\n", "\n")
  
}

