####################################################
# Boxplots for the CholeskyDecomp paper
# Author: Jami Jackson Mulgrave
#
#
######################################################


#clear the workspace
rm(list = ls())


fullTable <- read.csv(file = 'CholeskyDecomp_Boxplot_Edges_p25_n25.csv', header = TRUE, sep = ',')

fullTable <- as.data.frame(fullTable)

fullTable$Dimension <- as.factor(fullTable$Dimension)


#Read in the BDGraph data and add to the table

load('BDGraph_p25_n25_tenpercent.rdata')

SP_value <- c()
SE_value <- c()
MCC_value <- c()

for (i in 1: reps) {
SP_value[i] = result_list[[i]]$SP
SE_value[i] = result_list[[i]]$SE
MCC_value[i] = result_list[[i]]$MCC

}

SE_Type = rep('Sensitivity', reps)
MCC_Type = rep('MCC', reps)

SP_Type = rep('Specificity', reps)


Edges_value = c(SP_value, SE_value, MCC_value)
Edges_Type = c(SP_Type, SE_Type, MCC_Type)

number_elements  = length(Edges_value)

Sparsity = rep('Percent', number_elements)
Method = rep('BDGraph', number_elements)
Dimension = rep(p, number_elements)

BDGraph_table = data.frame(Edges_value, Method, Edges_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, BDGraph_table)



rm(BDGraph_table)


#Read in the BDGraph data and add to the table

load('BDGraph_p25_n25_AR2.rdata')

SP_value <- c()
SE_value <- c()
MCC_value <- c()

for (i in 1: reps) {
  SP_value[i] = result_list[[i]]$SP
  SE_value[i] = result_list[[i]]$SE
  MCC_value[i] = result_list[[i]]$MCC
  
}

SE_Type = rep('Sensitivity', reps)
MCC_Type = rep('MCC', reps)

SP_Type = rep('Specificity', reps)


Edges_value = c(SP_value, SE_value, MCC_value)
Edges_Type = c(SP_Type, SE_Type, MCC_Type)

number_elements  = length(Edges_value)

Sparsity = rep('AR(2)', number_elements)
Method = rep('BDGraph', number_elements)
Dimension = rep(p, number_elements)

BDGraph_table = data.frame(Edges_value, Method, Edges_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, BDGraph_table)

rm(BDGraph_table)



#Read in the BDGraph data and add to the table

load('BDGraph_p25_n25_circle.rdata')

SP_value <- c()
SE_value <- c()
MCC_value <- c()

for (i in 1: reps) {
  SP_value[i] = result_list[[i]]$SP
  SE_value[i] = result_list[[i]]$SE
  MCC_value[i] = result_list[[i]]$MCC
  
}

SE_Type = rep('Sensitivity', reps)
MCC_Type = rep('MCC', reps)

SP_Type = rep('Specificity', reps)


Edges_value = c(SP_value, SE_value, MCC_value)
Edges_Type = c(SP_Type, SE_Type, MCC_Type)

number_elements  = length(Edges_value)

Sparsity = rep('Circle', number_elements)
Method = rep('BDGraph', number_elements)
Dimension = rep(p, number_elements)

BDGraph_table = data.frame(Edges_value, Method, Edges_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, BDGraph_table)

rm(BDGraph_table)

#Read in the Frequentist data and add to the table

load('Frequentist_p25_n25_circle.rdata')

#stars model

SP_value_stars <- c()
SE_value_stars <- c()
MCC_value_stars <- c()

for (i in 1: reps) {
  SP_value_stars[i] = result_list[[i]]$SP_freq_stars
  SE_value_stars[i] = result_list[[i]]$SE_freq_stars
  MCC_value_stars[i] = result_list[[i]]$MCC_freq_stars
  
}

SE_Type = rep('Sensitivity', reps)
MCC_Type = rep('MCC', reps)

SP_Type = rep('Specificity', reps)


Edges_value =  c(SP_value_stars, SE_value_stars, MCC_value_stars)
Edges_Type = c(SP_Type, SE_Type, MCC_Type)

number_elements  = length(Edges_value)

Sparsity = rep('Circle', number_elements)
Method = rep('StARS', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Edges_value, Method, Edges_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(Frequentist_table)


#ric model

SP_value_ric <- c()
SE_value_ric <- c()
MCC_value_ric <- c()

for (i in 1: reps) {
  SP_value_ric[i] = result_list[[i]]$SP_freq_ric
  SE_value_ric[i] = result_list[[i]]$SE_freq_ric
  MCC_value_ric[i] = result_list[[i]]$MCC_freq_ric
  
}

SE_Type = rep('Sensitivity', reps)
MCC_Type = rep('MCC', reps)

SP_Type = rep('Specificity', reps)


Edges_value =  c(SP_value_ric, SE_value_ric, MCC_value_ric)
Edges_Type = c(SP_Type, SE_Type, MCC_Type)

number_elements  = length(Edges_value)

Sparsity = rep('Circle', number_elements)
Method = rep('RIC', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Edges_value, Method, Edges_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(Frequentist_table)


#ebic model

SP_value_ebic <- c()
SE_value_ebic <- c()
MCC_value_ebic <- c()

for (i in 1: reps) {
  SP_value_ebic[i] = result_list[[i]]$SP_freq_ebic
  SE_value_ebic[i] = result_list[[i]]$SE_freq_ebic
  MCC_value_ebic[i] = result_list[[i]]$MCC_freq_ebic
  
}

SE_Type = rep('Sensitivity', reps)
MCC_Type = rep('MCC', reps)

SP_Type = rep('Specificity', reps)


Edges_value =  c(SP_value_ebic, SE_value_ebic, MCC_value_ebic)
Edges_Type = c(SP_Type, SE_Type, MCC_Type)

number_elements  = length(Edges_value)

Sparsity = rep('Circle', number_elements)
Method = rep('EBIC', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Edges_value, Method, Edges_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(Frequentist_table)


#Read in the Frequentist data and add to the table

load('Frequentist_p25_n25_AR2.rdata')

#stars model

SP_value_stars <- c()
SE_value_stars <- c()
MCC_value_stars <- c()

for (i in 1: reps) {
  SP_value_stars[i] = result_list[[i]]$SP_freq_stars
  SE_value_stars[i] = result_list[[i]]$SE_freq_stars
  MCC_value_stars[i] = result_list[[i]]$MCC_freq_stars
  
}

SE_Type = rep('Sensitivity', reps)
MCC_Type = rep('MCC', reps)

SP_Type = rep('Specificity', reps)


Edges_value =  c(SP_value_stars, SE_value_stars, MCC_value_stars)
Edges_Type = c(SP_Type, SE_Type, MCC_Type)

number_elements  = length(Edges_value)

Sparsity = rep('AR(2)', number_elements)
Method = rep('StARS', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Edges_value, Method, Edges_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(Frequentist_table)


#ric model

SP_value_ric <- c()
SE_value_ric <- c()
MCC_value_ric <- c()

for (i in 1: reps) {
  SP_value_ric[i] = result_list[[i]]$SP_freq_ric
  SE_value_ric[i] = result_list[[i]]$SE_freq_ric
  MCC_value_ric[i] = result_list[[i]]$MCC_freq_ric
  
}

SE_Type = rep('Sensitivity', reps)
MCC_Type = rep('MCC', reps)

SP_Type = rep('Specificity', reps)


Edges_value =  c(SP_value_ric, SE_value_ric, MCC_value_ric)
Edges_Type = c(SP_Type, SE_Type, MCC_Type)

number_elements  = length(Edges_value)

Sparsity = rep('AR(2)', number_elements)
Method = rep('RIC', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Edges_value, Method, Edges_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(Frequentist_table)


#ebic model

SP_value_ebic <- c()
SE_value_ebic <- c()
MCC_value_ebic <- c()

for (i in 1: reps) {
  SP_value_ebic[i] = result_list[[i]]$SP_freq_ebic
  SE_value_ebic[i] = result_list[[i]]$SE_freq_ebic
  MCC_value_ebic[i] = result_list[[i]]$MCC_freq_ebic
  
}

SE_Type = rep('Sensitivity', reps)
MCC_Type = rep('MCC', reps)

SP_Type = rep('Specificity', reps)


Edges_value =  c(SP_value_ebic, SE_value_ebic, MCC_value_ebic)
Edges_Type = c(SP_Type, SE_Type, MCC_Type)

number_elements  = length(Edges_value)

Sparsity = rep('AR(2)', number_elements)
Method = rep('EBIC', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Edges_value, Method, Edges_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(Frequentist_table)

#Read in the Frequentist data and add to the table

load('Frequentist_p25_n25_tenpercent.rdata')

#stars model

SP_value_stars <- c()
SE_value_stars <- c()
MCC_value_stars <- c()

for (i in 1: reps) {
  SP_value_stars[i] = result_list[[i]]$SP_freq_stars
  SE_value_stars[i] = result_list[[i]]$SE_freq_stars
  MCC_value_stars[i] = result_list[[i]]$MCC_freq_stars
  
}

SE_Type = rep('Sensitivity', reps)
MCC_Type = rep('MCC', reps)

SP_Type = rep('Specificity', reps)


Edges_value =  c(SP_value_stars, SE_value_stars, MCC_value_stars)
Edges_Type = c(SP_Type, SE_Type, MCC_Type)

number_elements  = length(Edges_value)

Sparsity = rep('tenpercent', number_elements)
Method = rep('StARS', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Edges_value, Method, Edges_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(Frequentist_table)


#ric model

SP_value_ric <- c()
SE_value_ric <- c()
MCC_value_ric <- c()

for (i in 1: reps) {
  SP_value_ric[i] = result_list[[i]]$SP_freq_ric
  SE_value_ric[i] = result_list[[i]]$SE_freq_ric
  MCC_value_ric[i] = result_list[[i]]$MCC_freq_ric
  
}

SE_Type = rep('Sensitivity', reps)
MCC_Type = rep('MCC', reps)

SP_Type = rep('Specificity', reps)


Edges_value =  c(SP_value_ric, SE_value_ric, MCC_value_ric)
Edges_Type = c(SP_Type, SE_Type, MCC_Type)

number_elements  = length(Edges_value)

Sparsity = rep('tenpercent', number_elements)
Method = rep('RIC', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Edges_value, Method, Edges_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(Frequentist_table)


#ebic model

SP_value_ebic <- c()
SE_value_ebic <- c()
MCC_value_ebic <- c()

for (i in 1: reps) {
  SP_value_ebic[i] = result_list[[i]]$SP_freq_ebic
  SE_value_ebic[i] = result_list[[i]]$SE_freq_ebic
  MCC_value_ebic[i] = result_list[[i]]$MCC_freq_ebic
  
}

SE_Type = rep('Sensitivity', reps)
MCC_Type = rep('MCC', reps)

SP_Type = rep('Specificity', reps)


Edges_value =  c(SP_value_ebic, SE_value_ebic, MCC_value_ebic)
Edges_Type = c(SP_Type, SE_Type, MCC_Type)

number_elements  = length(Edges_value)

Sparsity = rep('tenpercent', number_elements)
Method = rep('EBIC', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Edges_value, Method, Edges_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(Frequentist_table)

#Read in the p=50 data


fullTable_p50 <- read.csv(file = 'CholeskyDecomp_Boxplot_Edges_p50_n100.csv', header = TRUE, sep = ',')

fullTable_p50 <- as.data.frame(fullTable_p50)

fullTable_p50$Dimension <- as.factor(fullTable_p50$Dimension)

#Combine this to the other table

fullTable = rbind(fullTable, fullTable_p50)



#Read in the BDGraph data and add to the table

load('BDGraph_p50_n100_fivepercent.rdata')

SP_value <- c()
SE_value <- c()
MCC_value <- c()

for (i in 1: reps) {
  SP_value[i] = result_list[[i]]$SP
  SE_value[i] = result_list[[i]]$SE
  MCC_value[i] = result_list[[i]]$MCC
  
}

SE_Type = rep('Sensitivity', reps)
MCC_Type = rep('MCC', reps)

SP_Type = rep('Specificity', reps)


Edges_value = c(SP_value, SE_value, MCC_value)
Edges_Type = c(SP_Type, SE_Type, MCC_Type)

number_elements  = length(Edges_value)

Sparsity = rep('Percent', number_elements)
Method = rep('BDGraph', number_elements)
Dimension = rep(p, number_elements)

BDGraph_table = data.frame(Edges_value, Method, Edges_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, BDGraph_table)



rm(BDGraph_table)



#Read in the BDGraph data and add to the table

load('BDGraph_p50_n100_AR2.rdata')

SP_value <- c()
SE_value <- c()
MCC_value <- c()

for (i in 1: reps) {
  SP_value[i] = result_list[[i]]$SP
  SE_value[i] = result_list[[i]]$SE
  MCC_value[i] = result_list[[i]]$MCC
  
}

SE_Type = rep('Sensitivity', reps)
MCC_Type = rep('MCC', reps)

SP_Type = rep('Specificity', reps)


Edges_value = c(SP_value, SE_value, MCC_value)
Edges_Type = c(SP_Type, SE_Type, MCC_Type)

number_elements  = length(Edges_value)

Sparsity = rep('AR(2)', number_elements)
Method = rep('BDGraph', number_elements)
Dimension = rep(p, number_elements)

BDGraph_table = data.frame(Edges_value, Method, Edges_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, BDGraph_table)

rm(BDGraph_table)



#Read in the BDGraph data and add to the table

load('BDGraph_p50_n100_circle.rdata')

SP_value <- c()
SE_value <- c()
MCC_value <- c()

for (i in 1: reps) {
  SP_value[i] = result_list[[i]]$SP
  SE_value[i] = result_list[[i]]$SE
  MCC_value[i] = result_list[[i]]$MCC
  
}

SE_Type = rep('Sensitivity', reps)
MCC_Type = rep('MCC', reps)

SP_Type = rep('Specificity', reps)


Edges_value = c(SP_value, SE_value, MCC_value)
Edges_Type = c(SP_Type, SE_Type, MCC_Type)

number_elements  = length(Edges_value)

Sparsity = rep('Circle', number_elements)
Method = rep('BDGraph', number_elements)
Dimension = rep(p, number_elements)

BDGraph_table = data.frame(Edges_value, Method, Edges_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, BDGraph_table)

rm(BDGraph_table)

#Read in the Frequentist data and add to the table

load('Frequentist_p50_n100_circle.rdata')

#stars model

SP_value_stars <- c()
SE_value_stars <- c()
MCC_value_stars <- c()

for (i in 1: reps) {
  SP_value_stars[i] = result_list[[i]]$SP_freq_stars
  SE_value_stars[i] = result_list[[i]]$SE_freq_stars
  MCC_value_stars[i] = result_list[[i]]$MCC_freq_stars
  
}

SE_Type = rep('Sensitivity', reps)
MCC_Type = rep('MCC', reps)

SP_Type = rep('Specificity', reps)


Edges_value =  c(SP_value_stars, SE_value_stars, MCC_value_stars)
Edges_Type = c(SP_Type, SE_Type, MCC_Type)

number_elements  = length(Edges_value)

Sparsity = rep('Circle', number_elements)
Method = rep('StARS', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Edges_value, Method, Edges_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(Frequentist_table)


#ric model

SP_value_ric <- c()
SE_value_ric <- c()
MCC_value_ric <- c()

for (i in 1: reps) {
  SP_value_ric[i] = result_list[[i]]$SP_freq_ric
  SE_value_ric[i] = result_list[[i]]$SE_freq_ric
  MCC_value_ric[i] = result_list[[i]]$MCC_freq_ric
  
}

SE_Type = rep('Sensitivity', reps)
MCC_Type = rep('MCC', reps)

SP_Type = rep('Specificity', reps)


Edges_value =  c(SP_value_ric, SE_value_ric, MCC_value_ric)
Edges_Type = c(SP_Type, SE_Type, MCC_Type)

number_elements  = length(Edges_value)

Sparsity = rep('Circle', number_elements)
Method = rep('RIC', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Edges_value, Method, Edges_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(Frequentist_table)


#ebic model

SP_value_ebic <- c()
SE_value_ebic <- c()
MCC_value_ebic <- c()

for (i in 1: reps) {
  SP_value_ebic[i] = result_list[[i]]$SP_freq_ebic
  SE_value_ebic[i] = result_list[[i]]$SE_freq_ebic
  MCC_value_ebic[i] = result_list[[i]]$MCC_freq_ebic
  
}

SE_Type = rep('Sensitivity', reps)
MCC_Type = rep('MCC', reps)

SP_Type = rep('Specificity', reps)


Edges_value =  c(SP_value_ebic, SE_value_ebic, MCC_value_ebic)
Edges_Type = c(SP_Type, SE_Type, MCC_Type)

number_elements  = length(Edges_value)

Sparsity = rep('Circle', number_elements)
Method = rep('EBIC', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Edges_value, Method, Edges_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(Frequentist_table)



#Read in the Frequentist data and add to the table

load('Frequentist_p50_n100_AR2.rdata')

#stars model

SP_value_stars <- c()
SE_value_stars <- c()
MCC_value_stars <- c()

for (i in 1: reps) {
  SP_value_stars[i] = result_list[[i]]$SP_freq_stars
  SE_value_stars[i] = result_list[[i]]$SE_freq_stars
  MCC_value_stars[i] = result_list[[i]]$MCC_freq_stars
  
}

SE_Type = rep('Sensitivity', reps)
MCC_Type = rep('MCC', reps)

SP_Type = rep('Specificity', reps)


Edges_value =  c(SP_value_stars, SE_value_stars, MCC_value_stars)
Edges_Type = c(SP_Type, SE_Type, MCC_Type)

number_elements  = length(Edges_value)

Sparsity = rep('AR(2)', number_elements)
Method = rep('StARS', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Edges_value, Method, Edges_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(Frequentist_table)


#ric model

SP_value_ric <- c()
SE_value_ric <- c()
MCC_value_ric <- c()

for (i in 1: reps) {
  SP_value_ric[i] = result_list[[i]]$SP_freq_ric
  SE_value_ric[i] = result_list[[i]]$SE_freq_ric
  MCC_value_ric[i] = result_list[[i]]$MCC_freq_ric
  
}

SE_Type = rep('Sensitivity', reps)
MCC_Type = rep('MCC', reps)

SP_Type = rep('Specificity', reps)


Edges_value =  c(SP_value_ric, SE_value_ric, MCC_value_ric)
Edges_Type = c(SP_Type, SE_Type, MCC_Type)

number_elements  = length(Edges_value)

Sparsity = rep('AR(2)', number_elements)
Method = rep('RIC', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Edges_value, Method, Edges_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(Frequentist_table)


#ebic model

SP_value_ebic <- c()
SE_value_ebic <- c()
MCC_value_ebic <- c()

for (i in 1: reps) {
  SP_value_ebic[i] = result_list[[i]]$SP_freq_ebic
  SE_value_ebic[i] = result_list[[i]]$SE_freq_ebic
  MCC_value_ebic[i] = result_list[[i]]$MCC_freq_ebic
  
}

SE_Type = rep('Sensitivity', reps)
MCC_Type = rep('MCC', reps)

SP_Type = rep('Specificity', reps)


Edges_value =  c(SP_value_ebic, SE_value_ebic, MCC_value_ebic)
Edges_Type = c(SP_Type, SE_Type, MCC_Type)

number_elements  = length(Edges_value)

Sparsity = rep('AR(2)', number_elements)
Method = rep('EBIC', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Edges_value, Method, Edges_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(Frequentist_table)

#Read in the Frequentist data and add to the table

load('Frequentist_p50_n100_fivepercent.rdata')

#stars model

SP_value_stars <- c()
SE_value_stars <- c()
MCC_value_stars <- c()

for (i in 1: reps) {
  SP_value_stars[i] = result_list[[i]]$SP_freq_stars
  SE_value_stars[i] = result_list[[i]]$SE_freq_stars
  MCC_value_stars[i] = result_list[[i]]$MCC_freq_stars
  
}

SE_Type = rep('Sensitivity', reps)
MCC_Type = rep('MCC', reps)

SP_Type = rep('Specificity', reps)


Edges_value =  c(SP_value_stars, SE_value_stars, MCC_value_stars)
Edges_Type = c(SP_Type, SE_Type, MCC_Type)

number_elements  = length(Edges_value)

Sparsity = rep('fivepercent', number_elements)
Method = rep('StARS', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Edges_value, Method, Edges_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(Frequentist_table)


#ric model

SP_value_ric <- c()
SE_value_ric <- c()
MCC_value_ric <- c()

for (i in 1: reps) {
  SP_value_ric[i] = result_list[[i]]$SP_freq_ric
  SE_value_ric[i] = result_list[[i]]$SE_freq_ric
  MCC_value_ric[i] = result_list[[i]]$MCC_freq_ric
  
}

SE_Type = rep('Sensitivity', reps)
MCC_Type = rep('MCC', reps)

SP_Type = rep('Specificity', reps)


Edges_value =  c(SP_value_ric, SE_value_ric, MCC_value_ric)
Edges_Type = c(SP_Type, SE_Type, MCC_Type)

number_elements  = length(Edges_value)

Sparsity = rep('fivepercent', number_elements)
Method = rep('RIC', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Edges_value, Method, Edges_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(Frequentist_table)


#ebic model

SP_value_ebic <- c()
SE_value_ebic <- c()
MCC_value_ebic <- c()

for (i in 1: reps) {
  SP_value_ebic[i] = result_list[[i]]$SP_freq_ebic
  SE_value_ebic[i] = result_list[[i]]$SE_freq_ebic
  MCC_value_ebic[i] = result_list[[i]]$MCC_freq_ebic
  
}

SE_Type = rep('Sensitivity', reps)
MCC_Type = rep('MCC', reps)

SP_Type = rep('Specificity', reps)


Edges_value =  c(SP_value_ebic, SE_value_ebic, MCC_value_ebic)
Edges_Type = c(SP_Type, SE_Type, MCC_Type)

number_elements  = length(Edges_value)

Sparsity = rep('fivepercent', number_elements)
Method = rep('EBIC', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Edges_value, Method, Edges_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(Frequentist_table)


#Read in the p=100 data


fullTable_p100 <- read.csv(file = 'CholeskyDecomp_Boxplot_Edges_p100_n300.csv', header = TRUE, sep = ',')

fullTable_p100 <- as.data.frame(fullTable_p100)

fullTable_p100$Dimension <- as.factor(fullTable_p100$Dimension)

#Combine this to the other table

fullTable = rbind(fullTable, fullTable_p100)



#Read in the BDGraph data and add to the table

load('BDGraph_p100_n300_twopercent.rdata')

SP_value <- c()
SE_value <- c()
MCC_value <- c()

for (i in 1: reps) {
  SP_value[i] = result_list[[i]]$SP
  SE_value[i] = result_list[[i]]$SE
  MCC_value[i] = result_list[[i]]$MCC
  
}

SE_Type = rep('Sensitivity', reps)
MCC_Type = rep('MCC', reps)

SP_Type = rep('Specificity', reps)


Edges_value = c(SP_value, SE_value, MCC_value)
Edges_Type = c(SP_Type, SE_Type, MCC_Type)

number_elements  = length(Edges_value)

Sparsity = rep('Percent', number_elements)
Method = rep('BDGraph', number_elements)
Dimension = rep(p, number_elements)

BDGraph_table = data.frame(Edges_value, Method, Edges_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, BDGraph_table)



rm(BDGraph_table)



#Read in the BDGraph data and add to the table

load('BDGraph_p100_n300_AR2.rdata')

SP_value <- c()
SE_value <- c()
MCC_value <- c()

for (i in 1: reps) {
  SP_value[i] = result_list[[i]]$SP
  SE_value[i] = result_list[[i]]$SE
  MCC_value[i] = result_list[[i]]$MCC
  
}

SE_Type = rep('Sensitivity', reps)
MCC_Type = rep('MCC', reps)

SP_Type = rep('Specificity', reps)


Edges_value = c(SP_value, SE_value, MCC_value)
Edges_Type = c(SP_Type, SE_Type, MCC_Type)

number_elements  = length(Edges_value)

Sparsity = rep('AR(2)', number_elements)
Method = rep('BDGraph', number_elements)
Dimension = rep(p, number_elements)

BDGraph_table = data.frame(Edges_value, Method, Edges_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, BDGraph_table)

rm(BDGraph_table)



#Read in the BDGraph data and add to the table

load('BDGraph_p100_n300_circle.rdata')

SP_value <- c()
SE_value <- c()
MCC_value <- c()

for (i in 1: reps) {
  SP_value[i] = result_list[[i]]$SP
  SE_value[i] = result_list[[i]]$SE
  MCC_value[i] = result_list[[i]]$MCC
  
}

SE_Type = rep('Sensitivity', reps)
MCC_Type = rep('MCC', reps)

SP_Type = rep('Specificity', reps)


Edges_value = c(SP_value, SE_value, MCC_value)
Edges_Type = c(SP_Type, SE_Type, MCC_Type)

number_elements  = length(Edges_value)

Sparsity = rep('Circle', number_elements)
Method = rep('BDGraph', number_elements)
Dimension = rep(p, number_elements)

BDGraph_table = data.frame(Edges_value, Method, Edges_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, BDGraph_table)

rm(BDGraph_table)

#Read in the Frequentist data and add to the table

load('Frequentist_p100_n300_circle.rdata')

#stars model

SP_value_stars <- c()
SE_value_stars <- c()
MCC_value_stars <- c()

for (i in 1: reps) {
  SP_value_stars[i] = result_list[[i]]$SP_freq_stars
  SE_value_stars[i] = result_list[[i]]$SE_freq_stars
  MCC_value_stars[i] = result_list[[i]]$MCC_freq_stars
  
}

SE_Type = rep('Sensitivity', reps)
MCC_Type = rep('MCC', reps)

SP_Type = rep('Specificity', reps)


Edges_value =  c(SP_value_stars, SE_value_stars, MCC_value_stars)
Edges_Type = c(SP_Type, SE_Type, MCC_Type)

number_elements  = length(Edges_value)

Sparsity = rep('Circle', number_elements)
Method = rep('StARS', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Edges_value, Method, Edges_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(Frequentist_table)


#ric model

SP_value_ric <- c()
SE_value_ric <- c()
MCC_value_ric <- c()

for (i in 1: reps) {
  SP_value_ric[i] = result_list[[i]]$SP_freq_ric
  SE_value_ric[i] = result_list[[i]]$SE_freq_ric
  MCC_value_ric[i] = result_list[[i]]$MCC_freq_ric
  
}

SE_Type = rep('Sensitivity', reps)
MCC_Type = rep('MCC', reps)

SP_Type = rep('Specificity', reps)


Edges_value =  c(SP_value_ric, SE_value_ric, MCC_value_ric)
Edges_Type = c(SP_Type, SE_Type, MCC_Type)

number_elements  = length(Edges_value)

Sparsity = rep('Circle', number_elements)
Method = rep('RIC', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Edges_value, Method, Edges_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(Frequentist_table)


#ebic model

SP_value_ebic <- c()
SE_value_ebic <- c()
MCC_value_ebic <- c()

for (i in 1: reps) {
  SP_value_ebic[i] = result_list[[i]]$SP_freq_ebic
  SE_value_ebic[i] = result_list[[i]]$SE_freq_ebic
  MCC_value_ebic[i] = result_list[[i]]$MCC_freq_ebic
  
}

SE_Type = rep('Sensitivity', reps)
MCC_Type = rep('MCC', reps)

SP_Type = rep('Specificity', reps)


Edges_value =  c(SP_value_ebic, SE_value_ebic, MCC_value_ebic)
Edges_Type = c(SP_Type, SE_Type, MCC_Type)

number_elements  = length(Edges_value)

Sparsity = rep('Circle', number_elements)
Method = rep('EBIC', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Edges_value, Method, Edges_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(Frequentist_table)

#Read in the Frequentist data and add to the table

load('Frequentist_p100_n300_AR2.rdata')

#stars model

SP_value_stars <- c()
SE_value_stars <- c()
MCC_value_stars <- c()

for (i in 1: reps) {
  SP_value_stars[i] = result_list[[i]]$SP_freq_stars
  SE_value_stars[i] = result_list[[i]]$SE_freq_stars
  MCC_value_stars[i] = result_list[[i]]$MCC_freq_stars
  
}

SE_Type = rep('Sensitivity', reps)
MCC_Type = rep('MCC', reps)

SP_Type = rep('Specificity', reps)


Edges_value =  c(SP_value_stars, SE_value_stars, MCC_value_stars)
Edges_Type = c(SP_Type, SE_Type, MCC_Type)

number_elements  = length(Edges_value)

Sparsity = rep('AR(2)', number_elements)
Method = rep('StARS', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Edges_value, Method, Edges_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(Frequentist_table)


#ric model

SP_value_ric <- c()
SE_value_ric <- c()
MCC_value_ric <- c()

for (i in 1: reps) {
  SP_value_ric[i] = result_list[[i]]$SP_freq_ric
  SE_value_ric[i] = result_list[[i]]$SE_freq_ric
  MCC_value_ric[i] = result_list[[i]]$MCC_freq_ric
  
}

SE_Type = rep('Sensitivity', reps)
MCC_Type = rep('MCC', reps)

SP_Type = rep('Specificity', reps)


Edges_value =  c(SP_value_ric, SE_value_ric, MCC_value_ric)
Edges_Type = c(SP_Type, SE_Type, MCC_Type)

number_elements  = length(Edges_value)

Sparsity = rep('AR(2)', number_elements)
Method = rep('RIC', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Edges_value, Method, Edges_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(Frequentist_table)


#ebic model

SP_value_ebic <- c()
SE_value_ebic <- c()
MCC_value_ebic <- c()

for (i in 1: reps) {
  SP_value_ebic[i] = result_list[[i]]$SP_freq_ebic
  SE_value_ebic[i] = result_list[[i]]$SE_freq_ebic
  MCC_value_ebic[i] = result_list[[i]]$MCC_freq_ebic
  
}

SE_Type = rep('Sensitivity', reps)
MCC_Type = rep('MCC', reps)

SP_Type = rep('Specificity', reps)


Edges_value =  c(SP_value_ebic, SE_value_ebic, MCC_value_ebic)
Edges_Type = c(SP_Type, SE_Type, MCC_Type)

number_elements  = length(Edges_value)

Sparsity = rep('AR(2)', number_elements)
Method = rep('EBIC', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Edges_value, Method, Edges_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(Frequentist_table)
#Read in the Frequentist data and add to the table

load('Frequentist_p100_n300_twopercent.rdata')

#stars model

SP_value_stars <- c()
SE_value_stars <- c()
MCC_value_stars <- c()

for (i in 1: reps) {
  SP_value_stars[i] = result_list[[i]]$SP_freq_stars
  SE_value_stars[i] = result_list[[i]]$SE_freq_stars
  MCC_value_stars[i] = result_list[[i]]$MCC_freq_stars
  
}

SE_Type = rep('Sensitivity', reps)
MCC_Type = rep('MCC', reps)

SP_Type = rep('Specificity', reps)


Edges_value =  c(SP_value_stars, SE_value_stars, MCC_value_stars)
Edges_Type = c(SP_Type, SE_Type, MCC_Type)

number_elements  = length(Edges_value)

Sparsity = rep('Percent', number_elements)
Method = rep('StARS', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Edges_value, Method, Edges_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(Frequentist_table)


#ric model

SP_value_ric <- c()
SE_value_ric <- c()
MCC_value_ric <- c()

for (i in 1: reps) {
  SP_value_ric[i] = result_list[[i]]$SP_freq_ric
  SE_value_ric[i] = result_list[[i]]$SE_freq_ric
  MCC_value_ric[i] = result_list[[i]]$MCC_freq_ric
  
}

SE_Type = rep('Sensitivity', reps)
MCC_Type = rep('MCC', reps)

SP_Type = rep('Specificity', reps)


Edges_value =  c(SP_value_ric, SE_value_ric, MCC_value_ric)
Edges_Type = c(SP_Type, SE_Type, MCC_Type)

number_elements  = length(Edges_value)

Sparsity = rep('Percent', number_elements)
Method = rep('RIC', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Edges_value, Method, Edges_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(Frequentist_table)


#ebic model

SP_value_ebic <- c()
SE_value_ebic <- c()
MCC_value_ebic <- c()

for (i in 1: reps) {
  SP_value_ebic[i] = result_list[[i]]$SP_freq_ebic
  SE_value_ebic[i] = result_list[[i]]$SE_freq_ebic
  MCC_value_ebic[i] = result_list[[i]]$MCC_freq_ebic
  
}

SE_Type = rep('Sensitivity', reps)
MCC_Type = rep('MCC', reps)

SP_Type = rep('Specificity', reps)


Edges_value =  c(SP_value_ebic, SE_value_ebic, MCC_value_ebic)
Edges_Type = c(SP_Type, SE_Type, MCC_Type)

number_elements  = length(Edges_value)

Sparsity = rep('Percent', number_elements)
Method = rep('EBIC', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Edges_value, Method, Edges_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(Frequentist_table)

save(fullTable, file = 'CholeskyDecomp_Boxplot_Edges_table.rdata')
