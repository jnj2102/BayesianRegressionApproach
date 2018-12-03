####################################################
# Boxplots for the CholeskyDecomp paper
# Author: Jami Jackson Mulgrave
#
#
######################################################


#clear the workspace
rm(list = ls())


fullTable <- read.csv(file = 'CholeskyDecomp_Boxplot_Loss_p25_n25_notransform.csv', header = TRUE, sep = ',')

fullTable <- as.data.frame(fullTable)

fullTable$Dimension <- as.factor(fullTable$Dimension)


#Read in the BDGraph data and add to the table

load('BDGraph_p25_n25_tenpercent_nocopula.rdata')

ELoss_value <- c()
BLoss_value <- c()
FrobLoss_value <- c()

for (i in 1: reps) {
ELoss_value[i] = result_list[[i]]$entropy_loss_correlation
BLoss_value[i] = result_list[[i]]$bounded_loss
FrobLoss_value[i] = result_list[[i]]$FrobeniusLoss_precision

}

BLoss_Type = rep('Bounded Loss', reps)
FrobLoss_Type = rep('Frobenius Loss', reps)

ELoss_Type = rep('ELoss', reps)


Loss_value = c(ELoss_value, BLoss_value, FrobLoss_value)
Loss_Type = c(ELoss_Type, BLoss_Type, FrobLoss_Type)

number_elements  = length(Loss_value)

Sparsity = rep('Percent', number_elements)
Method = rep('BDGraph', number_elements)
Dimension = rep(p, number_elements)

BDGraph_table = data.frame(Loss_value, Method, Loss_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, BDGraph_table)



rm(BDGraph_table)


#Read in the BDGraph data and add to the table

load('BDGraph_p25_n25_AR2_nocopula.rdata')

ELoss_value <- c()
BLoss_value <- c()
FrobLoss_value <- c()

for (i in 1: reps) {
  ELoss_value[i] = result_list[[i]]$entropy_loss_correlation
  BLoss_value[i] = result_list[[i]]$bounded_loss
  FrobLoss_value[i] = result_list[[i]]$FrobeniusLoss_precision
  
}

BLoss_Type = rep('Bounded Loss', reps)
FrobLoss_Type = rep('Frobenius Loss', reps)

ELoss_Type = rep('ELoss', reps)


Loss_value = c(ELoss_value, BLoss_value, FrobLoss_value)
Loss_Type = c(ELoss_Type, BLoss_Type, FrobLoss_Type)

number_elements  = length(Loss_value)

Sparsity = rep('AR(2)', number_elements)
Method = rep('BDGraph', number_elements)
Dimension = rep(p, number_elements)

BDGraph_table = data.frame(Loss_value, Method, Loss_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, BDGraph_table)

rm(BDGraph_table)



#Read in the BDGraph data and add to the table

load('BDGraph_p25_n25_circle_nocopula.rdata')

ELoss_value <- c()
BLoss_value <- c()
FrobLoss_value <- c()

for (i in 1: reps) {
  ELoss_value[i] = result_list[[i]]$entropy_loss_correlation
  BLoss_value[i] = result_list[[i]]$bounded_loss
  FrobLoss_value[i] = result_list[[i]]$FrobeniusLoss_precision
  
}

BLoss_Type = rep('Bounded Loss', reps)
FrobLoss_Type = rep('Frobenius Loss', reps)

ELoss_Type = rep('ELoss', reps)


Loss_value = c(ELoss_value, BLoss_value, FrobLoss_value)
Loss_Type = c(ELoss_Type, BLoss_Type, FrobLoss_Type)

number_elements  = length(Loss_value)

Sparsity = rep('Circle', number_elements)
Method = rep('BDGraph', number_elements)
Dimension = rep(p, number_elements)

BDGraph_table = data.frame(Loss_value, Method, Loss_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, BDGraph_table)

rm(BDGraph_table)

#Read in the Frequentist data and add to the table

load('Frequentist_p25_n25_circle_notransform.rdata')

#stars model

ELoss_value_stars <- c()
BLoss_value_stars <- c()
FrobLoss_value_stars <- c()

for (i in 1: reps) {
  ELoss_value_stars[i] = result_list[[i]]$entropy_loss_stars
  BLoss_value_stars[i] = result_list[[i]]$bounded_loss_stars
  FrobLoss_value_stars[i] = result_list[[i]]$FrobeniusLoss_precision_stars
  
}

BLoss_Type = rep('Bounded Loss', reps)
FrobLoss_Type = rep('Frobenius Loss', reps)

ELoss_Type = rep('ELoss', reps)


Loss_value =  c(ELoss_value_stars, BLoss_value_stars, FrobLoss_value_stars)
Loss_Type = c(ELoss_Type, BLoss_Type, FrobLoss_Type)

number_elements  = length(Loss_value)

Sparsity = rep('Circle', number_elements)
Method = rep('StARS', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Loss_value, Method, Loss_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(Frequentist_table)


#ric model

ELoss_value_ric <- c()
BLoss_value_ric <- c()
FrobLoss_value_ric <- c()

for (i in 1: reps) {
  ELoss_value_ric[i] = result_list[[i]]$entropy_loss_ric
  BLoss_value_ric[i] = result_list[[i]]$bounded_loss_ric
  FrobLoss_value_ric[i] = result_list[[i]]$FrobeniusLoss_precision_ric
  
}

BLoss_Type = rep('Bounded Loss', reps)
FrobLoss_Type = rep('Frobenius Loss', reps)

ELoss_Type = rep('ELoss', reps)


Loss_value =  c(ELoss_value_ric, BLoss_value_ric, FrobLoss_value_ric)
Loss_Type = c(ELoss_Type, BLoss_Type, FrobLoss_Type)

number_elements  = length(Loss_value)

Sparsity = rep('Circle', number_elements)
Method = rep('RIC', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Loss_value, Method, Loss_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(Frequentist_table)


#ebic model

ELoss_value_ebic <- c()
BLoss_value_ebic <- c()
FrobLoss_value_ebic <- c()

for (i in 1: reps) {
  ELoss_value_ebic[i] = result_list[[i]]$entropy_loss_ebic
  BLoss_value_ebic[i] = result_list[[i]]$bounded_loss_ebic
  FrobLoss_value_ebic[i] = result_list[[i]]$FrobeniusLoss_precision_ebic
  
}

BLoss_Type = rep('Bounded Loss', reps)
FrobLoss_Type = rep('Frobenius Loss', reps)

ELoss_Type = rep('ELoss', reps)


Loss_value =  c(ELoss_value_ebic, BLoss_value_ebic, FrobLoss_value_ebic)
Loss_Type = c(ELoss_Type, BLoss_Type, FrobLoss_Type)

number_elements  = length(Loss_value)

Sparsity = rep('Circle', number_elements)
Method = rep('EBIC', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Loss_value, Method, Loss_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(Frequentist_table)


#Read in the Frequentist data and add to the table

load('Frequentist_p25_n25_AR2_notransform.rdata')

#stars model

ELoss_value_stars <- c()
BLoss_value_stars <- c()
FrobLoss_value_stars <- c()

for (i in 1: reps) {
  ELoss_value_stars[i] = result_list[[i]]$entropy_loss_stars
  BLoss_value_stars[i] = result_list[[i]]$bounded_loss_stars
  FrobLoss_value_stars[i] = result_list[[i]]$FrobeniusLoss_precision_stars
  
}

BLoss_Type = rep('Bounded Loss', reps)
FrobLoss_Type = rep('Frobenius Loss', reps)

ELoss_Type = rep('ELoss', reps)


Loss_value =  c(ELoss_value_stars, BLoss_value_stars, FrobLoss_value_stars)
Loss_Type = c(ELoss_Type, BLoss_Type, FrobLoss_Type)

number_elements  = length(Loss_value)

Sparsity = rep('AR(2)', number_elements)
Method = rep('StARS', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Loss_value, Method, Loss_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(Frequentist_table)


#ric model

ELoss_value_ric <- c()
BLoss_value_ric <- c()
FrobLoss_value_ric <- c()

for (i in 1: reps) {
  ELoss_value_ric[i] = result_list[[i]]$entropy_loss_ric
  BLoss_value_ric[i] = result_list[[i]]$bounded_loss_ric
  FrobLoss_value_ric[i] = result_list[[i]]$FrobeniusLoss_precision_ric
  
}

BLoss_Type = rep('Bounded Loss', reps)
FrobLoss_Type = rep('Frobenius Loss', reps)

ELoss_Type = rep('ELoss', reps)


Loss_value =  c(ELoss_value_ric, BLoss_value_ric, FrobLoss_value_ric)
Loss_Type = c(ELoss_Type, BLoss_Type, FrobLoss_Type)

number_elements  = length(Loss_value)

Sparsity = rep('AR(2)', number_elements)
Method = rep('RIC', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Loss_value, Method, Loss_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(Frequentist_table)


#ebic model

ELoss_value_ebic <- c()
BLoss_value_ebic <- c()
FrobLoss_value_ebic <- c()

for (i in 1: reps) {
  ELoss_value_ebic[i] = result_list[[i]]$entropy_loss_ebic
  BLoss_value_ebic[i] = result_list[[i]]$bounded_loss_ebic
  FrobLoss_value_ebic[i] = result_list[[i]]$FrobeniusLoss_precision_ebic
  
}

BLoss_Type = rep('Bounded Loss', reps)
FrobLoss_Type = rep('Frobenius Loss', reps)

ELoss_Type = rep('ELoss', reps)


Loss_value =  c(ELoss_value_ebic, BLoss_value_ebic, FrobLoss_value_ebic)
Loss_Type = c(ELoss_Type, BLoss_Type, FrobLoss_Type)

number_elements  = length(Loss_value)

Sparsity = rep('AR(2)', number_elements)
Method = rep('EBIC', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Loss_value, Method, Loss_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(Frequentist_table)

#Read in the Frequentist data and add to the table

load('Frequentist_p25_n25_tenpercent_notransform.rdata')

#stars model

ELoss_value_stars <- c()
BLoss_value_stars <- c()
FrobLoss_value_stars <- c()

for (i in 1: reps) {
  ELoss_value_stars[i] = result_list[[i]]$entropy_loss_stars
  BLoss_value_stars[i] = result_list[[i]]$bounded_loss_stars
  FrobLoss_value_stars[i] = result_list[[i]]$FrobeniusLoss_precision_stars
  
}

BLoss_Type = rep('Bounded Loss', reps)
FrobLoss_Type = rep('Frobenius Loss', reps)

ELoss_Type = rep('ELoss', reps)


Loss_value =  c(ELoss_value_stars, BLoss_value_stars, FrobLoss_value_stars)
Loss_Type = c(ELoss_Type, BLoss_Type, FrobLoss_Type)

number_elements  = length(Loss_value)

Sparsity = rep('tenpercent', number_elements)
Method = rep('StARS', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Loss_value, Method, Loss_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(Frequentist_table)


#ric model

ELoss_value_ric <- c()
BLoss_value_ric <- c()
FrobLoss_value_ric <- c()

for (i in 1: reps) {
  ELoss_value_ric[i] = result_list[[i]]$entropy_loss_ric
  BLoss_value_ric[i] = result_list[[i]]$bounded_loss_ric
  FrobLoss_value_ric[i] = result_list[[i]]$FrobeniusLoss_precision_ric
  
}

BLoss_Type = rep('Bounded Loss', reps)
FrobLoss_Type = rep('Frobenius Loss', reps)

ELoss_Type = rep('ELoss', reps)


Loss_value =  c(ELoss_value_ric, BLoss_value_ric, FrobLoss_value_ric)
Loss_Type = c(ELoss_Type, BLoss_Type, FrobLoss_Type)

number_elements  = length(Loss_value)

Sparsity = rep('tenpercent', number_elements)
Method = rep('RIC', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Loss_value, Method, Loss_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(Frequentist_table)


#ebic model

ELoss_value_ebic <- c()
BLoss_value_ebic <- c()
FrobLoss_value_ebic <- c()

for (i in 1: reps) {
  ELoss_value_ebic[i] = result_list[[i]]$entropy_loss_ebic
  BLoss_value_ebic[i] = result_list[[i]]$bounded_loss_ebic
  FrobLoss_value_ebic[i] = result_list[[i]]$FrobeniusLoss_precision_ebic
  
}

BLoss_Type = rep('Bounded Loss', reps)
FrobLoss_Type = rep('Frobenius Loss', reps)

ELoss_Type = rep('ELoss', reps)


Loss_value =  c(ELoss_value_ebic, BLoss_value_ebic, FrobLoss_value_ebic)
Loss_Type = c(ELoss_Type, BLoss_Type, FrobLoss_Type)

number_elements  = length(Loss_value)

Sparsity = rep('tenpercent', number_elements)
Method = rep('EBIC', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Loss_value, Method, Loss_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(Frequentist_table)

#Read in the p=50 data


fullTable_p50 <- read.csv(file = 'CholeskyDecomp_Boxplot_Loss_p50_n100_notransform.csv', header = TRUE, sep = ',')

fullTable_p50 <- as.data.frame(fullTable_p50)

fullTable_p50$Dimension <- as.factor(fullTable_p50$Dimension)

#Combine this to the other table

fullTable = rbind(fullTable, fullTable_p50)



#Read in the BDGraph data and add to the table

load('BDGraph_p50_n100_fivepercent_nocopula.rdata')

ELoss_value <- c()
BLoss_value <- c()
FrobLoss_value <- c()

for (i in 1: reps) {
  ELoss_value[i] = result_list[[i]]$entropy_loss_correlation
  BLoss_value[i] = result_list[[i]]$bounded_loss
  FrobLoss_value[i] = result_list[[i]]$FrobeniusLoss_precision
  
}

BLoss_Type = rep('Bounded Loss', reps)
FrobLoss_Type = rep('Frobenius Loss', reps)

ELoss_Type = rep('ELoss', reps)


Loss_value = c(ELoss_value, BLoss_value, FrobLoss_value)
Loss_Type = c(ELoss_Type, BLoss_Type, FrobLoss_Type)

number_elements  = length(Loss_value)

Sparsity = rep('Percent', number_elements)
Method = rep('BDGraph', number_elements)
Dimension = rep(p, number_elements)

BDGraph_table = data.frame(Loss_value, Method, Loss_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, BDGraph_table)



rm(BDGraph_table)



#Read in the BDGraph data and add to the table

load('BDGraph_p50_n100_AR2_nocopula.rdata')

ELoss_value <- c()
BLoss_value <- c()
FrobLoss_value <- c()

for (i in 1: reps) {
  ELoss_value[i] = result_list[[i]]$entropy_loss_correlation
  BLoss_value[i] = result_list[[i]]$bounded_loss
  FrobLoss_value[i] = result_list[[i]]$FrobeniusLoss_precision
  
}

BLoss_Type = rep('Bounded Loss', reps)
FrobLoss_Type = rep('Frobenius Loss', reps)

ELoss_Type = rep('ELoss', reps)


Loss_value = c(ELoss_value, BLoss_value, FrobLoss_value)
Loss_Type = c(ELoss_Type, BLoss_Type, FrobLoss_Type)

number_elements  = length(Loss_value)

Sparsity = rep('AR(2)', number_elements)
Method = rep('BDGraph', number_elements)
Dimension = rep(p, number_elements)

BDGraph_table = data.frame(Loss_value, Method, Loss_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, BDGraph_table)

rm(BDGraph_table)



#Read in the BDGraph data and add to the table

load('BDGraph_p50_n100_circle_nocopula.rdata')

ELoss_value <- c()
BLoss_value <- c()
FrobLoss_value <- c()

for (i in 1: reps) {
  ELoss_value[i] = result_list[[i]]$entropy_loss_correlation
  BLoss_value[i] = result_list[[i]]$bounded_loss
  FrobLoss_value[i] = result_list[[i]]$FrobeniusLoss_precision
  
}

BLoss_Type = rep('Bounded Loss', reps)
FrobLoss_Type = rep('Frobenius Loss', reps)

ELoss_Type = rep('ELoss', reps)


Loss_value = c(ELoss_value, BLoss_value, FrobLoss_value)
Loss_Type = c(ELoss_Type, BLoss_Type, FrobLoss_Type)

number_elements  = length(Loss_value)

Sparsity = rep('Circle', number_elements)
Method = rep('BDGraph', number_elements)
Dimension = rep(p, number_elements)

BDGraph_table = data.frame(Loss_value, Method, Loss_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, BDGraph_table)

rm(BDGraph_table)

#Read in the Frequentist data and add to the table

load('Frequentist_p50_n100_circle_notransform.rdata')

#stars model

ELoss_value_stars <- c()
BLoss_value_stars <- c()
FrobLoss_value_stars <- c()

for (i in 1: reps) {
  ELoss_value_stars[i] = result_list[[i]]$entropy_loss_stars
  BLoss_value_stars[i] = result_list[[i]]$bounded_loss_stars
  FrobLoss_value_stars[i] = result_list[[i]]$FrobeniusLoss_precision_stars
  
}

BLoss_Type = rep('Bounded Loss', reps)
FrobLoss_Type = rep('Frobenius Loss', reps)

ELoss_Type = rep('ELoss', reps)


Loss_value =  c(ELoss_value_stars, BLoss_value_stars, FrobLoss_value_stars)
Loss_Type = c(ELoss_Type, BLoss_Type, FrobLoss_Type)

number_elements  = length(Loss_value)

Sparsity = rep('Circle', number_elements)
Method = rep('StARS', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Loss_value, Method, Loss_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(Frequentist_table)


#ric model

ELoss_value_ric <- c()
BLoss_value_ric <- c()
FrobLoss_value_ric <- c()

for (i in 1: reps) {
  ELoss_value_ric[i] = result_list[[i]]$entropy_loss_ric
  BLoss_value_ric[i] = result_list[[i]]$bounded_loss_ric
  FrobLoss_value_ric[i] = result_list[[i]]$FrobeniusLoss_precision_ric
  
}

BLoss_Type = rep('Bounded Loss', reps)
FrobLoss_Type = rep('Frobenius Loss', reps)

ELoss_Type = rep('ELoss', reps)


Loss_value =  c(ELoss_value_ric, BLoss_value_ric, FrobLoss_value_ric)
Loss_Type = c(ELoss_Type, BLoss_Type, FrobLoss_Type)

number_elements  = length(Loss_value)

Sparsity = rep('Circle', number_elements)
Method = rep('RIC', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Loss_value, Method, Loss_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(Frequentist_table)


#ebic model

ELoss_value_ebic <- c()
BLoss_value_ebic <- c()
FrobLoss_value_ebic <- c()

for (i in 1: reps) {
  ELoss_value_ebic[i] = result_list[[i]]$entropy_loss_ebic
  BLoss_value_ebic[i] = result_list[[i]]$bounded_loss_ebic
  FrobLoss_value_ebic[i] = result_list[[i]]$FrobeniusLoss_precision_ebic
  
}

BLoss_Type = rep('Bounded Loss', reps)
FrobLoss_Type = rep('Frobenius Loss', reps)

ELoss_Type = rep('ELoss', reps)


Loss_value =  c(ELoss_value_ebic, BLoss_value_ebic, FrobLoss_value_ebic)
Loss_Type = c(ELoss_Type, BLoss_Type, FrobLoss_Type)

number_elements  = length(Loss_value)

Sparsity = rep('Circle', number_elements)
Method = rep('EBIC', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Loss_value, Method, Loss_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(Frequentist_table)



#Read in the Frequentist data and add to the table

load('Frequentist_p50_n100_AR2_notransform.rdata')

#stars model

ELoss_value_stars <- c()
BLoss_value_stars <- c()
FrobLoss_value_stars <- c()

for (i in 1: reps) {
  ELoss_value_stars[i] = result_list[[i]]$entropy_loss_stars
  BLoss_value_stars[i] = result_list[[i]]$bounded_loss_stars
  FrobLoss_value_stars[i] = result_list[[i]]$FrobeniusLoss_precision_stars
  
}

BLoss_Type = rep('Bounded Loss', reps)
FrobLoss_Type = rep('Frobenius Loss', reps)

ELoss_Type = rep('ELoss', reps)


Loss_value =  c(ELoss_value_stars, BLoss_value_stars, FrobLoss_value_stars)
Loss_Type = c(ELoss_Type, BLoss_Type, FrobLoss_Type)

number_elements  = length(Loss_value)

Sparsity = rep('AR(2)', number_elements)
Method = rep('StARS', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Loss_value, Method, Loss_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(Frequentist_table)


#ric model

ELoss_value_ric <- c()
BLoss_value_ric <- c()
FrobLoss_value_ric <- c()

for (i in 1: reps) {
  ELoss_value_ric[i] = result_list[[i]]$entropy_loss_ric
  BLoss_value_ric[i] = result_list[[i]]$bounded_loss_ric
  FrobLoss_value_ric[i] = result_list[[i]]$FrobeniusLoss_precision_ric
  
}

BLoss_Type = rep('Bounded Loss', reps)
FrobLoss_Type = rep('Frobenius Loss', reps)

ELoss_Type = rep('ELoss', reps)


Loss_value =  c(ELoss_value_ric, BLoss_value_ric, FrobLoss_value_ric)
Loss_Type = c(ELoss_Type, BLoss_Type, FrobLoss_Type)

number_elements  = length(Loss_value)

Sparsity = rep('AR(2)', number_elements)
Method = rep('RIC', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Loss_value, Method, Loss_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(Frequentist_table)


#ebic model

ELoss_value_ebic <- c()
BLoss_value_ebic <- c()
FrobLoss_value_ebic <- c()

for (i in 1: reps) {
  ELoss_value_ebic[i] = result_list[[i]]$entropy_loss_ebic
  BLoss_value_ebic[i] = result_list[[i]]$bounded_loss_ebic
  FrobLoss_value_ebic[i] = result_list[[i]]$FrobeniusLoss_precision_ebic
  
}

BLoss_Type = rep('Bounded Loss', reps)
FrobLoss_Type = rep('Frobenius Loss', reps)

ELoss_Type = rep('ELoss', reps)


Loss_value =  c(ELoss_value_ebic, BLoss_value_ebic, FrobLoss_value_ebic)
Loss_Type = c(ELoss_Type, BLoss_Type, FrobLoss_Type)

number_elements  = length(Loss_value)

Sparsity = rep('AR(2)', number_elements)
Method = rep('EBIC', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Loss_value, Method, Loss_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(Frequentist_table)

#Read in the Frequentist data and add to the table

load('Frequentist_p50_n100_fivepercent_notransform.rdata')

#stars model

ELoss_value_stars <- c()
BLoss_value_stars <- c()
FrobLoss_value_stars <- c()

for (i in 1: reps) {
  ELoss_value_stars[i] = result_list[[i]]$entropy_loss_stars
  BLoss_value_stars[i] = result_list[[i]]$bounded_loss_stars
  FrobLoss_value_stars[i] = result_list[[i]]$FrobeniusLoss_precision_stars
  
}

BLoss_Type = rep('Bounded Loss', reps)
FrobLoss_Type = rep('Frobenius Loss', reps)

ELoss_Type = rep('ELoss', reps)


Loss_value =  c(ELoss_value_stars, BLoss_value_stars, FrobLoss_value_stars)
Loss_Type = c(ELoss_Type, BLoss_Type, FrobLoss_Type)

number_elements  = length(Loss_value)

Sparsity = rep('fivepercent', number_elements)
Method = rep('StARS', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Loss_value, Method, Loss_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(Frequentist_table)


#ric model

ELoss_value_ric <- c()
BLoss_value_ric <- c()
FrobLoss_value_ric <- c()

for (i in 1: reps) {
  ELoss_value_ric[i] = result_list[[i]]$entropy_loss_ric
  BLoss_value_ric[i] = result_list[[i]]$bounded_loss_ric
  FrobLoss_value_ric[i] = result_list[[i]]$FrobeniusLoss_precision_ric
  
}

BLoss_Type = rep('Bounded Loss', reps)
FrobLoss_Type = rep('Frobenius Loss', reps)

ELoss_Type = rep('ELoss', reps)


Loss_value =  c(ELoss_value_ric, BLoss_value_ric, FrobLoss_value_ric)
Loss_Type = c(ELoss_Type, BLoss_Type, FrobLoss_Type)

number_elements  = length(Loss_value)

Sparsity = rep('fivepercent', number_elements)
Method = rep('RIC', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Loss_value, Method, Loss_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(Frequentist_table)


#ebic model

ELoss_value_ebic <- c()
BLoss_value_ebic <- c()
FrobLoss_value_ebic <- c()

for (i in 1: reps) {
  ELoss_value_ebic[i] = result_list[[i]]$entropy_loss_ebic
  BLoss_value_ebic[i] = result_list[[i]]$bounded_loss_ebic
  FrobLoss_value_ebic[i] = result_list[[i]]$FrobeniusLoss_precision_ebic
  
}

BLoss_Type = rep('Bounded Loss', reps)
FrobLoss_Type = rep('Frobenius Loss', reps)

ELoss_Type = rep('ELoss', reps)


Loss_value =  c(ELoss_value_ebic, BLoss_value_ebic, FrobLoss_value_ebic)
Loss_Type = c(ELoss_Type, BLoss_Type, FrobLoss_Type)

number_elements  = length(Loss_value)

Sparsity = rep('fivepercent', number_elements)
Method = rep('EBIC', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Loss_value, Method, Loss_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(Frequentist_table)


#Read in the p=100 data


fullTable_p100 <- read.csv(file = 'CholeskyDecomp_Boxplot_Loss_p100_n300_notransform.csv', header = TRUE, sep = ',')

fullTable_p100 <- as.data.frame(fullTable_p100)

fullTable_p100$Dimension <- as.factor(fullTable_p100$Dimension)

#Combine this to the other table

fullTable = rbind(fullTable, fullTable_p100)



#Read in the BDGraph data and add to the table

load('BDGraph_p100_n300_twopercent_nocopula.rdata')

ELoss_value <- c()
BLoss_value <- c()
FrobLoss_value <- c()

for (i in 1: reps) {
  ELoss_value[i] = result_list[[i]]$entropy_loss_correlation
  BLoss_value[i] = result_list[[i]]$bounded_loss
  FrobLoss_value[i] = result_list[[i]]$FrobeniusLoss_precision
  
}

BLoss_Type = rep('Bounded Loss', reps)
FrobLoss_Type = rep('Frobenius Loss', reps)

ELoss_Type = rep('ELoss', reps)


Loss_value = c(ELoss_value, BLoss_value, FrobLoss_value)
Loss_Type = c(ELoss_Type, BLoss_Type, FrobLoss_Type)

number_elements  = length(Loss_value)

Sparsity = rep('Percent', number_elements)
Method = rep('BDGraph', number_elements)
Dimension = rep(p, number_elements)

BDGraph_table = data.frame(Loss_value, Method, Loss_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, BDGraph_table)



rm(BDGraph_table)



#Read in the BDGraph data and add to the table

load('BDGraph_p100_n300_AR2_nocopula.rdata')

ELoss_value <- c()
BLoss_value <- c()
FrobLoss_value <- c()

for (i in 1: reps) {
  ELoss_value[i] = result_list[[i]]$entropy_loss_correlation
  BLoss_value[i] = result_list[[i]]$bounded_loss
  FrobLoss_value[i] = result_list[[i]]$FrobeniusLoss_precision
  
}

BLoss_Type = rep('Bounded Loss', reps)
FrobLoss_Type = rep('Frobenius Loss', reps)

ELoss_Type = rep('ELoss', reps)


Loss_value = c(ELoss_value, BLoss_value, FrobLoss_value)
Loss_Type = c(ELoss_Type, BLoss_Type, FrobLoss_Type)

number_elements  = length(Loss_value)

Sparsity = rep('AR(2)', number_elements)
Method = rep('BDGraph', number_elements)
Dimension = rep(p, number_elements)

BDGraph_table = data.frame(Loss_value, Method, Loss_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, BDGraph_table)

rm(BDGraph_table)



#Read in the BDGraph data and add to the table

load('BDGraph_p100_n300_circle_nocopula.rdata')

ELoss_value <- c()
BLoss_value <- c()
FrobLoss_value <- c()

for (i in 1: reps) {
  ELoss_value[i] = result_list[[i]]$entropy_loss_correlation
  BLoss_value[i] = result_list[[i]]$bounded_loss
  FrobLoss_value[i] = result_list[[i]]$FrobeniusLoss_precision
  
}

BLoss_Type = rep('Bounded Loss', reps)
FrobLoss_Type = rep('Frobenius Loss', reps)

ELoss_Type = rep('ELoss', reps)


Loss_value = c(ELoss_value, BLoss_value, FrobLoss_value)
Loss_Type = c(ELoss_Type, BLoss_Type, FrobLoss_Type)

number_elements  = length(Loss_value)

Sparsity = rep('Circle', number_elements)
Method = rep('BDGraph', number_elements)
Dimension = rep(p, number_elements)

BDGraph_table = data.frame(Loss_value, Method, Loss_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, BDGraph_table)

rm(BDGraph_table)

#Read in the Frequentist data and add to the table

load('Frequentist_p100_n300_circle_notransform.rdata')

#stars model

ELoss_value_stars <- c()
BLoss_value_stars <- c()
FrobLoss_value_stars <- c()

for (i in 1: reps) {
  ELoss_value_stars[i] = result_list[[i]]$entropy_loss_stars
  BLoss_value_stars[i] = result_list[[i]]$bounded_loss_stars
  FrobLoss_value_stars[i] = result_list[[i]]$FrobeniusLoss_precision_stars
  
}

BLoss_Type = rep('Bounded Loss', reps)
FrobLoss_Type = rep('Frobenius Loss', reps)

ELoss_Type = rep('ELoss', reps)


Loss_value =  c(ELoss_value_stars, BLoss_value_stars, FrobLoss_value_stars)
Loss_Type = c(ELoss_Type, BLoss_Type, FrobLoss_Type)

number_elements  = length(Loss_value)

Sparsity = rep('Circle', number_elements)
Method = rep('StARS', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Loss_value, Method, Loss_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(Frequentist_table)


#ric model

ELoss_value_ric <- c()
BLoss_value_ric <- c()
FrobLoss_value_ric <- c()

for (i in 1: reps) {
  ELoss_value_ric[i] = result_list[[i]]$entropy_loss_ric
  BLoss_value_ric[i] = result_list[[i]]$bounded_loss_ric
  FrobLoss_value_ric[i] = result_list[[i]]$FrobeniusLoss_precision_ric
  
}

BLoss_Type = rep('Bounded Loss', reps)
FrobLoss_Type = rep('Frobenius Loss', reps)

ELoss_Type = rep('ELoss', reps)


Loss_value =  c(ELoss_value_ric, BLoss_value_ric, FrobLoss_value_ric)
Loss_Type = c(ELoss_Type, BLoss_Type, FrobLoss_Type)

number_elements  = length(Loss_value)

Sparsity = rep('Circle', number_elements)
Method = rep('RIC', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Loss_value, Method, Loss_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(Frequentist_table)


#ebic model

ELoss_value_ebic <- c()
BLoss_value_ebic <- c()
FrobLoss_value_ebic <- c()

for (i in 1: reps) {
  ELoss_value_ebic[i] = result_list[[i]]$entropy_loss_ebic
  BLoss_value_ebic[i] = result_list[[i]]$bounded_loss_ebic
  FrobLoss_value_ebic[i] = result_list[[i]]$FrobeniusLoss_precision_ebic
  
}

BLoss_Type = rep('Bounded Loss', reps)
FrobLoss_Type = rep('Frobenius Loss', reps)

ELoss_Type = rep('ELoss', reps)


Loss_value =  c(ELoss_value_ebic, BLoss_value_ebic, FrobLoss_value_ebic)
Loss_Type = c(ELoss_Type, BLoss_Type, FrobLoss_Type)

number_elements  = length(Loss_value)

Sparsity = rep('Circle', number_elements)
Method = rep('EBIC', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Loss_value, Method, Loss_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(Frequentist_table)

#Read in the Frequentist data and add to the table

load('Frequentist_p100_n300_AR2_notransform.rdata')

#stars model

ELoss_value_stars <- c()
BLoss_value_stars <- c()
FrobLoss_value_stars <- c()

for (i in 1: reps) {
  ELoss_value_stars[i] = result_list[[i]]$entropy_loss_stars
  BLoss_value_stars[i] = result_list[[i]]$bounded_loss_stars
  FrobLoss_value_stars[i] = result_list[[i]]$FrobeniusLoss_precision_stars
  
}

BLoss_Type = rep('Bounded Loss', reps)
FrobLoss_Type = rep('Frobenius Loss', reps)

ELoss_Type = rep('ELoss', reps)


Loss_value =  c(ELoss_value_stars, BLoss_value_stars, FrobLoss_value_stars)
Loss_Type = c(ELoss_Type, BLoss_Type, FrobLoss_Type)

number_elements  = length(Loss_value)

Sparsity = rep('AR(2)', number_elements)
Method = rep('StARS', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Loss_value, Method, Loss_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(Frequentist_table)


#ric model

ELoss_value_ric <- c()
BLoss_value_ric <- c()
FrobLoss_value_ric <- c()

for (i in 1: reps) {
  ELoss_value_ric[i] = result_list[[i]]$entropy_loss_ric
  BLoss_value_ric[i] = result_list[[i]]$bounded_loss_ric
  FrobLoss_value_ric[i] = result_list[[i]]$FrobeniusLoss_precision_ric
  
}

BLoss_Type = rep('Bounded Loss', reps)
FrobLoss_Type = rep('Frobenius Loss', reps)

ELoss_Type = rep('ELoss', reps)


Loss_value =  c(ELoss_value_ric, BLoss_value_ric, FrobLoss_value_ric)
Loss_Type = c(ELoss_Type, BLoss_Type, FrobLoss_Type)

number_elements  = length(Loss_value)

Sparsity = rep('AR(2)', number_elements)
Method = rep('RIC', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Loss_value, Method, Loss_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(Frequentist_table)


#ebic model

ELoss_value_ebic <- c()
BLoss_value_ebic <- c()
FrobLoss_value_ebic <- c()

for (i in 1: reps) {
  ELoss_value_ebic[i] = result_list[[i]]$entropy_loss_ebic
  BLoss_value_ebic[i] = result_list[[i]]$bounded_loss_ebic
  FrobLoss_value_ebic[i] = result_list[[i]]$FrobeniusLoss_precision_ebic
  
}

BLoss_Type = rep('Bounded Loss', reps)
FrobLoss_Type = rep('Frobenius Loss', reps)

ELoss_Type = rep('ELoss', reps)


Loss_value =  c(ELoss_value_ebic, BLoss_value_ebic, FrobLoss_value_ebic)
Loss_Type = c(ELoss_Type, BLoss_Type, FrobLoss_Type)

number_elements  = length(Loss_value)

Sparsity = rep('AR(2)', number_elements)
Method = rep('EBIC', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Loss_value, Method, Loss_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(Frequentist_table)
#Read in the Frequentist data and add to the table

load('Frequentist_p100_n300_twopercent_notransform.rdata')

#stars model

ELoss_value_stars <- c()
BLoss_value_stars <- c()
FrobLoss_value_stars <- c()

for (i in 1: reps) {
  ELoss_value_stars[i] = result_list[[i]]$entropy_loss_stars
  BLoss_value_stars[i] = result_list[[i]]$bounded_loss_stars
  FrobLoss_value_stars[i] = result_list[[i]]$FrobeniusLoss_precision_stars
  
}

BLoss_Type = rep('Bounded Loss', reps)
FrobLoss_Type = rep('Frobenius Loss', reps)

ELoss_Type = rep('ELoss', reps)


Loss_value =  c(ELoss_value_stars, BLoss_value_stars, FrobLoss_value_stars)
Loss_Type = c(ELoss_Type, BLoss_Type, FrobLoss_Type)

number_elements  = length(Loss_value)

Sparsity = rep('Percent', number_elements)
Method = rep('StARS', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Loss_value, Method, Loss_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(Frequentist_table)


#ric model

ELoss_value_ric <- c()
BLoss_value_ric <- c()
FrobLoss_value_ric <- c()

for (i in 1: reps) {
  ELoss_value_ric[i] = result_list[[i]]$entropy_loss_ric
  BLoss_value_ric[i] = result_list[[i]]$bounded_loss_ric
  FrobLoss_value_ric[i] = result_list[[i]]$FrobeniusLoss_precision_ric
  
}

BLoss_Type = rep('Bounded Loss', reps)
FrobLoss_Type = rep('Frobenius Loss', reps)

ELoss_Type = rep('ELoss', reps)


Loss_value =  c(ELoss_value_ric, BLoss_value_ric, FrobLoss_value_ric)
Loss_Type = c(ELoss_Type, BLoss_Type, FrobLoss_Type)

number_elements  = length(Loss_value)

Sparsity = rep('Percent', number_elements)
Method = rep('RIC', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Loss_value, Method, Loss_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(Frequentist_table)


#ebic model

ELoss_value_ebic <- c()
BLoss_value_ebic <- c()
FrobLoss_value_ebic <- c()

for (i in 1: reps) {
  ELoss_value_ebic[i] = result_list[[i]]$entropy_loss_ebic
  BLoss_value_ebic[i] = result_list[[i]]$bounded_loss_ebic
  FrobLoss_value_ebic[i] = result_list[[i]]$FrobeniusLoss_precision_ebic
  
}

BLoss_Type = rep('Bounded Loss', reps)
FrobLoss_Type = rep('Frobenius Loss', reps)

ELoss_Type = rep('ELoss', reps)


Loss_value =  c(ELoss_value_ebic, BLoss_value_ebic, FrobLoss_value_ebic)
Loss_Type = c(ELoss_Type, BLoss_Type, FrobLoss_Type)

number_elements  = length(Loss_value)

Sparsity = rep('Percent', number_elements)
Method = rep('EBIC', number_elements)
Dimension = rep(p, number_elements)

Frequentist_table = data.frame(Loss_value, Method, Loss_Type, Sparsity, Dimension)

fullTable = rbind(fullTable, Frequentist_table)

rm(Frequentist_table)

save(fullTable, file = 'CholeskyDecomp_Boxplot_Loss_table_notransform.rdata')
