#day8
# BACKGROUNDCanonical methods are used to determine thecommon structure,or correspondence,
# between two sets of variables (matrices) measured on the same
# sampling units (Legendre andLegendre 1998). 
#A typical application involves determiningthe degree and nature of the relationship
# between speciesand environmental variables or the examination of the relationship betweentwo 
# species data sets (e.g., plants and insects). However, the possibilities are endless. 
# In all cases, the objective is toexamine and summarize the common structure of the twodata
# sets (Dray et al. 2003).In constrained (canonical) ordination analysesonly the variation 
# in the Y matrix (e.g., object-by-species)that can be explained 
# by theX matrix (e.g., object-by-environment)is displayed and analyzed, and
# notall the variation in the Y matrix.In other words, constrained ordination provides 
# a method of explicitly testing whether an explanatory matrix X significantly
# explains variation in a response matrix Y.  
# The method is best conceptualized as an extension of linear regression in which our 
# Y matrix corresponds to our response variable and our predictor variables are contained in matrix X.  
# Linear combinations of the explanatory variables are generated in a manner that explains 
# as much variance as possible in the response matrix. 

#load multiple pkgs
lapply(c("vegan"), require, character.only = TRUE)

#source
source('biostats.R')

#Let’s import the data, by typing:
speabu <-read.csv('./data/MAHA_speciesabu.csv',header=TRUE, row.names=1)
envdata <-read.csv('./data/MAHA_environment.csv',header=TRUE, row.names=1)

#Now let’s transform both datasets before diving into the analysis. First, 
#let’s conduct a ‘row normalization’ (i.e., to rescale each row so that the sums 
#equal1, or in other words to calculate relative species abundances), 
#type the following:
speabu.tran<-data.stand(speabu,method='total',margin='row',plot=F)

#Next, let’s log (based 10) transform the environmental dataset.
envdata.tran<-data.trans(envdata,method='log',plot=F)

spe.cca<-cca(speabu.tran~.,envdata.tran) 
summary(spe.cca)


