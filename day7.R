#Day 7
#RSource:biostats, evplot

#load multiple pkgs
lapply(c("vegan", "FactoMineR", "factoextra"), require, character.only = TRUE)

#source
source('biostats.R')
source('evplot.R')

speabu <- read.csv('./data/MAHA_speciesabu.csv',header=TRUE, row.names=1)
envdata <- read.csv('./data/MAHA_environment.csv',header=TRUE, row.names=1) 

#speocc <- data.trans(speabu,method='power',exp=0,plot=F)
#Alternatively we can directly import the data, by typing:
speocc <- read.csv('./data/MAHA_speciesocc.csv',header=TRUE, row.names=1)

#Let’s perform CA on the species presence-absence matrix, by typing:
  spe.ca<-cca(speocc)

#   The results of a CA are complex and stored in a list, including:
#     • eigenvalues (inertia) associated with each axis
#   • eigenvectors (sample and species scores)
# #igenvectors (sample and species scores)
# to view the summary of the eigenvalues. Your results should be as below.
summary(spe.ca)

# Recall that in PCA eigenvalues represent ‘variance’ explained. 
# In CA, eigenvalues represent ‘inertia’, where total inertia equals
# the chi-squared statistic of the data matrix standardized to unit total. 
# Thus, the eigenvalues are something akin to variance, but they are not variances exactly. 
# The chi-square statistic is a measure of association between samples and species.

#To see the inertia (eigenvalue) of each axis, type:
spe.ca$CA$eig
#The total inertia is equal to the sum of all the eigenvalues, and 
#you can get this by typing either one of the commands below:
    sum(spe.ca$CA$eig)
    spe.ca$CA$tot.chi
#Recall that you divide each eigenvalue by the sum of 
#eigenvalues to calculate the proportion of variation accounted for by each axis. 
#To calculate this as a percentage you can type:
spe.ca$CA$eig/sum(spe.ca$CA$eig)*100
#We may wish to test the statistical significance of the first several eigenvalues 
#      using a randomization test, as follows:
ordi.monte(speocc,ord='ca',dim=5,perm=500)
par(mar=c(1,1,1,1))     
evplot(spe.ca$CA$eig)      

 #We can look at this same information graphically with the ordi.scree() function:
ordi.scree(spe.ca,ord='ca')      

spe.ca$CA$u[,1:2] # sample scores (‘u’)
spe.ca$CA$v[,1:2] # species scores (‘v’)   
par(mar=c(1,1,1,1))
#Here, let’s create a CA joint-plot using scaling type 2 by typing:
  ordiplot(spe.ca,choices=c(1,2),scaling=2)
#You might be interested in labeling the plot above with object (site) names to help with interpretation.
  ordiplot(spe.ca,type='t', scaling=2)
#You might also interested in comparing the two scaling approaches. To do this you can type:  
  par(mfrow=c(1,2))
  ordiplot(spe.ca, scaling=1) 
  ordiplot(spe.ca, scaling=2)  
  
#DactoMineR plot
  
#First type the following to conduct the CA:
    spe.ca2<-CA(speocc)  

#ype the following to obtain the other results
    summary(spe.ca2)        
    