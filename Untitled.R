
source('biostats.R')
#load multiple pkgs vegan, simba, cluster, ecodist, gclus
lapply(c("vegan", "cluster", "simba", "ecodist", "gclus"), require, character.only = TRUE)
envdata <- read.csv('data/MAHA_environment.csv',header=TRUE, row.names=1)
speabu <- read.csv('data/MAHA_speciesabu.csv',header=TRUE, row.names=1) 
spetrait <- read.csv('data/MAHA_speciestrait.csv',header=TRUE, row.names=1)

# Let’s calculate the similarity among sites according to their fish community composition based 
# on species presence/absence. First, let’s transform the species abundances into presence/absence 
# (i.e., binary transformation) using the power method with an exponent equal to zero, by typing:
speocc <- data.trans(speabu,method='power',exp=0,plot=F)
  
# Remember to include plot=F so that you are not required to press enter for each plot.
#   If you do not you must press enter for all plots before speocc is created. 
#   Alternatively you could directly import the data (assuming that you have prepared it already), by typing:
    
speocc <- read.csv('data/MAHA_speciesocc.csv',header=TRUE, row.names=1)  

# Let’s try calculating Jaccard’s similarity matrix, by typing:
sp.jac <- sim(speocc, method = "jaccard")
#sp.jac
sp.sim <- sim(speocc, method = "simplematching") 
sp.sor <- sim(speocc, method = "soerensen")

# How do these matrices compare to each other? You should be able to note and explain the differences. 
# We can do this graphically by plotting the pair-wise similarities based on two different coefficients 
# against one another, by typing:
plot(sp.jac,sp.sim,xlab="Jaccard’s coefficient",ylab="Simple Matching coefficient")

# higher j , almost higher simple matching

#We can add a 1:1 line to help with the comparisons, by typing:
  abline(0,1,col="darkgray")
  