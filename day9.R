#day 9
#BACKGROUND
# When the effects of a particular predictor dataset need to be tested after elimination of
# possible effects due to other variables, partial ordination may be used (e.g. partial CCA, 
#                                                                          partial RDA). 
# Such an approach is also referred to as ‘partialling out’ or ‘controlling for’ the effects
# of specific variables, which are specified as covariables in the constrained analysis. 
# In contrast to the previous exercise in which a single set of explanatory variables (or constraints)
# were used in the constrained ordination, here we are going to subdivide or partition the explanatory 
# variables into logical subsets of variables and determine how much of the species variance can be 
# explained independently by each subset or jointly by two or more subsets. Specifically, we seek to 
# understand the independent and joint (or confounded) explanatory power of the logical subsets of 
# variables.
# The variance decompositions that result from canonical partitioning provide a comprehensive picture 
# of the relative importance, independent effects, and confounding of the factors included in the
# analysis. The canonical partitioning method was originally used to partition variation in community
# data sets among environmental and spatial components (Borcard et al. 1992). Anderson and Gribble (1998)
# extended the technique to include the effects of temporal variation, so that variation in community data 
# is partitioned into the components explainable by environmental, spatial, and temporal factors, and their
# overlap. Cushman and McGarigal (2002 and 2004) extended the method to address the specific challenges of 
# hierarchically structured data.
# Importantly, the canonical partitioning approach can be applied to a variety of ordination procedures,
# including the linear methods of redundancy analysis (RDA) and distance-based RDA (dbRDA) 
# (also referred to as constrained analysis of principal coordinates, CAP), 
# and the unimodal method of canonical correspondence analysis (CCA).
# Before diving into the analyses, let’s explore the general concept of partial CCA (or RDA or CAP)
# using a simple example. In doing so you will quickly see that partial direct ordination operates in a 
# fashion similar to partial regression. In this example we have
# the species response matrix Y and two predictor
# matrices X1 and X2. First, both sets of predictor matrices (X1 and X2) are included in a 
# constrained ordination. The amount of variance explained (total inertia for CCA) corresponds to the 
# total variance explained by X1 and X2. Next, a partial CCA is performed where the effect of X1 
# is cleansed from the species matrix Y. CCA is performed using only X1 and residuals from the analysis
# are submitted to a second CCA using X2 as the constraining matrix. Effectively, the effect of X1 is 
# removed from the species community matrix Y and the unique proportion of
# variance explained by X2 can be determined. The process is repeated a second time, but now the 
# species community data set is cleansed of the effects of X2 in order to isolate the variance explained by X1.
# The figure above conceptualizes the process of variance partitioning, 
# where X1X2 is the proportion of variance explained by both matrices X1 and X2.
#load multiple pkgs
lapply(c("vegan"), require, character.only = TRUE)

#source
source('biostats.R')

#Let’s import the data, by typing:
speabu <-read.csv('./data/MAHA_speciesabu.csv',header=TRUE, row.names=1)
envdata <-read.csv('./data/MAHA_environment.csv',header=TRUE, row.names=1)


# Now let’s transform both datasets before diving into the analysis. First, let’s conduct a ‘row normalization’ 
# (i.e., to rescale each row so that the sums equal 1, or in other words to calculate relative species abundances), type the following:
  speabu.tran<-data.stand(speabu,method='total',margin='row',plot=F)
# Next, let’s log (based 10) transform the environmental dataset.
envdata.tran<-data.trans(envdata,method='log',plot=F)

# PARTIAL CANONICAL CORRESPONDENCE ANALYSIS
# Here, we will examine the independent and joint variation in species composition explained by 
# environmental descriptors that depict natural habitat processes versus human-induced impacts on the river.
# Natural questions to address are whether natural and human variables explain unique or shared variance 
# in the fish community data and whether one set of variables explains more variance than the other. 
# To help with the analysis below, let’s divide the environmental variables into 2 matrices by combining 
# columns from envdata.tran using the function cbind(). Do this by typing:
  attach(envdata.tran)
env.human<-cbind(HabQual, RoadDen,Agricult,HumanUse) 
env.natural<-cbind(Sinuosity,Slope,WDRatio,SubEmbed,Elev,BasinAre) 
detach(envdata.tran)

# Conduct the variance partitioning
# In this lab we will use two different ways of conducting pCCA in R – both of which utilize the cca 
# function and thus result in the same output. The cca() function in the vegan library allows for both
# unconstrained correspondence analysis (CA) as well as the constrained version, canonical 
# correspondence analysis (CCA). Recall that we used this function in the previous chapter. Its usage is:
#   cca (X, Y, Z,...)

spe.pcca1<-cca(speabu.tran,env.human,env.natural) 
summary(spe.pcca1)

# Let’s repeat the same variance partitioning, but this time we will isolate the variation in 
# species composition explained uniquely by naturals factors. Notice that we simply switched the Y and Z 
# matrices! Do this by typing:
  spe.pcca2<-cca(speabu.tran,env.natural,env.human) 
  summary(spe.pcca2)
  
  # Perhaps a better way of conducting the above partitioning is to utilize the BIOSTATS library. 
  # To compute the variance partitioning, type:
spe.pcca<-ordi.part(speabu.tran,env.natural,env.human,method='cca')

# Venn Diagrams
# The entire variance partitioning can be graphically portrayed in a Venn diagram using the plot.ordi.part() 
# function in biostats, by typing:
  plot.ordi.part(spe.pcca,which='total')    

  #Variance explained after controlling for env facyors
  
  spe.pcca_ce<-ordi.part(speabu.tran,env.human,method='cca',p=env.natural)
  