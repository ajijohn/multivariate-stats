#day11 - ANOSIM
# Ecologists are often interested in identifying factors that help explain the
# differences between two or more groups of sampling entities. The purpose might
# be to identify ecological or environmental factors that best explain
# differences in the distribution patterns of one or more organisms. Here, the
# emphasis is on explaining the differences among pre-specified, well defined
# groups of sampling entities based on a suite of discriminating variables.
# Alternatively, the purpose might be to predict group membership for samples of
# unknown membership. For example, this might involve predicting whether a
# species is likely to occur at a site based on a suite of ecological or
# environmental variables. Here, the emphasis is on prediction.


# ANOSIM (Analysis of Similarity) is a non-parametric multivariate procedure
# broadly analogous to ANOVA that has been widely used for testing whether or
# not groups of objects are statistically different in respect to their relative
# similarities according to a set of descriptors. It was originally proposed by
# Clarke (1993) to perform community analysis. ANOSIM tests priori-defined
# groups against random groups in ordination space by calculating the average of
# all rank similarities among objects within groups, and the average of rank
# similarities among objects between groups. A test statistic, R, lies between
# -1 and +1. A value of 1 indicates that all objects within groups are more
# similar to one another than any objects from different groups, a value of 0
# indicates that there is no difference among groups (i.e., representing the
# null hypothesis), and a value of -1 indicates that all objects within groups
# are less similar to one another than any objects from different groups (Clarke
# and Gorley 2001).


lapply(c("vegan"), require, character.only = TRUE)

#source
source('biostats.R')

#Let’s import the data, by typing:
speabu <-read.csv('./data/MAHA_speciesabu.csv',header=TRUE, row.names=1)
envdata <-read.csv('./data/MAHA_environment.csv',header=TRUE, row.names=1)

sitegroup <- read.csv('./data/MAHA_groups.csv',header=TRUE, row.names=1)

# ANALYSIS OF SIMILARITY (ANOSIM)
# We can breakdown ANOSIM into a series of steps. Below, we tackle each step in turn.'

# Calculate dissimilar matrix The first step in ANOSIM is the selection of a
# distance matrix to describe similarities among objects according to the
# descriptors. Again, any (dis) similarity coefficient can be selected. For
# today, let’s generate a distance matrix using an appropriate dissimilarity
# coefficient for species abundance – Bray- Curtis coefficient.We’ll use the
# vegdist() function (see Ecological Resemblace chapter).

y.anosim<-anosim(speabu.d,sitegroup[,1])

summary(y.anosim)

# The most important information in the summary is the value of the test
# statistic, R, and the p-value derived from the Monte Carlo permutation test.
# Based on this information, can we reject the null hypothesis of no differences
# in fish community composition among the major river drainages? If so, is the
# difference between groups ecological significant as well based on the
# magnitude of R? To facilitate interpretation, we can plot the dissimilarities
# in a grouped box-and-whisker plot using the plot.anosim() function in
# BIOSTATS, as follows:
plot.anosim(y.anosim)

# It is favorable to use non-metric multidimensional scaling (NMDS) to assist in
# the interpretation of ANOSIM results because both procedures are based on
# ranked similarities. Let’s conduct a NMDS by typing:
spe.nmds<-metaMDS(speabu, distance='bray', k=2, autotransform=FALSE, trymax=100)


plot(spe.nmds,type='n') 
text(spe.nmds,labels=sitegroup$BASIN,col=as.vector(sitegroup$COLOR)) 
ordiellipse(spe.nmds,sitegroup$BASIN,show.groups="MidAt",col=c("red")) 
ordiellipse(spe.nmds,sitegroup$BASIN,show.groups="SouthAt",col=c("blue")) 
ordiellipse(spe.nmds,sitegroup$BASIN,show.groups="Ohio",col=c("green")) 
ordiellipse(spe.nmds,sitegroup$BASIN,show.groups="Tenness",col=c("purple"))

data<-cbind(speabu,sitegroup) 
newdata<-data[which(data$BASIN=='SouthAt'| data$BASIN=='MidAt'), ] 
newdata$BASIN<-factor(newdata$BASIN) 
speabu.dd<-vegdist(newdata[,1:36], "bray") 
yy.anosim<-anosim(speabu.dd,newdata$BASIN)
summary(yy.anosim)


# Lastly, let’s determine which species (variables) best discriminate each of
# the river basins (groups). This provides a way to add biological intpretation
# to the ANOSIM results presented above. To accomplish this we will the simper()
# function which performs pairwise comparisons of groups of objects and finds
# the average contributions of each species to the average overall Bray-Curtis
# dissimilarity. Currently, this function only operates on Bray-Curtis
# dissimilarity, which limits it’s utility. Now type:
simper(speabu,sitegroup$BASIN)

