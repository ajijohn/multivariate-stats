#day11 - MRPP

# Multi-response permutation procedure (MRPP)is a nonparametric approach for
# testing the hypothesis of no difference between two or more groupsof entities
# (Mielke et al. 1976). The strategy of MRPP is to comparethe observed
# intra-group average distances with the averagedistances that would have
# resulted from all the other possible combinations of the data under thenull
# hypothesis. The test statistic, usually symbolized with a lower case delta, is
# the average ofthe observed intra-group distances weighted by relative group
# size. The observed delta iscompared to the possible deltas resulting from
# Monte Carlo permutations of the assignment ofsample observations to the
# groups. If the hypothesis that the groupsare not different (the
# nullhypothesis) is true, then each of the possible assignments (permutations)
# is equally likely. Theflexibility of MRPP stems from the use of any distance
# measure.

lapply(c("vegan"), require, character.only = TRUE)

#source
source('biostats.R')

#Let’s import the data, by typing:
speabu <-read.csv('./data/MAHA_speciesabu.csv',header=TRUE, row.names=1)
envdata <-read.csv('./data/MAHA_environment.csv',header=TRUE, row.names=1)

sitegroup <- read.csv('./data/MAHA_groups.csv',header=TRUE, row.names=1)

# Calculate dissimilar matrixThe firststep in MRPP is the selectionof a distance
# matrix to describe similarities among objects according to the descriptors.
# Again, any (dis) similarity coefficient can be selected.  For today, let’s
# generate a distance matrix using an appropriate dissimilarity coefficientfor
# speciesabundance –Bray-Curtis coefficient. We’ll use the vegdist() function
# (see Ecological Resemblace chapter). Type,
speabu.d<-vegdist(speabu, "bray")

sp.mrpp<-mrpp(speabu.d,sitegroup[,1])
sp.mrpp

# The resulting object is a list containing several components, butthe most
# important information is the value of the test statistic, delta, the p-value
# derived from the Monte Carlo permutation test,and the value of the ‘effect
# size’, A.From the above output we can see that the river basins are
# significantly different in their fish community composition. If several groups
# are included in MRPP, subsequent paired tests can be performed to identify
# which groups in particular significantly differ.Now type:
sp.md <-meandist(speabu.d, sitegroup[,1])
sp.md

summary(sp.md)

#We can plot the distribution of random delta values by typing:
hist(sp.mrpp$boot.deltas)


spe.nmds<-metaMDS(speabu, distance='bray', k=2, autotransform=FALSE,trymax=100)
plot(spe.nmds,type='n')
text(spe.nmds,labels=sitegroup[,1])
ordihull(spe.nmds,sitegroup[,1],show.groups="MidAt",col=c("red"))
ordihull(spe.nmds,sitegroup[,1],show.groups="SouthAt",col=c("blue"))
ordihull(spe.nmds,sitegroup[,1],show.groups="Ohio",col=c("green"))
ordihull(spe.nmds,sitegroup[,1],show.groups="Tenness",col=c("purple"))


#Finally, to produce a plot without different group colors simply type:
plot(spe.nmds,type='n')
text(spe.nmds,labels=sitegroup[,1])
ordihull(spe.nmds,sitegroup[,1])


