#day12-Permanova

# BACKGROUND Another approach to comparing two or groups in multivariate space
# is known as permutational (or nonparametric) multivariate analysis of variance
# (perMANOVA or npMANOVA) (Anderson 2001). perMANOVA is a closely related
# procedure for the analysis and partitioning sums of squares using semi- metric
# and metric distance matrices. Insofar as it partitions sums of squares of a
# multivariate data set, it is directly analogous to MANOVA (multivariate
# analysis of variance). This method allows the use of any semi-metric (e.g.
# Bray-Curtis, aka Steinhaus, Czekanowski, and Sørensen) or metric (e.g.
# Euclidean) distance matrix (Anderson 2001, McArdle and Anderson 2001).
# Significance tests are done using F-tests based on sequential sums of squares
# from permutations of the raw data, and not permutations of residuals.
# perMANOVA is quite flexible and can accommodate a wide variety of experimental
# designs, including nested and factorial designs consisting of both factors and
# continuous variables. If a significant difference between groups is detected
# using perMANOVA, then this could be due to differences in location in
# multivariate space, differences in spread, or a combination of the two. This
# leads us into the next approach. Test of multivariate homogeneity of group
# dispersions can be used to test if the dispersions (variances) of one or more
# groups are different (Anderson 2006). One measure of multivariate dispersion
# (variance) for a group of samples is to calculate the average distance of
# group members to the group centroid or spatial median in multivariate space.
# To test if the dispersions among groups are different, the distances of group
# members to the group centroid are subject to ANOVA. This is a multivariate
# analogue of Levene's test for homogeneity of variances if the distances
# between group members and group centroids is the Euclidean distance. However,
# better measures of distance than the Euclidean distance are available for
# ecological data. These can be accommodated by reducing the distances produced
# using any dissimilarity coefficient to principal coordinates, which embeds
# them within a Euclidean space. The analysis then proceeds by calculating the
# Euclidean distances between group members and the group centroid on the basis
# of the principal coordinate axes rather than the original distances.

lapply(c("vegan"), require, character.only = TRUE)

#source
source('biostats.R')

#Let’s import the data, by typing:
speabu <-read.csv('./data/MAHA_speciesabu.csv',header=TRUE, row.names=1)
envdata <-read.csv('./data/MAHA_environment.csv',header=TRUE, row.names=1)

sitegroup <- read.csv('./data/MAHA_groups.csv',header=TRUE, row.names=1)

spe.perm<-adonis(speabu~BASIN,data=sitegroup,permutations=1000,method='bray') 
spe.perm

hist(spe.perm$f.perms,main="Histogram of pseudo-F values under the null model",xlab="pseudo-F values",xlim=c(0.4,2.0),col='gray')
abline(v=1.89,col="red")

# Now, let’s test for differences in fish community composition among river
# basins and human disturbance according to 4 levels that are classified as low,
# moderate, high and very high. Try typing:
adonis(speabu~BASIN+DISTURB,data=sitegroup,perm=1000,method='bray')


# The results present the 2-way perMANOVA comparing differences in fish
# community composition among river basins and human impact levels. You see that
# we find little evidence for fish community differences between rivers in
# different categories of human disturbance.


# Let’s start by calculating the multivariate dispersions in fish community
# composition between major river basins. The first step is the selection of a
# distance matrix to describe similarities among objects according to the
# descriptors. Again, any (dis) similarity coefficient can be selected. For
# today, let’s generate a distance matrix using an appropriate dissimilarity
# coefficient for species abundance – Bray-Curtis coefficient.We’ll use the
# vegdist() function (see Ecological Resemblace chapter). Type:
speabu.d<-vegdist(speabu,'bray')
#Next, let’s conduct the test of multivariate dispersion by typing:
  sp.bdp<-betadisper(speabu.d,sitegroup$BASIN)
 sp.bdp
 #Next, let’s perform the ANOVA test by typing:
   anova(sp.bdp)

   # You will notice that we find no significant difference in the between-site
   # variation in fish community composition among river basins. We can also
   # calculate significant levels using a permutational approach by typing:
     permutest(sp.bdp,pairwise=TRUE)   
#
#      The ANOVA and permutation test demonstrate little differences between
#      river basins with respect to the degree of variation in fish community
#      composition. Pairwise comparisons are also presented. We can plot a
#      boxplot of distances by typing:
       boxplot(sp.bdp)     

       # Finally, let’s select an appropriate ordination approach to visual our
       # data, namely PCoA. Let’s use the code from Chapter 7 to perform PCoA.
       spe.pcoa<-cmdscale(speabu.d, eig=TRUE, add=T)
       ordiplot(spe.pcoa, choices = c(1, 2), type="text", display="sites", xlab="PC 1", ylab="PC 2")
       ordihull(spe.pcoa,sitegroup$BASIN,show.groups="MidAt",col=c("red")) 
       ordihull(spe.pcoa,sitegroup$BASIN,show.groups="SouthAt",col=c("blue")) 
       ordihull(spe.pcoa,sitegroup$BASIN,show.groups="Ohio",col=c("green")) 
       ordihull(spe.pcoa,sitegroup$BASIN,show.groups="Tenness",col=c("purple"))
              