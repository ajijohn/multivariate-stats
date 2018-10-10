#Day 4

#load multiple pkgs
lapply(c("vegan", "cluster", "pvclust"), require, character.only = TRUE)

#source
source('biostats.R')

#Non-Hierarchical clustering
envdata <- read.csv('./data/MAHA_environment.csv',header=TRUE, row.names=1) 
speabu <- read.csv('./data/MAHA_speciesabu.csv',header=TRUE, row.names=1)

#CONDUCT NON-HIERARCHICAL DIVISIVE CLUSTERING
#where the sole criterion for clustering samples is maximizing within-cluster similarity

#K-MEANS CLUSTER ANALYSIS IN R
# Although most cluster functions in R operate on a dissimilarity/distance matrix, 
# K-means is based directly on the object-by-descriptor dataset. 
# Recall from your lecture that K-means clustering is based on the Euclidean distance among objects.
# However, we still must account for the situations where the descriptors (variables) 
# are measured on different scales. In this example, we need to 
# standardize the environment matrix (mean = 0, variance = 1) in order to control for this fact,
# by typing:
envdata.tra<-data.stand(envdata,method='standardize',margin='column',plot=F)
#Now, let’s calculate the Euclidean distance matrix, by typing:
site.eucd<-vegdist(envdata.tra,method='euclidean')

#We can produce a scree plot of both objective criteria using the nhclus.scree() 
#function in biostats,
#as follows:
nhclus.scree(site.eucd,max.k=44)


# Compute K-means clustering
# The K-means algorithm is a non-hierarchical method to cluster objects based on attributes 
# into K partitions. In other words, given n objects in a p-dimensional space, 
# this algorithm determines a partition of objects into K groups, or clusters, 
# such that the objects within each cluster are more similar to each other than to objects in other
# clusters.  

sitecl.kmeans<-kmeans(envdata.tra,centers=8,iter.max=10000,nstart=25)
#Let’s look at the results. Type:
sitecl.kmeans

# you can obtain the site classifications by typing:
sitecl.kmeans$cluster

#Let’s begin by typing:
sitecl.kmeans.cas<-cascadeKM(envdata.tra,inf.gr=3,sup.gr=10,iter=100) 
plot(sitecl.kmeans.cas,sortg=T)

#Silhoutte plot 
plot(silhouette(sitecl.kmeans$cluster,site.eucd))
#Milligon and Cooper 
#K-means  Strongest evidence for 8 (Use the Calinkski Criterion)

#TODO - Evaluating the number of clusters using randomization
# no.clus<-10
# no.ran<-999
# no.var<-10
# sse<-matrix(,nrow=no.ran+1,ncol=no.clus-1)
# randata<-matrix(,nrow=nrow(envdata.tra),ncol=ncol(envdata.tra))
# for(j in 3:no.clus) 
#   {
#   sitecl.kmeans<-kmeans(envdata.tra,centers=j,iter.max=1000,nstart=5)
#   sse[1,j-1]<-sum(sitecl.kmeans$withinss)
#   #Maximum number of clusters to consider
#   #Number of randomizations to performNumberof descriptors in the dataset
#   For loop that repeats all commands between { and } a total of j times, where j = the number of clusters to consider starting with 3
#   Creates an empty matrix to store within SSE values from the k-means on the original data (first row) and randomized data(rows 2, 3 ... no.ran), where each column is a different number of clustersAlso creates an empty matrix to store the randomized datasetPerforming k-means using the original data, calculating the SSE and storing results in the first row and the j-1 column (e.g., column 1 for a two cluster solution)
#   
#   for(k in 1:no.ran) 
#     {
#     for (c in 1:no.var) {randata[,c]<-envdata.tra[,c][sample(nrow(envdata.tra))] }ran.kmeans<-kmeans(randata,centers=j,iter.max=1000,nstart=5)sse[k+1,j-1]<-sum(ran.kmeans$withinss)}}diff<-colMeans(sse)-sse[1,]plot(seq(2,no.clus,by=1),diff,xlab="No. clusters",ylab="mean(ranSSE)-actualSSE")Based on this information, what level of support do we have for the effectiveness of k-means based on different numbers of clusters? How many clusters would you pick for subsequent interpretation?For loop that repeats all commands between { and } a total of k times, where k = the number of randomizationsPerforming k-means using the random data, calculating the SSE and storing results in the k+1 row and j-1 columnClosing the cluster (j) and randomization (k) for loopsCalculating the difference between the average SSE based on the random k-means and the true SSE based on the actual data. Then the re
#PARTITIONING AROUND MEDOIDS (PAM)
# A number of other non-hierarchical divisive clustering algorithms exist,
# although these approaches are not commonly used in ecology. For examzple, 
# the pam() or clara() functions in the Cluster library are 
# methods for partitioning (clustering) of the data into k clusters ‘around medoids’, 
# which is a more robust version of K-means clustering (Kaufman and Rousseeuw 1990).
# A particularly nice property is that PAM (Partitioning around medoids) allows 
# clustering with respect to any specified distance metric (not just Euclidean distance,
# which k-means is based on). In addition, the medoids are robust representations of the 
# cluster centers, which is particularly important in the common context that 
# many elements donot belong well to any cluster. 
# The clara()function is similar,but is designed for larger data sets with say >200 observations.


sitecl.pam<-pam(site.eucd,k=8) 
sitecl.pam

#the available components (this is true of any object in R). 
#For example, you can obtain the site classifications by typing:
sitecl.pam$cluster
#There is a special plotting function for outputs from the cluster library,
#including pam(). In R, just plot the pam object.
plot(sitecl.pam)

#Interpretation
# Is it good for continuous 
#silhoutte width uses withingroup variation and how far from other groups


# strong affinitues vs the outliers


# Look at the similarity between two classification systemds
# Do homework from this week and past week.

  