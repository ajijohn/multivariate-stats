#load multiple pkgs
lapply(c("vegan", "cluster", "pvclust", "NbClust", "clusteval"), require, character.only = TRUE)

#source
source('biostats.R')
envdata <- read.csv('./data/MAHA_environment.csv',header=TRUE, row.names=1) 
speabu <- read.csv('./data/MAHA_speciesabu.csv',header=TRUE, row.names=1)

#1. Compute an appropriate dissimilarity/distance matrix

#
#distance matrix

envdata.tra<-data.stand(envdata,method='standardize',margin='column',plot=F)

site.eucd<-vegdist(envdata.tra,method='euclidean')

#compute 2. Compute hierarchical clustering

sitecl.ave<-hclust(site.eucd,method='average') 
sitecl.cen<-hclust(site.eucd,method='centroid') 
sitecl.med<-hclust(site.eucd,method='median') 
sitecl.sin<-hclust(site.eucd,method='single') 
sitecl.com<-hclust(site.eucd,method='complete') 
sitecl.ward<-hclust(site.eucd,method='ward.D2')


#examining clustering results
hclus.table(sitecl.ave)

#plot the dendrogrsam

plot(sitecl.ave)
#fancy version
plot(sitecl.ave,main='Average-linkage Dendrogram',xlab='Sites',ylab='Euclidean distance',hang=-1)

speocc <- data.trans(speabu,method='power',exp=0,plot=F) 
sp.jacd<-vegdist(speocc,method="jaccard") 
specl.ave<-hclust(sp.jacd,method='average') 
plot(specl.ave)

plot(sitecl.ave) 
rect.hclust(sitecl.ave,k=9)

#identifies 9 clusters in the dendrogram.
plot(sitecl.ave) 
rect.hclust(sitecl.ave,h=6)

#“cut” the dendrogram into clusters, by typing:
sitecl.class<-cutree(sitecl.ave,k=9)

#4. Cluster Diagnostics: Evaluating the cluster solution
coef.hclust(sitecl.ave)

#The cophenetic correlation can be computed directly by typing:
  cor(site.eucd,cophenetic(sitecl.ave))
  
hclus.cophenetic(site.eucd, sitecl.ave)  

#5. Cluster Diagnostics: Determining the number and stability of clusters
hclus.scree(sitecl.sin)
hclus.scree(sitecl.ave)  
