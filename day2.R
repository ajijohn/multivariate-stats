envdata <-read.csv('./data/MAHA_environment.csv',header=TRUE, row.names=1)
source('./biostats.R')


speabu <-read.csv('./data/MAHA_speciesabu.csv',header=TRUE, row.names=1)
spetrait <-read.csv('./data/MAHA_speciestrait.csv',header=TRUE, row.names=1)

#convert species abundance to species occurence
speocc <-data.trans(speabu,method='power',exp=0,plot=F)

#similarity matrix
library(simba)
sp.jac <-sim(speocc, method = "jaccard")

sp.sim <-sim(speocc, method = "simplematching")
sp.sor <-sim(speocc, method = "soerensen")

#Simple Matchinh - uses double zerors
#Jaccards - excludes 


#difference in coefficients (pairwise coefficients)
plot(sp.jac,sp.sim,xlab="Jaccard’s coefficient",ylab="Simple Matching coefficient")

abline(0,1,col="darkgray")


plot(sp.jac[1:45],sp.sim[1:45],xlab="Jaccard’s coefficient",ylab="Simple Matching coefficient",type="n")
text(sp.jac[1:45],sp.sim[1:45],row.names(speocc))
abline(0,1,col="darkgray")

#comparison with Sor

#why sorr is (both don;t look double zeroes)
# double weighingting where there us 1 i.e occur\

#difference in coefficients
plot(sp.jac,sp.sor,xlab="Jaccard’s coefficient",ylab="Sor coefficient")

abline(0,1,col="darkgray")

#calculatio=ng coeeficient for continouns data

#
#Let’s try calculating a dissimilarity (distance) matrix based on the Bray-Curtis coefficient:  
sp.bray <-vegdist(speabu, method="bray")
head(sp.bray)

sp.bray

library(gclus)
source('./coldiss.R')
#display
coldiss(sp.bray,nc=4,byrank=FALSE,diag=TRUE)

#does the data needs to be paramtrzed ?


sp.jacd<-vegdist(speocc, method="jaccard")
#You should find thatsp.jacd= 1-sp.jacbecause the only difference is that the 
#first matrix expresses the values as distances and the second matrix contains similarity value

plot(sp.jac,1-sp.jacd)
abline(0,1,col="darkgray")

#There are a number of alternative dissimilarity and distance metrics that are
#appropriate for asymmetrical species data (abundance) that are worth exploring. 
# For example, we can use Chord distanceor Hellinger distanceon 
# the raw species abundance datasets by first transforming the dataset and 
# then calculating the Euclidean distance matrix. 
# This is discussed in detail in the lecture.

#For Chord’s distance:
  speabu.norm<-decostand(speabu, method='nor')
  sp.chordd<-dist(speabu.norm)
#For Hellinger’s distance:
  speabu.hel<-decostand(speabu, method='hel')
  sp.held<-dist(speabu.hel)

#why this plot??
  plot(sp.chordd,sp.held)
  abline(0,1,col="darkgray")
  
#CALCULATING COEFFICIENTS OF SIMILARITY FOR MIXED DATA TYPES
#Use Fowers for mixed
str(spetrait)  
  
#Let’s try calculating a dissimilarity (distance) matrix based on Gower’s coefficient: 
sptr.gower <-daisy(spetrait,metric="gower")  
#looking at similarity if the species
#We can quickly produce a heat map to look at similarities/diferences between 
#species according to their traits. Here is just the ordered heat map.
coldiss(sptr.gower,nc=4,byrank=FALSE,diag=TRUE)
# The resulting plot (below) displays the raw and ordered dissimilarity matrices, 
# where magenta are dissimilarities close to 0 and cyan are dissimilaries close to 1
#cynan-awqua blue

#CALCULATING COEFFICIENTS OF (DIS)SIMILARITY FOR CONTINUOUS DATA
#For continuous data in which double-zeros are meaningful
#(i.e., environment data, species measurements or traits, etc ...) 
#we can calculate dissimilarities between objects based on any number of distance coefficients.
#Let’s calculate environmental dissimilarity 
#among the MAHA sites according to Euclidean distance (the most commonly used coefficient

env.euc <-vegdist(envdata, method="euclidean")
coldiss(env.euc,nc=4,byrank=FALSE,diag=TRUE)

#ONVERTING SIMILARITY TO DISSIMILARITY (DISTANCE) OR VICE VERSAR
#Remember that when dealing with matrices, it is possible to change asimilarity matrix (S) 
#into a dissimilarity matrix (D) byapplying the following transformations:D = 1 –S, D = √(1 –S)

#CONVERTING CORRELATION TO DISTANCE
#Attimes you might want to convert a correlation matrix to a distance matrix 
#in order to perform subsequent multivariate analysis. Of course, 
#you cannot use the transformations listed above because correlation 
#values will vary between -1 and 1, thus resulting in negative distances when R<0.
#Not possible! Therefore I recommend the following transformation:
#D = √(2 –2*correlation value)
#Let’s do this using the environmental data by typing:

env.dis <-sqrt(2-2*cor(envdata))
