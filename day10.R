#day 10

lapply(c("vegan",'MASS', 'caret', 'e1071'), require, character.only = TRUE)

#source
source('biostats.R')

#Let’s import the data, by typing:
speabu <-read.csv('./data/MAHA_speciesabu.csv',header=TRUE, row.names=1)
envdata <-read.csv('./data/MAHA_environment.csv',header=TRUE, row.names=1)

sitegroup <- read.csv('./data/MAHA_groups.csv',header=TRUE, row.names=1)
# Recall that sitegroup contains two column vectors depicting the major basin 
# and degree of human impact for each site (object).
# Let’s log (based 10) transform the environmental dataset.
envdata.tran<-data.trans(envdata,method='log',plot=F)

#Now, let’s perform LDA by typing:
  env.lda<-lda(envdata.tran, sitegroup[,1])
#To return the results from the analysis, type:
  env.lda

#Canonical Scores
# We can compute the canonical scores, i.e., the scores for
#  each sample observation on each canonical axis, using the predict() function, as follows:
env.lda.pred<-predict(env.lda) 
env.lda.pred


# Note, the resulting object contains the predicted group membership (which we will
# use below to evaluate the classification criterion), the posterior probability of membership in 
# each group (which is the basis for the prediction of group membership), and the canonical scores. 
# For convenience, we can extract the scores from this object (the list component named ‘x’) and 
# store them in a separate object, by typing:
  scores<-env.lda.pred$x
 scores 

 par(mar=c(1,1,1,1))
 result<-as.data.frame(cbind(sitegroup[,1],scores))
 #Then we can produce grouped box-and-whisker plots by typing:
 box.plots(result,var='LD1',by='V1',notch=TRUE,varwidth=TRUE) 
 
 env.table<-table(sitegroup[,1],env.lda.pred$class) 
 env.table
 #MidAt Ohio MidAt 21 2 Ohio 2 9 SouthAt 4 0 Tenness 2 0
 #SouthAt Tenness 0 0 0 0 2 0 0 3
 
 # The rows of the table are the ‘from’ group (i.e. the true or known group of each observation) 
 # and the columns are the ‘to’ group (i.e., the predicted group membership). 
 # The diagonals of the table give the correct predictions (frequencies) and the off-diagonals 
 # give the errors of omission and commission. To compute the overall correct classification rate (CCR), type:
   sum(diag(env.table))/sum(env.table)
   # This shows that the LDA model correctly classifies sites
   # Better yet, let’s calculate a whole suite of indices based on this table 
   # (which is called a confusion matrix), including measures of CCR, sensitivity, 
   # specificity, Cohen’s Kappa index, etc., by typing:
     confusionMatrix(env.table)   
     # Interpreting the Canonical Functions
     # Assuming that we achieve good group separation on the canonical axes, 
     # our next task is to determine if we can provide a meaningful ecological interpretation
     # of the group differences. First, we might suspect that the eigenvector coefficients or
     # the canonical coefficients of the linear discriminant functions might provide the answer. 
     # We can view the raw canonical coefficients of the linear functions by typing:
       env.lda$scaling     
#       These structure coefficients can be computed using the lda.structure() function in biostats, by typing:
         lda.structure(scores,envdata.tran)
         # LDA Bioplot
         # Finally, let’s visualize the ordination of the canonical axes from the LDA.
         par(mar=c(1,5,5,1))         
          plot(scores,type='n') 
          text(scores,labels=sitegroup[,1])       
          