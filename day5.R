#day 5

#load multiple pkgs
lapply(c("vegan", "FactoMineR", "factoextra", "dplyr"), require, character.only = TRUE)

#source
source('biostats.R')
envdata <- read.csv('./data/MAHA_environment.csv',header=TRUE, row.names=1) 
speabu <- read.csv('./data/MAHA_speciesabu.csv',header=TRUE, row.names=1)
sitegroup <- read.csv('./data/MAHA_Groups.csv',header=TRUE, row.names=1)

#Let’s perform PCA by typing:
env.pca<-prcomp(envdata, scale=TRUE)
str(env.pca)  

#The results of a PCA are complex and stored in a list, including:
#  • variance (eigenvalue) explained by each eigenvector
#  • variable loadings by eigenvector
#  • sample scores and ordination plots

summary(env.pca)

#Note that this summary presents the standard deviations for each PC and not the variances 
#(i.e., the eigenvalues). Therefore, technically you have to square these values to obtain 
#the eigenvalues. You can do this by calling the values from the env.pca object 
#(recall you can type names(env.pca) to see all available results) by typing:
env.eig<-env.pca$sdev^2 
env.eig

#Next, we can create a screeplot where each eigenvalue for each successive PC is depicted
screeplot(env.pca, bstick=TRUE)


ordi.monte(envdata,ord='pca',dim=5)
# Randomization Test of Eigenvalues:
#   PC1   PC2   PC3   PC4   PC5
# Eigenvalue 2.739 1.895 1.447 1.364 0.795
# P-value    0.000 0.001 0.075 0.002 1.000

#explore the loadings
#Eigenvectors (variable loadings):
#We can examine the eigenvectors (i.e., variable loadings) on each PC (only PC1-5 are presented), 
#by typing:
env.pca$rotation
# PC1        PC2         PC3         PC4         PC5         PC6         PC7
# Sinuosity  0.00429546 -0.4735237 -0.22940933  0.27361829 -0.47648306 -0.54765114  0.20127845
# Slope      0.26319558  0.1303184 -0.34968154 -0.46323791 -0.25803639  0.29083445  0.59479588
# WDRatio    0.37836949 -0.1633701  0.31612993  0.23288139  0.30258611 -0.13536114  0.58552567
# SubEmbed  -0.32966590  0.2478083  0.33469255  0.25747281 -0.39866023  0.36140006  0.29771065
# HabQual    0.35595505 -0.3516462 -0.24667911 -0.05692156 -0.31449155  0.36957851 -0.33295999
# Elev       0.23866884 -0.2113364  0.35123065 -0.53998940  0.19735782 -0.13571147 -0.05421031
# RoadDen   -0.26905123 -0.2427914 -0.49904214  0.15497594  0.52202538  0.23305323  0.16700010
# Agricult  -0.41161178 -0.2821745  0.25449588 -0.38762312 -0.18790521 -0.05668155  0.03804884
# HumanUse  -0.48328053 -0.3592271 -0.02179364 -0.23257076  0.10093767  0.07015571  0.16962229
# BasinAre   0.14080781 -0.4845218  0.33858118  0.25851631  0.00551955  0.49910580 -0.04820080
# PC8          PC9         PC10
# Sinuosity -0.2719200  0.080274347  0.020146213
# Slope     -0.1512735 -0.220912714  0.015587286
# WDRatio    0.4627663  0.094548224  0.012708657
# SubEmbed  -0.1163315  0.508732787  0.032960738
# HabQual    0.4398123  0.381031841 -0.023971395
# Elev      -0.4369820  0.481459480  0.018964299
# RoadDen   -0.1547203  0.280097315  0.372151790
# Agricult   0.3390379 -0.190299509  0.587200145
# HumanUse   0.1585031  0.006536016 -0.716787140
# BasinAre  -0.3565021 -0.430471328  0.009176327

# The loadings indicate the correlation of the original variables with each PC. 
# Strong correlations (either positive or negative) 
# indicate a high contribution towards the linear combination comprising each PC.


# Alternatively, we can use the pca.eigenvec() function as follows:
#   This function suppresses small values below the specified cutoff value (the default is 0, 
# so all coefficients are printed) 
# to emphasize the more important ones.
pca.eigenvec(env.pca,dim=5,digits=3,cutoff=.1)
  # interpretable to convert the 
  # eigenvector coefficients into simple correlation coefficients
  # (i.e., Pearson product-moment correlations), as follows:
pca.structure(env.pca,envdata,dim=5,cutoff=.4)

# Sample scores and ordination plots:
#   Each sample has a score or location on each principal component axis. 
# To see these sample scores, which are stored in the component named ‘x’ in the list object
# output from the prcomp() function, simply type:
  env.scores<-env.pca$x 
  env.scores
  
  #visualize the plot
  
  biplot(env.pca)

  
  #redo the visualization
  
#First type the following to conduct the PCA:
    env2.pca<-PCA(envdata,graph=F)
#Type the following to produce the scree-plot of eigenvalues:
    fviz_screeplot(env2.pca, addlabels = TRUE, ylim = c(0, 30))  
    
#Now, let’s look at the variable loadings (i.e., eigenvectors) by typing:
      fviz_pca_var(env2.pca, col.var = "black")
#and a fancy looking one where total contributions of the variables are color coded:
fviz_pca_var(env2.pca, col.var="contrib",gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE)

#Now let’s look at the relative contributions of the variables to the first PC, by typing:
fviz_contrib(env2.pca, choice = "var", axes = 1, top = 10)

#Producethe PCA biplot by typing:
  fviz_pca_biplot(env2.pca, repel = TRUE)

#Finally, let’s now interpret the ordination plot by symbolizing each site (point) according to their level of human disturbance contained in the sitegroup dataset. We will also draw concentration ellispes. Try typing:
    fviz_pca(env2.pca,label = "var", habillage = sitegroup$DISTURB,palette = c("orange", "green", "yellow","red"), addEllipses = TRUE)
    