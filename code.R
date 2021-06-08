library(R.matlab)
library(FactoMineR)
library(NbClust)
library(aricode)
library(StatMeasures)
library(MASS)
library(ggplot2)
library(factoextra)


setwd("Chemin du workspace")


#Récupération des données
#Nous récupérons l'un des 4 Datasets. Ceci est un exemple avec Data1 de BERT
modelData <- readMat("word_representations__bert-base-cased__Data1.Mat")

#Regroupement des couches
ly = cbind(modelData$layer.01,modelData$layer.02,modelData$layer.03,modelData$layer.04
           ,modelData$layer.05,modelData$layer.06,modelData$layer.07,modelData$layer.08
           ,modelData$layer.09,modelData$layer.10,modelData$layer.11,modelData$layer.12)

data = as.data.frame(ly)
dim(data)




# PCA sur les couches
# Ex. sur la couche 12
i=12
l=eval(parse(text=paste0("modelData$layer.",i)))
res.pcai = PCA(l, graph = FALSE)
fviz_pca_ind(res.pcai, title=paste("Layer",i),geom.ind = "point", 
             pointshape = 21, fill = "#E7B800",repel = TRUE)





# AFM sur toutes les couches
v=c(768,768,768,768,768,768,768,768,768,768,768,768)
names=c("l1","l2","l3","l4","l5","l6","l7","l8","l9","l10","l11","l12")
res = MFA(data, group=v, ncp=12, name.group=names)

fviz_mfa_var(res, title ="AFM - Groupes de variables", "group")
fviz_mfa_var(res, title = "AFM - Cercle de correlation", repel = TRUE)
fviz_mfa_ind(res, title = "AFM - Plan factoriel des Individus", fill = "#E7B800", repel = TRUE)




# PCA sur les 20 premières composantes 
res.pca = PCA(data, ncp=20)
plot(res.pca)




# Application de la CAH avec différents critères et mesurer leurs performances
criteres = c("ward.D", "single", "complete", "average")
realClass <- as.numeric(factor(modelData$real.class, levels = unique(modelData$real.class), exclude = NULL))
ClusterPurity <- function(clusters, classes) {
  sum(apply(table(classes, clusters), 2, max)) / length(clusters)
}

par(mfrow=c(2,2))
for (c in criteres) {
  res<-NbClust(res.pca$ind$coord, distance = "euclidean", min.nc=2, max.nc=6, method = c, index = "ch")
  bestnc = res$Best.nc[1]
  pred = res$Best.partition
  plot(res.pca$ind$coord,col=res$Best.partition, main = c)
  print(c)
  print("confusion table")
  print(table(pred,realClass))
  print("ARI")
  print(ARI(pred,realClass))
  print("NMI")
  print(NMI(pred,realClass))
  print("Purity")
  print(ClusterPurity(pred,realClass))
}





# Application of kmeans and visualisation of clusters
# Et mesurer sa performance
res<-NbClust(res.pca$ind$coord, distance = "euclidean", min.nc=2, max.nc=6, method = "complete", index = "ch")
z.kmeans <- kmeans(res.pca$ind$coord, res$Best.nc[1], nstart = 100)
par(mfrow=c(1,1))
plot(res.pca$ind$coord,col=z.kmeans$cluster,pch="x", xlab = "")

print("confusion matrix")
table(z.kmeans$cluster,realClass)
print("ARI")
print(ARI(z.kmeans$cluster,realClass))
print("NMI")
print(NMI(z.kmeans$cluster,realClass))
print("Purity")
print(ClusterPurity(z.kmeans$cluster,realClass))




# ACP sur le premier plan factoriel
fviz_pca_ind(res.pca, col.ind = modelData$real.class,geom.ind = "point", 
             pointshape = 21,repel = TRUE)



# AFD
d = cbind(data, realClass)
df = as.data.frame(d)
Lda <- lda(realClass~., data = df)
lda.data <- cbind(df, predict(Lda)$x)
p <- ggplot(lda.data, aes(LD1, LD2))+geom_point(aes(color = realClass))
p + ggtitle("AFD")
