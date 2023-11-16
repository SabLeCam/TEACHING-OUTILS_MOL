# TP Analyses de différenciation entre populations

( France Dufresne et Sabrina Le Cam, d'après un document produit par Tom Jenkins source https://tomjenkins.netlify.app/tutorials/r-popgen-getting-started/)

## Ressources et packages nécessaires

Nous allons effectuer ces analyses sous R. 
Installez les packages suivants:

```r
source("http://bioconductor.org/biocLite.R")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.16")
`%>%` <- magrittr::`%>%`

source("http://membres-timc.imag.fr/Olivier.Francois/POPSutilities.R")

source("http://membres-timc.imag.fr/Olivier.Francois/Conversion.R")

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("LEA")

install.packages("tidyverse")
```


Ouverture des librairies
```r
library("adegenet")
library("poppr")
library("dplyr")
library("hierfstat")
library("reshape2")
library("ggplot2")
library("RColorBrewer")
library("scales")
library("LEA")
library("tidyverse")
```

Nous repartons du jeux de données de l'étude de [Jenkins et al. 2019](https://onlinelibrary.wiley.com/doi/10.1111/eva.12849) sur le Homard bleu avec lequel nous avons fait le [TP Hardy Weinberg](https://github.com/SabLeCam/OUTILS_MOL/tree/main/HardyWeinberg)
Récupérer les données en suivant le début du [TPHW](https://github.com/SabLeCam/OUTILS_MOL/blob/main/HardyWeinberg/readme.md) jusqu'à l'étape suivante:

```r
lobster_gen
```
```r
lobster_gen_sub = popsub(lobster_gen, sublist = c("Ale","Ber","Brd","Pad","Sar17","Vig"))
```
```r
 dataset.hfstat <- genind2hierfstat(lobster_gen_sub)
```
On compare maintenant les patrons de diversité génétique entre les populations par rapport à la diversité globale.

## Matrice de Fst par paire de population (estimateur du Fst de Weir de Cockerham(1984))

```r
lobster_fst = genet.dist(lobster_gen_sub, method = "WC84") %>% round(digits = 3)
lobster_fst
```
Représenter cette matice avec une heatmap

```r
mat.obs <- pairwise.WCfst(dataset.hfstat)
#garder que la partie supérieure de la matrice
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }
mat.obs_upper <-get_upper_tri(mat.obs)
```
```r
melted_matobs <- melt(mat.obs_upper, na.rm = TRUE)
colnames(melted_matobs)<-c("pop1","pop2","value")

ggheatmap <- ggplot2::ggplot(melted_matobs, aes(pop1, pop2, fill = value)) +
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red",  
                       midpoint = 0.075, limit = c(0,0.15), space = "Lab" ) +
  geom_text(aes(label = round(value,digits = 3)), color="black", size = 3)+
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))
ggheatmap
```

<p float="left">
 <img width="450" alt="image" src="https://github.com/SabLeCam/OUTILS_MOL/assets/20643860/24783241-4586-465d-b6fc-d7c0c942f514">
 <img width="450" alt="image" src="https://github.com/SabLeCam/OUTILS_MOL/assets/20643860/8710515c-b3ac-482f-a86f-84e643239c01.png">
</p>

#Quelles populations possèdent les valeurs de Fst  les plus divergentes ?
Quelles populations possèdent les valeurs de Fst les plus faibles ?

## Analyses ACP (analyse en composantes principales)

On replace les données manquantes par la moyenne des fréquences allèliques
```r
x = tab(lobster_gen_sub, NA.method = "mean")
```

On réalise l'ACP
```r
pca<-dudi.pca(x, scannf = FALSE, scale = FALSE, nf = 3)
```
En représentant les vecteurs eigen, on peut analyser quel pourcentage de la variance génétique est expliqué par chacun des axes
```r
percent = pca$eig/sum(pca$eig)*100
barplot(percent, ylab = "Genetic variance explained by eigenvectors (%)", ylim = c(0,12),
        names.arg = round(percent, 1))
```
![image](https://github.com/SabLeCam/OUTILS_MOL/assets/20643860/826543ae-de37-4f55-a420-10eeed27a2f7)

On représente graphiquement les résultats de l'ACP
```r
dfpca<-data.frame(a1=pca$li[,1],a2=pca$li[,2],pop=as.vector(lobster_gen_sub$pop),ind=rownames(pca$li))

centroids <- aggregate(cbind(a1,a2)~pop,data=dfpca,mean)

gpca<-ggplot(data = dfpca, aes(x=a1,y=a2, color=pop)) +
  theme(panel.background=element_blank(),
        axis.text=element_text(family="Arial Narrow", size=14,color='grey40'),
        axis.title=element_text(family="Arial Narrow", size=14,color='grey40'),
        panel.border=element_rect(fill=NA,colour="grey40")) +
  xlab(paste0('Axis',1))+
  ylab(paste0('Axis',2))+
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") +
  geom_point( position= "jitter", alpha = 0.8, size = 4) +
  geom_label(data=centroids,aes(label=pop), size=5,position="jitter",
             fill= 'white',show.legend = FALSE)

gpca
```
#Comment peut-on décrire les résultats de cette analyse ?
>On peut refaire ce graphique en considérant les composantes 1 et 3. qu'est ce que ce nouveau graphique nous indique?

![image](https://github.com/SabLeCam/OUTILS_MOL/assets/20643860/bc254c2c-df8f-4263-98af-4cbb1c264ba3)
![image](https://github.com/SabLeCam/OUTILS_MOL/assets/20643860/5729add6-54fa-4be8-9401-e93e31c7962c)




## Inférence bayésienne de la structure de population

Formater les données en format 'geno' supporté par le package LEA
```r
# extraire les données allèliques
allelecount <- as.data.frame(tab(lobster_gen_sub))
# remplacer les données manquante par 9
allelecount[is.na(allelecount)] <- 9
# écrire un fichier en format lfmm
write.lfmm(allelecount, "lobster.lfmm")
#écrire un fichier en format geno
geno <- lfmm2geno(input.file = "lobster.lfmm",
                  output.file = "lobster.geno")
```
### Faire l'analyse d'inférence bayésienne
calculer coefficient de coancestrie, c'est à dire la probablilité de tirer deux gènes identiques dans deux individus d'une même population. (nombre de pop testées (K) entre 1 et 10 et 10 répétition par K)
```r
snmf.obj1 <- snmf(input.file = "lobster.geno",
                  K = 1:10,
                  project = "new",
                  repetitions = 10,
                  CPU = 65,
                  entropy = TRUE,
                  seed = 42)
```
Identifier le K qui est statistiquement le plus supporté à partire du criètre de coss-entropie pour chaque K
représentation graphique

```r
K <- summary(snmf.obj1)$crossEntropy %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("temp") %>%
  mutate(K = as.numeric(str_extract(temp, "(\\d)+"))) %>%
  select(-temp)

ggplot(K, aes(x = K, y = mean)) +
  geom_line(color = "black", size = 0.25 ) +
  geom_segment(aes(x = K, y = min, xend = K, yend = max)) +
  geom_point(shape = 21, size = 4, color = "black", fill = "darkorange") +
  scale_x_continuous(breaks = seq(0, 10, by = 1)) +
  labs(x = "number of ancestral populations", y = "cross-entropy criterion")
```
Choisir le K pour lequel la fonction atteint un plateau ou présente une forte augmentation (K=3)

![image](https://github.com/SabLeCam/OUTILS_MOL/assets/20643860/9c6e0cae-3acd-42db-84ec-cd8228888a01)

Identifier le meilleur réplicat pour le K donné
```r
run <- which.min(cross.entropy(snmf.obj1, K = 3))
run
```

Extraire la matrice de valeur q et faire un graphique
```r
qmatrix = Q(snmf.obj1, K = 3, run=1)
```
```r
par(mar=c(4,4,0.5,0.5))
pops<-levels(lobster_gen_sub@pop)
barplot(t(qmatrix), col=RColorBrewer::brewer.pal(9,"Paired"), 
        border=NA, space=0, xlab="Individuals", 
        ylab="Admixture coefficients")
#Add population labels to the axis:
for (i in 1:length(pops)){
  axis(1, at=median(which(lobster_gen_sub@pop==pops[i])), labels=pops[i])}
```
ajouter les délimitations des populations
```r
table(lobster_gen_sub$pop)
abline(v=c(35,68,104,140,155))
```
![image](https://github.com/SabLeCam/OUTILS_MOL/assets/20643860/98084149-ee1c-4c2d-b233-c158a48bdc03)
>Vérifier les résulats avec K=2, qu'est que ça nous apprend?
![image](https://github.com/SabLeCam/OUTILS_MOL/assets/20643860/d8c8d144-ae0b-4fe9-ab7c-833e5814d548)
>Quels facteurs peuvent expliquer cette structure de populations ?
>Essayez maintenant de refaire ces analyses avec le jeux de données de gorgones.

seafan_gen_sub = popsub(seafan_gen, sublist = c("ArmI","Far","Bla","Thu","Bre","Lio"))


```r
####gorgones pop diff #####


seafan_gen = import2genind("Pinkseafan_13MicrosatLoci.gen", ncode = 3, quiet = TRUE)
set.seed(1)
tab(seafan_gen[loc = "Ever002"])[runif(3, 1, nInd(seafan_gen)), ]
popNames(seafan_gen) = gsub("[^a-zA-Z]", "", popNames(seafan_gen))
popNames(seafan_gen)



#ici c'est les pop selectionnée par Jenkins
seafan_gen_sub = popsub(seafan_gen, sublist = c("Bla","Bov","Bre","Lun","PorI","Sko"))
dataset.hfstat <- genind2hierfstat(seafan_gen_sub)


###matrice de Fst
mat.obs <- pairwise.WCfst(dataset.hfstat)

get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}



mat.obs_upper <-get_upper_tri(mat.obs)
melted_matobs <- melt(mat.obs_upper, na.rm = TRUE)
colnames(melted_matobs)<-c("pop1","pop2","value")

ggheatmap <- ggplot2::ggplot(melted_matobs, aes(pop1, pop2, fill = value)) +
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red",  
                       midpoint = 0.05, limit = c(0,0.12), space = "Lab" ) +
  geom_text(aes(label = round(value,digits = 3)), color="black", size = 3)+
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))
ggheatmap


###ACP###
x = tab(seafan_gen_sub, NA.method = "mean")
pca<-dudi.pca(x, scannf = FALSE, scale = FALSE, nf = 3)

percent = pca$eig/sum(pca$eig)*100
barplot(percent, ylab = "Genetic variance explained by eigenvectors (%)", ylim = c(0,12),
        names.arg = round(percent, 1))


dfpca<-data.frame(a1=pca$li[,1],a2=pca$li[,3],pop=as.vector(seafan_gen_sub$pop),ind=rownames(pca$li))

centroids <- aggregate(cbind(a1,a2)~pop,data=dfpca,mean)

gpca<-ggplot(data = dfpca, aes(x=a1,y=a2, color=pop)) +
  theme(panel.background=element_blank(),
        axis.text=element_text(family="Arial Narrow", size=14,color='grey40'),
        axis.title=element_text(family="Arial Narrow", size=14,color='grey40'),
        panel.border=element_rect(fill=NA,colour="grey40")) +
  xlab(paste0('Axis',1))+
  ylab(paste0('Axis',3))+
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") +
  geom_point( position= "jitter", alpha = 0.8, size = 4) +
  geom_label(data=centroids,aes(label=pop), size=5,position="jitter",
             fill= 'white',show.legend = FALSE)

gpca

##STRUCTURE##

#input les données, fichier avec hearder différent
# FORMAT DATA ----
# extract allele counts
allelecount <- as.data.frame(tab(seafan_gen_sub))
# replace missing values with 9
allelecount[is.na(allelecount)] <- 9
# write lfmm input file
write.lfmm(allelecount, "seafan.lfmm")
geno <- lfmm2geno(input.file = "seafan.lfmm",
                  output.file = "seafan.geno")
# RUN ANALYSIS ----
# calculate coancestry coefficents (K = 1-10 for 10 iterations)
snmf.obj1 <- snmf(input.file = "seafan.geno",
                  K = 1:10,
                  project = "new",
                  repetitions = 10,
                  CPU = 65,
                  entropy = TRUE,
                  seed = 42)
# IDENTIFY STATISTICALLY MOST SUPPORTED K ----
# plot cross-entroypy criterion for each K
K <- summary(snmf.obj1)$crossEntropy %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("temp") %>%
  mutate(K = as.numeric(str_extract(temp, "(\\d)+"))) %>%
  select(-temp)
# choose K for which the function plateaus or increases sharply (K = 2)
ggplot(K, aes(x = K, y = mean)) +
  geom_line(color = "black", size = 0.25 ) +
  geom_segment(aes(x = K, y = min, xend = K, yend = max)) +
  geom_point(shape = 21, size = 4, color = "black", fill = "darkorange") +
  scale_x_continuous(breaks = seq(0, 10, by = 1)) +
  labs(x = "number of ancestral populations", y = "cross-entropy criterion")

# identify best run for chosen K value
run <- which.min(cross.entropy(snmf.obj1, K = 3))

# extract q.matrix

qmatrix = Q(snmf.obj1, K = 3, run=6)
barplot(t(qmatrix), col = c("orange","violet","lightgreen","gray","red","darkblue"), border = NA, space = 0,
        xlab = "Individuals", ylab = "Admixture coefficients")

par(mar=c(4,4,0.5,0.5))
pops<-levels(seafan_gen_sub@pop)
barplot(t(qmatrix), col=RColorBrewer::brewer.pal(9,"Paired"), 
        border=NA, space=0, xlab="Individuals", 
        ylab="Admixture coefficients")
#Add population labels to the axis:
for (i in 1:length(pops)){
  axis(1, at=median(which(seafan_gen_sub@pop==pops[i])), labels=pops[i])}

table(seafan_gen_sub$pop)

abline(v=c(29,49,92,114,153))
```
## DAPC 
```r
# Perform cross validation to find the optimal number of PCs to retain in DAPC
set.seed(123)
x = tab(seafan_gen_sub, NA.method = "mean")
crossval = xvalDapc(x, seafan_gen_sub$pop, result = "groupMean", xval.plot = TRUE)

# Number of PCs with best stats (lower score = better)
crossval$`Root Mean Squared Error by Number of PCs of PCA`
##        10        20        30        40        50        60        70 
## 0.6252777 0.6326131 0.6380681 0.6057849 0.6587395 0.6412447 0.6113320
crossval$`Number of PCs Achieving Highest Mean Success`
## [1] "40"
crossval$`Number of PCs Achieving Lowest MSE`
## [1] "40"
numPCs = as.numeric(crossval$`Number of PCs Achieving Lowest MSE`)

# Run a DAPC using site IDs as priors
dapc1 = dapc(seafan_gen_sub, seafan_gen_sub$pop, n.pca = numPCs, n.da = 3)

# Analyse how much percent of genetic variance is explained by each axis
percent = dapc1$eig/sum(dapc1$eig)*100
barplot(percent, ylab = "Genetic variance explained by eigenvectors (%)", ylim = c(0,60),
        names.arg = round(percent, 1))


# Create a data.frame containing individual coordinates
ind_coords = as.data.frame(dapc1$ind.coord)

# Rename columns of dataframe
colnames(ind_coords) = c("Axis1","Axis2","Axis3")

# Add a column containing individuals
ind_coords$Ind = indNames(seafan_gen_sub)

# Add a column with the site IDs
ind_coords$Site = seafan_gen_sub$pop

# Calculate centroid (average) position for each population
centroid = aggregate(cbind(Axis1, Axis2, Axis3) ~ Site, data = ind_coords, FUN = mean)

# Add centroid coordinates to ind_coords dataframe
ind_coords = left_join(ind_coords, centroid, by = "Site", suffix = c("",".cen"))

# Define colour palette
cols = brewer.pal(nPop(seafan_gen_sub), "Set2")

# Custom x and y labels
xlab = paste("Axis 1 (", format(round(percent[1], 1), nsmall=1)," %)", sep="")
ylab = paste("Axis 2 (", format(round(percent[2], 1), nsmall=1)," %)", sep="")

# Scatter plot axis 1 vs. 2
ggplot(data = ind_coords, aes(x = Axis1, y = Axis2))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  # spider segments
  geom_segment(aes(xend = Axis1.cen, yend = Axis2.cen, colour = Site), show.legend = FALSE)+
  # points
  geom_point(aes(fill = Site), shape = 21, size = 3, show.legend = FALSE)+
  # centroids
  geom_label(data = centroid, aes(label = Site, fill = Site), size = 4, show.legend = FALSE)+
  # colouring
  scale_fill_manual(values = cols)+
  scale_colour_manual(values = cols)+
  # custom labels
  labs(x = xlab, y = ylab)+
  ggtitle("Pink sea fan DAPC")
```









