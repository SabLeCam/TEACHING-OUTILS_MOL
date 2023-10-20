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

li<-pca$scores
dfpca<-data.frame(a1=li[,1],a2=li[,2],pop=lobster_gen_sub$pop,ind=rownames(pca$scores))

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
>Vérifier les résulats avec K=2, qu'est que ça nous apprend?>
![image](https://github.com/SabLeCam/OUTILS_MOL/assets/20643860/d8c8d144-ae0b-4fe9-ab7c-833e5814d548)

>Essayez maintenant de refaire ces analyses avec le jeux de données de gorgones.
>Travaillez sur une sous-echantillons comprenant les populations ("Bla","Bov","Bre","Lun","PorI","Sko")







