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

## Matrice de Fst par paire de population (estimateur du Fst de Weir de Cockerham(1984)

```r
lobster_fst = genet.dist(lobster_gen_sub, method = "WC84") %>% round(digits = 3)
lobster_fst
```
Représenter cette matice avec une heatmap

```r
#garder que la partie suppérieure de la matrice
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
 <img width="300" alt="image" src="https://github.com/SabLeCam/OUTILS_MOL/assets/20643860/a44d655f-81bb-4b59-b665-4df44f844f7f.png">
 <img width="450" alt="image" src="https://github.com/SabLeCam/OUTILS_MOL/assets/20643860/8710515c-b3ac-482f-a86f-84e643239c01.png">
</p>

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
![image](https://github.com/SabLeCam/OUTILS_MOL/assets/20643860/bc254c2c-df8f-4263-98af-4cbb1c264ba3)


## Inférence bayésienne de la structure de population

Formater les données en format 'geno' surppoté par le package LEA
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

  
# identify best run for chosen K value
run <- which.min(cross.entropy(snmf.obj1, K = 2))
# format individual names
temp <- as.data.frame(indNames(gen)) %>%
  left_join(strata) %>%
  select(LIB_ID, REGION)
# extract q.matrix
q.matrix <- Q(snmf.obj1, K = 2, run = run) %>% 
  as.data.frame() %>%
  rename(GRP1 = V1,
         GRP2 = V2) %>%
  bind_cols(temp) %>%
  gather(key = GRP, value = MEMBSHIP, 1:2) %>%
  mutate(REGION = ordered(REGION, levels = reg))
  
# plot ancestry coefficients
ggplot(q.matrix, aes(x = LIB_ID, y = MEMBSHIP, fill = GRP, color = GRP)) +
  geom_bar(stat = "identity") +
  labs(x = "INDV", y = "memb. prob") +
  facet_grid(. ~ REGION, scales = "free", space = "free") +
  scale_fill_viridis_d() +
  scale_color_viridis_d() +
  theme_standard +
  theme(axis.text.x = element_blank())
input.file2 = "/Users/slecam/Desktop/COURS/BIODIV&CONS/4-TP structure tigres/169_tigers_structure.stru.txt"

struct2geno(file = input.file2, TESS = FALSE, diploid = TRUE, FORMAT = 2,
            extra.row = 1, extra.col = 2, output = "tiger.geno")

obj.snmf = snmf("tiger.geno", K=1:8, ploidy=2, entropy=T, repetition=10, alpha=100, project="new")
plot(obj.snmf, cex = 1.2, col = "lightblue", pch = 19)

Comment choisir le modèle avec la meilleure probabilité postérieur?
"We use SNMF’s cross-entropy criterion to infer the best estimate of K. The lower the cross-entropy, the better our model accounts for population structure. Sometimes cross-entropy continues to decline, so we might choose K where cross entropy first decreases the most."
Souvent la solution n'est pas de choisir une K mais bien comparer les informations amenées par plusierus K Ici la première diminution de K a lieu à K=4 mais ça diminue également à K=6 Comparons les 2.
Comment représenter les données?
qmatrix = Q(obj.snmf, K = 4, run=27)

barplot(t(qmatrix), col = c("orange","violet","lightgreen","gray","red","darkblue"), border = NA, space = 0,
        xlab = "Individuals", ylab = "Admixture coefficients")

par(mar=c(4,4,0.5,0.5))
pops<-levels(dataset@pop)
barplot(t(qmatrix), col=RColorBrewer::brewer.pal(9,"Paired"), 
        border=NA, space=0, xlab="Individuals", 
        ylab="Admixture coefficients")
#Add population labels to the axis:
for (i in 1:length(pops)){
  axis(1, at=median(which(dataset@pop==pops[i])), labels=pops[i])}
