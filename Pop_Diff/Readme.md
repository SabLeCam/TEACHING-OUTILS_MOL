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
percent = pca1$eig/sum(pca1$eig)*100
barplot(percent, ylab = "Genetic variance explained by eigenvectors (%)", ylim = c(0,12),
        names.arg = round(percent, 1))
```

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
  geom_text(data=centroids,aes(label=pop), size=5,position="jitter")

gpca

# Create a data.frame containing individual coordinates
ind_coords = as.data.frame(pca1$li)

# Rename columns of dataframe
colnames(ind_coords) = c("Axis1","Axis2","Axis3")

# Add a column containing individuals
ind_coords$Ind = indNames(lobster_gen_sub)

# Add a column with the site IDs
ind_coords$Site = lobster_gen_sub$pop

# Calculate centroid (average) position for each population
centroid = aggregate(cbind(Axis1, Axis2, Axis3) ~ Site, data = ind_coords, FUN = mean)

# Add centroid coordinates to ind_coords dataframe
ind_coords = left_join(ind_coords, centroid, by = "Site", suffix = c("",".cen"))

# Define colour palette
cols = brewer.pal(nPop(lobster_gen_sub), "Set1")

# Custom x and y labels
xlab = paste("Axis 1 (", format(round(percent[1], 1), nsmall=1)," %)", sep="")
ylab = paste("Axis 2 (", format(round(percent[2], 1), nsmall=1)," %)", sep="")

# Custom theme for ggplot2
ggtheme = theme(axis.text.y = element_text(colour="black", size=12),
                axis.text.x = element_text(colour="black", size=12),
                axis.title = element_text(colour="black", size=12),
                panel.border = element_rect(colour="black", fill=NA, size=1),
                panel.background = element_blank(),
                plot.title = element_text(hjust=0.5, size=15) 
)

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
  ggtitle("Lobster PCA")+
  # custom theme
  ggtheme
```


## Inférence bayésienne de la structure de population
#input les données, fichier avec hearder différent
# FORMAT DATA ----
# extract allele counts
allelecount <- as.data.frame(tab(gen_n))
# replace missing values with 9
allelecount[is.na(allelecount)] <- 9
# write lfmm input file
write.lfmm(allelecount, "data/POPGEN/GTF.lfmm")
geno <- lfmm2geno(input.file = "data/POPGEN/GTF.lfmm",
                  output.file = "data/POPGEN/GTF.geno")
# RUN ANALYSIS ----
# calculate coancestry coefficents (K = 1-10 for 10 iterations)
snmf.obj1 <- snmf(input.file = "data/POPGEN/GTF.geno",
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
  labs(x = "number of ancestral populations", y = "cross-entropy criterion") +
  theme_standard
  
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
