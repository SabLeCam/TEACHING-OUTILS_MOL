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
Récupérer les données (RAPPEL)

```r
#setwd("C:/Users/marie/Documents/Deuxieme_cycleUQAR/")
lobster = read.csv("Lobster_SNP_Genotypes.csv")
```
 ```r
 dataset.hfstat <- genind2hierfstat(dataset)
basicstat <- basic.stats(dataset, diploid = TRUE, digits = 2) 
names(basicstat)
```

On compare maintenant les patrons de diversité génétique entre les populations par rapport à la diversité globale.

## Matrice de Fst par paire de population (estimateur du Fst de Weir de Cockerham(1984)

```r
mat.obs <- pairwise.WCfst(dataset.hfstat)
```
Représenter cette matice avec une heatmap
```r
melted_matobs <- melt(mat.obs, na.rm = TRUE)
colnames(melted_matobs)<-c("pop1","pop2","value")

ggheatmap <- ggplot2::ggplot(melted_matobs, aes(pop1, pop2, fill = value)) +
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red",  
                       midpoint = 0.15, limit = c(0,0.3), space = "Lab" ) +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))
ggheatmap
```

## Analyses PCA (principle components analysis) sur le jeux de homard

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
