# TP Analyses de différenciation entre populations

( France Dufresne et Sabrina Le Cam, d'après un document produit par Tom Jenkins source https://tomjenkins.netlify.app/tutorials/r-popgen-getting-started/)

Nous repartons du jeux de données de l'étude de [Jenkins et al. 2019](https://onlinelibrary.wiley.com/doi/10.1111/eva.12849) sur le Homard bleu avec lequel nous avons fait le [TP Hardy Weinberg](https://github.com/SabLeCam/OUTILS_MOL/tree/main/HardyWeinberg)
Récupérer les données
```r

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

Inférence bayésienne de la structure de population
#input les données, fichier avec hearder différent

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
