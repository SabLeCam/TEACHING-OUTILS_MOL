# Phylogénies de la lactate déshydrogénase des daphnies (Neighbor-Joining et Maximum de vraisemblance)

[source](https://static1.squarespace.com/static/5459da8ae4b042d9849b7a7b/t/57ea64eae58c62718aa34769/1474979059782/Nesin_Winternitz_Practical_1and2.pdf)(Introduction to the phylogenetic (comparative) method, J. Wintternitz 2016)


Objectif : Construire des arbres phylogénétiques pour observer l’évolution de la lactate déshydrogénase chez les daphnies à l’aide de différentes méthodes de reconstruction.


Les microcrustacés du genre *_Daphnia_* sont très répandus dans les écosystèmes d’eau douce de l’hémisphère nord. On y retrouve deux espèces écologiquement distinctes : *_Daphnia pulex_* dans les étangs et *_Daphnia pulicaria_* dans les lacs. Ces deux espèces sœurs ont divergé il y a environ 80,000 ans. Elles produisent des hybrides F1 qui peuvent se rétro-croiser avec les espèces parentales. Ces deux espèces sont morphologiquement très semblables et ne peuvent être distinguées qu’à l’aide du **gène ND5** du génome mitochondrial ainsi que du gène de la **lactate désydrogénase A**.  *_D. pulex_* possède l’allèle *‘slow’* pour cet enzyme et *_D. pulicaria_* la forme *‘fast’*. Leurs hybrides F1 sont hétérozygotes pour les deux allèles.
Le séquençage complet du génome de *_Daphnia pulex_* a révélé une deuxième forme de la lactate désydrogénase : *la B*. Ce gène est vraisemblablement apparu par duplication et nous ne connaissons pas sa fonction. Il est possible qu’il s’agisse d’un pseudogène. Les mutations sont plus susceptibles de s’accumuler sur les pseudogènes puisqu’ils ne sont pas sous sélection pour maintenir une fonction précise dans le génome. Ils sont alors libres d’accumuler plus de mutations. 

<p align="center">
<img width="149" alt="image" src="https://github.com/SabLeCam/OUTILS_MOL/assets/20643860/01a513af-ab2c-4323-bd0a-88097c74356f">
</p>


Ici vous trouverez comment réaliser l'enesemble des analyses de phylogénétique comparative sous R

## package R nécessaires pour ces analyses

```r
adegenet
ape
geiger
phytools
stringr
msa
pegas
phangorn
```
installer les packages necessaires et les appeler

```r
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("msa")

install.packages("adegenet") ##continue with following package names
library(adegenet) ##continue with following packages

adegenet
ape
geiger
phytools
stringr
msa
pegas
phangorn
```

## Importer des données



```r
seq_LDHA <-read.dna("LDHA_nuc_121.fas", format="fasta")
##check your sequences
image(seq_LDHA)
```

## Aligner les sequences

Ici les séquences sont déjà alignées mais pour information voici la marche à suivre pour les aligner

```r
seq_LDHA <-Biostrings::readDNAStringSet("LDHA_nuc.fas")# load file
seq_LDHA
Aln_LDHA <- msa(seq_LDHA, method = "ClustalOmega", order="aligned")
#export alignment as DNAbin object
seq_align_LDHA <- msaConvert(Aln_LDHA, type="seqinr::alignment")
```
## Calculer la diversité nucléotidique totale et par groupe

```r
# diversité nucleotidique totale
nuc.div(seq_LDHA)
```
```r
# diversité nucleotidique totale par groupe
LDHA_C <- seq_LDHA[grep("LDHA-C",rownames(seq_LDHA)),]
LDHA_C
nuc.div(LDHA_C)
```
>Faites de même avec tous les groupes (C, CX et X)


## Construire des phylogénies à partir des données moléculaires

## Partie I: Neighbor-joining



Calculer une matrice de distance 

La réduction des alignements de séquences dans une matrice de distances par paires permet une estimation rapide des arbres phylogénétiques (au prix des informations supplémentaires que ces séquences contiennent). Toutefois, pour produire une matrice de distance, vous devrez choisir le modèle d’évolution des nucléotides (ou des protéines) qui correspond le mieux à vos données, en effectuant un test du rapport de vraisemblance.

```r
seq_LDHA2<-as.phyDat(seq_LDHA)
mt <- modelTest(seq_LDHA2)
print(mt)
dna_dist <- dist.ml(seq_LDHA2, model="JC")
```

```r
?dist.dna()###regarder les différents modèles d'évolution moléculaire disponibles
distxj<-dist.dna(seq_LDHA, model="JC69")
```

**Neighbor-Joining**: méthode de reconstruction phylogénétique se basant sur les distances évolutives. À partir des distances mesurées entre chaque séquence on part d’un arbre en étoile (i.e les séquences sont placées arbitrairement les unes par rapport aux autres). Puis, on calcule pour chaque paire de séquence la longueur des branches entre celles-ci. On retient la paire de séquence pour laquelle la longueur minimale des branches est obtenue; puis on les regroupe. On recommence avec les N-2 séquences restantes et jusqu’à avoir placé toutes les séquences.

# Construire un arbre avec la méthode du neighbor-joining
```r
tree <- nj(distxj_LDHA)
#plot a basic tree
plot.phylo(tree, type="phylogram", cex=0.4)
```
<img width="932" alt="image" src="https://github.com/SabLeCam/OUTILS_MOL/assets/20643860/6b257e5e-c48f-4661-a6fb-11218006b2a1">

Exemples de représentations de l'arbre
```r
par(mfrow=c(2,2)) #plot 2 rows by 2 columns
plot.phylo(tree, type="phylogram", cex=0.3)
plot.phylo(x=tree, type="cladogram", edge.width=2, cex=0.3)
plot.phylo(x=tree, type="fan", edge.width=2, edge.lty=2, cex=0.3)
plot.phylo(x=tree, type="radial", edge.color="red", edge.width=2, edge.lty=3, cex=0.3)
par(mfrow=c(1,1))
```
<img width="1330" alt="image" src="https://github.com/SabLeCam/OUTILS_MOL/assets/20643860/ddea800e-3259-4176-a5f3-c429a11539e7">


Bootstrap pour evaluer la robustesse des noeuds: technique de ré-échantillonnage qui permet d’estimer la probabilité de l’existence de chaque branche interne.


```r
tree$tip.labeln# connaitre la position de notre outgroup
outgroup<-1  

#fonction pour enraciner l'arbre:
foo <- function(xx) root(nj(dist.dna(xx)), outgroup)
tr <- foo(seq_LDHA) 
bp <- boot.phylo(tr, seq_LDHA, foo, B=1000) #1000 bootstraps


#représentation graphique de l'arbre avec les valeur de boostrap pour les noeuds
plot(tr, cex=0.3)
nodelabels(round(bp/10), cex=0.5, adj=c(1,-0.2),frame="none")

```
<img width="1220" alt="image" src="https://github.com/SabLeCam/OUTILS_MOL/assets/20643860/c3af8203-ba90-490a-873f-e0ef0b7337a0">

mais ça n'est pas très lisible...


Représenter les groupes plutôt que les noms de séquence et utiliser des symbolesou des couleurs

```r
###creer un vecteur de nom  de groupe

gpname<-substr(tree$tip.label, 0, 7)
table(gpname)

#présenter les groupes par formes

plot.phylo(x=tr, type="phylogram", show.tip=FALSE, lwd=3, main="Neighbour-Joining tree LDHA")
#add axis with distances
axisPhylo()
#use symbols (pch) for tip labels
tiplabels(frame="none", pch=rep(x=0:3, times=c(2,24,36,59)), lwd=2, cex=0.8)
legend(x="bottomleft", legend=c("Outgroup","C","CX","X"), border="black", pch=0:7, pt.lwd=2, pt.cex=1.5, bty="o", bg="lightgrey", box.lwd=1, cex=1.2, title="Famille")
```
<img width="1376" alt="image" src="https://github.com/SabLeCam/OUTILS_MOL/assets/20643860/d68b5c14-0f5e-4040-a52c-8888a3a9f22d">

ou

```r
#présenter les groupes par couleur
plot.phylo(x=tr, type="phylogram", show.tip=FALSE, lwd=3, main="Neighbour-Joining tree LDHA")
#add axis with distances
axisPhylo()

cols<-setNames(c("Green","White","Blue","Red"), sort(unique(gpname)))
tiplabels(pie=to.matrix(gpname,sort(unique(gpname))),piecol=cols,cex=0.5)
nodelabels(round(bp/10), cex=0.5, adj=c(1,-0.2),frame="none")
legend(x="bottomleft", legend=c("C","CX","Outgroup","X"), border="black",
       fill=cols, pt.lwd=2, pt.cex=1.5, bty="o", bg="lightgrey", box.lwd=1, cex=1.2, title="Famille")

 ```
<img width="1339" alt="image" src="https://github.com/SabLeCam/OUTILS_MOL/assets/20643860/2657da71-3469-4346-939b-495655baf0d9">


>Refaire  le même exercice pour LDHB-nuc.fas ainsi que pour les deux fichiers de protéine.Utilisez le script générique à la fin du document


## Partie II: Maximum de vraisemblance

Maximum de vraisemblance: Méthode qui mathématise le processus évolutif. Cette méthode consiste à définir les probabilités de tous les évènements évolutifs possibles à partir de n’importe quelle séquence ancestrale jusqu’aux feuilles et à trouver  le scenario le plus probable qui a permis d’aboutir aux séquences analysées. Le scénario se matérialise sous la forme d’un arbre avec une certaine topologie et des longueurs de branches particulières. 

```r
#on recalcule l'arbre de NJ avec un autre package pour avoir le bon format
seq_LDHA_NJ<-NJ(seq_LDHA)

fit <- pml(seq_LDHA_NJ, seq_LDHA2)
print(fit)
fitJC <- optim.pml(fit, model = "JC", rearrangement = "stochastic")
logLik(fitJC)
bs <- bootstrap.pml(fitJC, bs=100, optNni=TRUE, multicore=TRUE, control = pml.control(trace=0))
```

```r
BStree<-plotBS(midpoint(fitJC$tree), bs, p = 50, type="p")
bst <- BStree$node.label


rootedtree <- root(as.phylo(fitJC$tree), outgroup)
#présenter les groupes par couleur
plot.phylo(x=rootedtree, type="phylogram", show.tip=FALSE, lwd=3, main="Max Vraisemblance tree LDHA")
#add axis with distances
axisPhylo()

cols<-setNames(c("Green","White","Blue","Red"), sort(unique(gpname)))
tiplabels(pie=to.matrix(gpname,sort(unique(gpname))),piecol=cols,cex=0.2)
nodelabels(bst, cex=0.8, adj=c(1,-0.6),frame="none")
```


<img width="1354" alt="image" src="https://github.com/SabLeCam/OUTILS_MOL/assets/20643860/bef4b665-63b3-4912-a865-d611e48455f6">
>Refaire pour le même exercice pour LDHB-nuc.fas ainsi que pour les deux fichiers de protéine


## Phylogénéitque comparative

Utilisez le "scritp générique dna" pour obtenir les arbres phylogénétiques des gènes ND5 et LDHB

Utilisez le "scritp générique aa" pour obtenir les arbres phylogénétiques à partir des séquences de protéines de LDHA et LDHB



>Décrivez les arbres produits par les deux méthodes. Est-ce-que les espèces sont bien séparées pour le gène mitochondrial ND5, le gène LDHA et LDHB ?

>Ou se situent les hybrides ?

>Rapporter les indices de diversité nucléotidiques (pi) pour chacune des espèces pour la LDHA et la LDHB ainsi que pour l’ensemble des séquences pour la LDHA et la LDHB.

>1.Quelle est la divergence nucléotidique entre les espèces pour les séquences nucléotidiques de LDHA et LDHB ?

>2. Quelle espèce possède le plus grand niveau de polymorphisme nucléotidique pour la LDHA et pour la LDHB ?

>3. Est-ce-que le gène LDHA est plus variable que le gène LDHB au niveau nucléotidique? Au niveau des acides aminés ?

>4. Est-ce-que les arbres pour la LDH sont congruents (concordent) avec la phylogénie du gène mitochondrial ND5 du complexe Daphnia pulex ?  Si non pourquoi ?


## Script générique pour les sequence nucléotidique pour tous les gènes:
```r
########generic script dna#####


dna<-"ND5_TP.fasta"
#import data
seq <-read.dna(dna, format="fasta")
gene_name<-substr(rownames(seq),0, 4)

#Model Testing and Distance Matrices
seq2<-as.phyDat(seq)
mt <- modelTest(seq2)
print(mt)



distxj<-dist.dna(seq, model="JC69")
tree <- nj(distxj)

tree$tip.label
outgroup<-1 #"LDHB-EPX_CZR" 

#fonction pour enraciner l'arbre:
foo <- function(xx) root(nj(dist.dna(xx)), outgroup)
tr <- foo(seq) 
bp <- boot.phylo(tr, seq, foo, B=1000) #1000 bootstraps


###creer un vecteur de nom  de groupe
gpname<-substr(tree$tip.label, 0, 7)
table(gpname)

#présenter les groupes par couleur
plot.phylo(x=tr, type="phylogram", show.tip=FALSE, lwd=3, main= paste("NJ tree ", gene_name[1]))
#add axis with distances
axisPhylo()

cols<-setNames(c("Green","White","Blue","Red"), sort(unique(gpname)))
tiplabels(pie=to.matrix(gpname,sort(unique(gpname))),piecol=cols,cex=0.2)
nodelabels(round(bp/10), cex=0.5, adj=c(1,-0.2),frame="none")
legend(x="bottomleft", legend=c("C","CX","Outgroup","X"), border="black",
       fill=cols, pt.lwd=2, pt.cex=0.8, bty="o", bg="lightgrey", box.lwd=1, cex=0.8, title="Famille")

#maximum de vraisemblance
seq_NJ<-NJ(distxj)
fit <- pml(seq_NJ, seq2)
print(fit)
fitJC <- optim.pml(fit, model = "JC", rearrangement = "stochastic")
#pour ND5
#fitJC <- optim.pml(fit)
logLik(fitJC)
bs <- bootstrap.pml(fitJC, bs=100, optNni=TRUE, multicore=TRUE, control = pml.control(trace=0))

BStree<-plotBS(midpoint(fitJC$tree), bs, p = 50, type="p")
y <- BStree$node.label

gpname<-substr(fitJC$tree$tip.label, 0, 7)
table(gpname)
rootedtree <- root(as.phylo(fitJC$tree), outgroup)
#présenter les groupes par couleur
plot.phylo(x=rootedtree, type="phylogram", show.tip=FALSE, lwd=3, main=paste("Max vraisemblance", gene_name[1]))
#add axis with distances
axisPhylo()

cols<-setNames(c("Green","White","Blue","Red"), sort(unique(gpname)))
tiplabels(pie=to.matrix(gpname,sort(unique(gpname))),piecol=cols,cex=0.2)
nodelabels(y, cex=0.8, adj=c(1,-0.6),frame="none")
legend(x="bottomleft", legend=c("C","CX","Outgroup","X"), border="black",
       fill=cols, pt.lwd=2, pt.cex=0.8, bty="o", bg="lightgrey", box.lwd=1, cex=0.8, title="Famille")
```

































