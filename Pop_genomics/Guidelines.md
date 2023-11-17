# Structure génétique chez les zostères à l’aide de données génomiques

*_Objectif du laboratoire : Se familiariser avec la manipulation de données génomiques afin de réaliser des analyses statistiques dans R._*

Références pour approfondir le sujet : Gruber, B., Unmack, P., Berry, O., & Georges, A. (2019). Introduction to dartR. User Manual.

**DArTSeq** est une méthode qui extrait la variation génomique reproductible à travers les génomes de nombreux individus à un coût abordable. La technique
digère l'ADN génomique à l'aide de paires d'enzymes de restriction (coupeurs). Quand l'ADN est coupé en deux endroits situés à une distance raisonnable l'un de l'autre, le fragment est disponible pour le séquençage à l’aide des plateformes de lecture courte Illumina. Par conséquent, les données sont représentatifs dans le sens où ils sont générés pour une sélection aléatoire de petits fragments de séquence uniquement, fragments qui présentent une variation au niveau de paires de bases uniques (SNP).
<p align="center">
<img width="680" alt="image" src="https://github.com/SabLeCam/OUTILS_MOL/assets/20643860/d1aef0a9-b383-489e-aae8-a3c7f7589793">
</p>
Les **SNP**, ou single nucleotide polymorphisms, sont des mutations d'une seule paire de bases à un locus nucléaire. Ce locus nucléaire est représenté dans l'ensemble de données par deux séquences qui, sur un locus hétérozygote, prennent deux états alléliques, l'un appelé l'état de référence, l'autre comme état alternatif ou SNP.
<p align="center">
<img width="703" alt="image" src="https://github.com/SabLeCam/OUTILS_MOL/assets/20643860/862d6576-7d0a-4a29-80e7-2a5143f880f8">
</p>
Les données peuvent être représentées dans un tableau de bases SNP (A, T, C ou G), avec deux états pour chaque individu à chaque locus dans un organisme diploïde.Alternativement, comme les données sont bialléliques, il est pratique de coder les données comme 0 pour les homozyogotes pour un allèle, 1 pour les hétérozygotes et 2 pourhomozygotes de l’autre allèle.L'allèle de référence est arbitrairement choisi comme l'allèle le plus courant, donc 0 est le
score pour un homozygote à l'allèle de référence, et 2 est le score pour un homozygote à l'allèle alternatif. NA indique que le SNP n’a pas pu être évalué.

 <img width="538" alt="image" src="https://github.com/SabLeCam/OUTILS_MOL/assets/20643860/cbd625a7-d3a1-4f25-ad7e-0c829691d16b"> <img width="533" alt="image" src="https://github.com/SabLeCam/OUTILS_MOL/assets/20643860/6cf8379f-9540-438f-9b7a-193ae0b792e6">


Nous allons analyser un jeu de données de zostères echantillonnées dans l'Estaire et le Golfe du Saint Laurent
<p align="center">
<img width="642" alt="image" src="https://github.com/SabLeCam/OUTILS_MOL/assets/20643860/c2f4b779-d39f-4a2b-8d87-55cc428976f9">
</p>

## Exploration des données

Nous allons voir comment , au préalable de l'analyse de la diversité génétique en elle même, ces données méritent d'être explorées et filtrées afin de s'assurer de leur robustesse (taux de données manquantes, neutralité...)

Voici les packages R que nous allons utiliser dans ce TP. Les installer (si nécessaire) puis les appeler

```r
install.packages("dartR")
library(dartR)
install.packages("devtools")
library(devtools)
install.packages("BiocManager")
BiocManager::install(c("SNPRelate", "qvalue"))
#ATTENTION dites non aux mises à jour !
install_github("green-striped-gecko/dartR")
library(dartR)
install.packages("plotly")
library(plotly)
library(adegenet)
library(ggplot2)
```

Ouverture de fichier avec la fonction *_gl.read.dart()_*.

Les données dartSeq sont composée de 2 fichier: un fichier de génotype et un fichier de métadata, c'est à dire toutes les infos complémentaire sur les individus génotypes (ex :lieu/date d'échantillonnage, taille, sexe, traitement expérimental ....)
<p align="center">
<img width="478" alt="image" src="https://github.com/SabLeCam/OUTILS_MOL/assets/20643860/0480474a-cdc5-4481-8c19-62a7e337d2b2">
</p>

```r
# Définir le répertoire de travail par défaut 
setwd("PATH_TO_TOUR_FILE")
gl <- gl.read.dart(filename = "Zost_bonfichier.csv", ind.metafile= "Popidzostbonfichier(1).csv", nas = "-", topskip = 6, probar = TRUE)
gl
```
<img width="574" alt="image" src="https://github.com/SabLeCam/OUTILS_MOL/assets/20643860/e206dd28-322d-4261-b126-50419c288340">

Une fois la fonction effectuée il faut vérifier que le fichier a bien été lu correctement. Pour cela, il faut cliquer sur l’objet dans l’environnement de Rstudio et vérifier le bon nombre d’individus (ici 161), le bon nombre de loci et les libellés des populations. 
 
1)	Utilisez la fonction nLoc(gl) pour déterminer le nombre de loci

2)	Utilisez la fonction nInd(gl) pour déterminer le nombre d’individus

3)	Utilisez la fonction nPop(gl) pour déterminer le nombre de populations

4)	Utilisez la fonction  popNames(gl) pour obtenir le nom des populations


Examiner la qualité des données (nb d'individus, nb de locus, % de données manquantes)

```r
#par locus
gl.report.callrate(gl)
#par individu
gl.report.callrate(gl, method='ind')
```
<img width="702" alt="image" src="https://github.com/SabLeCam/OUTILS_MOL/assets/20643860/555a0df5-4ee2-4bf0-8f7f-9062c5dc3f75">
<img width="711" alt="image" src="https://github.com/SabLeCam/OUTILS_MOL/assets/20643860/3c9d96e3-0d68-40c2-9b66-b597852038bd">

On peut aussi visualiser rapidement les données avec un *_smear plot_*
```r
gl.smearplot(gl)
```
![image](https://github.com/SabLeCam/OUTILS_MOL/assets/20643860/39454085-da08-42dc-8867-ab725c66ff8e)

## Manipulation de données

Presque systématiquement, ce type de données demande d'entreprendre diverses manipulations avant l'analyse en supprimant certains individus ou locus. On veut garder le maximum d'individus tout en réduisant au maximum le nombre de données manquantes qui peuvent induire du biais dans la analyses. Les loci avec beaucoup de données manquantes ne sont pas informatifs (monomorphes, dupliqués…) et les individus avec trop de données manquantes peuvent altérer les résultats comme des analyses de diversité génétique par exemple.
On commence par filtrer les SNPs qui n'ont que des données manquantes ou ceux qui sont monomorphes.
*_Cette étape est otpionnelle ici car elle a déjà été effectué en amont de ce TP (le fichier initial était beaucoup plus volumineux)_*

```r
###SNP filtering####
#snp <- gl.filter.allna(gl)
#snp2<-gl.recalc.metrics(snp)
#snp3<-gl.filter.monomorphs(snp2)
#snp4<-gl.recalc.metrics(snp3)
#snp4
```


Filtres pour enlever les données manquantes (loci 99.5% et individu 97%) et les doublons
```r
gl2 <- gl.filter.callrate(gl, method = "loc", threshold = 0.95) #filtre les loci avec plus de 5% de donnees manquantes

gl3 <- gl.filter.secondaries(gl2) #filtres les loci doublons

gl4 <- gl.filter.callrate(gl3, method = "ind", threshold = 0.80) #filtre les individus avec plus de 20% de donnees manquantes

gl5 <- gl.filter.callrate(gl4, method = "loc", threshold = 0.97)#filtre les loci avec plus de 3% de donnees manquantes

gl6 <- gl.filter.callrate(gl5, method = "ind", threshold = 0.97)#filtre les individus avec plus de 3% de donnees manquantes

gl7 <- gl.filter.callrate(gl6, method = "loc", threshold = 0.995)#filtre les loci avec plus de 0,5% de donnees manquantes

```
Revisualiser les données

```r
gl.report.callrate(gl7)

gl.report.callrate(gl7, method='ind')

gl.smearplot(gl7)
```
<img width="702" alt="image" src="https://github.com/SabLeCam/OUTILS_MOL/assets/20643860/707316cb-457d-4fbe-9a70-b493c56f71d8">


## Analyse de la structure de la diversité génétique avec une PCOA
### Qu’est-ce qu’une PCOA ?
L’Analyse des Coordonnées Principales (PCoA) est une analyse de type multivariée, ça signifie qu’elle intègre une variable dépendante (VD) et plus d’une variable indépendante (VI). Une PCoA est une analyse semblable à une PCA et permet de simplifier des données complexes. L’objectif de cette analyse est de réaliser un graphique qui nous permettra d’interpréter les différences génétiques entre des groupes et des individus. La PCoA va créer autant de nouvelles variables (PC) qu’il y a de variables indépendantes initiales (individus) et les deux plus informatives seront projetées. Chaque nouvelle variable contiendra de l’information sur l’ensemble des variables indépendantes initiales (individus).

```r
pc7 <- gl.pcoa(gl7, nfactors=5)
names(pc7)
barplot(pc7$eig/sum(pc7$eig)*100, )#graphique de pourcentage de distribution des axes
gl.pcoa.plot(pc7, gl7, label = "ind", xaxis=1, yaxis=2)
```
<p align="center">
<img width="702" alt="image" src="https://github.com/SabLeCam/OUTILS_MOL/assets/20643860/a6abf8d9-9b5a-40a9-844e-f655808381ec">
</p>

>Comment interprétez-vous le graphique obtenu ?


>Reprenez les mêmes commandes pour les laminaires (filename= Report_DSacc21-6007_3_moreOrders_SNP_2. Indmetafile= laminaria2.csv)
>Pour ces fichiers, il y a des populations qui proviennent des USA et d'autres régions du monde et que nous n'utiliserons pas dans l'analyse. Il
faut donc les enlever avec cette fonction

```r
gl1<- gl.drop.pop(gl, pop.list=c("Denmark","Norway","russia","Iceland","Ireland","Korea","Alaska","Maine","Connecticut","Rhode_island","CapeCod","NewYork","NewHampshire","washington"))
```

Voici la carte d'échantillonnage des laminaires.
<p align="center">
<img width="752" alt="image" src="https://github.com/SabLeCam/OUTILS_MOL/assets/20643860/32974c21-b99c-4b86-911f-87190221821b">
</p>

>En quoi les zostères diffèrent-elles des laminaires au niveau de leur structure populationnelle ?



