# Structure génétique chez les zostères à l’aide de données génomiques

*_Objectif du laboratoire : Se familiariser avec la manipulation de données génomiques afin de réaliser des analyses statistiques dans R._*


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




## Manipulation de données

Nous allons voir comment , au préalable de l'analyse de la diversité génétique en elle même, ces données méritent d'être explorées et filtrées afin de s'assurer de leur robustesse (taux de données manquantes, neutralité...)

Voici les packages R que nous allons utiliser dans ce TP. Les installer (si nécessaire) puis les appeler

```r
BiocManager::install(c("SNPRelate", "qvalue"))
install.packages("plotly")
```
```r
library(dartR)
library(BiocManager)
library(devtools)
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

Examiner les données (nb d'individus, nb de locus, % de données manquantes)

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














 
## Qu’est-ce qu’une PCOA ?
	L’Analyse des Coordonnées Principales (PCoA) est une analyse de type multivariée, ça signifie qu’elle intègre une variable dépendante (VD) et plus d’une variable indépendante (VI). Une PCoA est une analyse semblable à une PCA et permet de simplifier des données complexes. L’objectif de cette analyse est de réaliser un graphique qui nous permettra d’interpréter les différences génétiques entre des groupes et des individus. La PCoA va créer autant de nouvelles variables (PC) qu’il y a de variables indépendantes initiales (individus) et les deux plus informatives seront projetées. Chaque nouvelle variable contiendra de l’information sur l’ensemble des variables indépendantes initiales (individus).


## Faire une PCoA sur Rstudio :

Dans un premier temps, il faut télécharger et activer les librairies que nous allons utiliser dans l’analyse.
 
Ensuite, il faut donner à Rstudio l’accès à notre dossier dans lequel se trouvent nos fichier Dart. Pour ce faire nous utilisons la fonction « setwd ». Il faut prendre le chemin qui mène à votre dossier et le coller dans cette fonction. Attention au sens des « / » : il se peut que par default ils soient collés dans ce sens : « \ », dans ce cas il faut le changer.
Voilà un exemple (à modifier selon la destination sur votre ordinateur) :
 
	Il faut maintenant créer l’objet qui correspondra à notre fichier. Pour ce faire, il faut s’assurer de bien avoir deux fichiers enregistrés en .csv : 1) un qui correspond au fichier snp reçut de la compagnie Dart, 2) un qui correspond à identification des individus et des populations qu’il faut créer soit même (ici nous avons déjà fait ce fichier pour vous).
Forme fichier 1 :
 
Forme fichier 2 :
 
La fonction « gl.read.dart » est utilisée dans ce cas :
 
Dans cette fonction,
« filename = » correspond au nom de fichier qu’il faut mettre le nom du fichier suivi de .csv ; « ind.metafile = » correspond au nom du fichier correspondant à l’identification des population .csv ; « nas = » ici on inscrit ce qui correspond aux données manquante ; « topskip = » il s’agit du nombre de ligne à passer dans le fichier pour arriver à la ligne où se trouve les identifiant des individus ici 6 lignes.
	Une fois la fonction effectuée il faut vérifier que le fichier a bien été lu correctement. Pour cela, il faut cliquer sur l’objet dans l’environnement de Rstudio et vérifier le bon nombre d’individus (ici 161), le bon nombre de loci et les libellés des populations. 
 
1)	Utilisez la fonction nLoc(gl) pour déterminer le nombre de loci

2)	Utilisez la fonction nInd(gl) pour déterminer le nombre d’individus

3)	Utilisez la fonction nPop(gl) pour déterminer le nombre de populations

4)	Utilisez la fonction  popNames(gl) pour obtenir le nom des populations

PCoA sans filtres
	Nous allons d’abord faire une analyse PCoA sans faire de filtre au préalable. Pour ce faire il faut performer les fonctions suivantes.
 
Dans ce cas ci (sans filtre) il est rare que la fonction fonctionne ou donne quelque chose de correct, car parmi les 6365 loci il y en a forcement quelques un qui on des problèmes (monomorphes, trop de données manquantes…). Il faut donc faire des filtres pour éliminer ces loci.
PCoA avec filtres :

Quand on fait des tests génétiques, surtout avec les snp il faut faire des filtres pour garder que les bons loci et les bons individus. Il est très important d’éliminer des loci avec beaucoup de données manquantes, ceux qui ne sont pas informatifs (monomorphes, dupliqués…). Aussi, les individus avec trop de données manquantes peuvent altérer les résultats (en étant rapprochés des individus ayant eux aussi des données manquantes au lieu de ceux lui ressemblant génétiquement), pour cette raison il faut aussi les enlever. Il est aussi important de les exclure car ils peuvent altérer les résultats d’autres analyses qui peuvent être faites comme des analyses de diversité génétique par exemple.
Le premier filtre que nous allons faire est le filtre pour enlever les loci monomorphe grâce à cette fonction :
 
	Pour chaque filtre effectué nous ne supprimons pas le fichier précédent, nous créons à chaque fois un nouvel objet (ici : gl1). Nous pouvons maintenant regarder combien de loci ont été enlevé en cliquant sur « gl1 » dans l’environnement de RStudio.
 
	Il reste 6293 loci sur les 6365 du départ, il y avait donc 72 loci monomorphes dans notre jeu de donnée. Nous pouvons nous réessayer à faire la PCoA pour voir si ce filtre à changer quelque chose. Pour ce faire, nous reprenons la formule de la PCoA et nous changeons « gl » par « gl1 ». On change aussi le nom de l’objet « pc » pour qu’il n’y ait pas de suppression d’information. Cela donne donc :
(Attention cette fonction est assez longue notamment pour l’affichage du graphique, il fait bien attendre.)
 
	Avec déjà ce premier filtre nous obtenons un graphique de groupement. Cela signifie donc bien que le problème précédent était donc bien lié aux loci monomorphes.
 
Comment interprétez-vous le graphique obtenu ?
	Nous allons ensuite filtrer les loci et les individus avec trop de données manquantes pour les loci on vise 0.05% maximum de données manquantes et pour les individus 3% maximum. Nous voulons aussi enlever les doublons. Nous faisons ces filtres en plusieurs pour maximiser le nombre de loci et d’individus tout en enlevant un maximum de données manquantes.

Après tous les filtres, nous obtenons donc un total de 152 individus et 1136 loci.  
Nous pouvons faire la PCoA finale en utilisant la même fonction et en changeant « gl » par « gl7 » et « pc » par « pc7 ». Cela donne donc :
 
Le résultat obtenu est donc ce graphique :
 


Comment interprétez-vous le graphique obtenu ?
Il y a-t-il des différences avec le graphique précédent ? si oui expliquez
Reprenez les mêmes commandes pour les laminaires (filename= Report_DSacc21-6007_3_moreOrders_SNP_2. Indmetafile= laminaria2.csv)

En quoi les zostères diffèrent-elles des laminaires au niveau de leur structure populationnelle ?

Références pour approfondir le sujet :
Gruber, B., Unmack, P., Berry, O., & Georges, A. (2019). Introduction to dartR. User Manual.
