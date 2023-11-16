# TP8 Structure génétique chez les zostères à l’aide de données génomiques

Objectif du laboratoire : Se familiariser avec la manipulation de données génomiques afin de réaliser des analyses statistiques dans R.
La PCoA :
## Qu’est-ce qu’une PCOA ?
	L’Analyse des Coordonnées Principales (PCoA) est une analyse de type multivariée, ça signifie qu’elle intègre une variable dépendante (VD) et plus d’une variable indépendante (VI). Une PCoA est une analyse semblable à une PCA et permet de simplifier des données complexes. L’objectif de cette analyse est de réaliser un graphique qui nous permettra d’interpréter les différences génétiques entre des groupes et des individus. La PCoA va créer autant de nouvelles variables (PC) qu’il y a de variables indépendantes initiales (individus) et les deux plus informatives seront projetées. Chaque nouvelle variable contiendra de l’information sur l’ensemble des variables indépendantes initiales (individus).
	Les données qui devront être importées dans R seront sous la forme d’une variable dépendante (les loci) et de plusieurs variables indépendantes (les individus).
 
	1 allèle référence et 0 allèle alternative :0
	0 allèle de référence et 1 allèle alternative : 2
	1 allèle référence et 1 allèle alternative : 1


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
