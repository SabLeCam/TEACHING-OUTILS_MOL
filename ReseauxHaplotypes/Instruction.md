# TP Les réseaux d’haplotypes

(France Dufresne)

<ins>Lecture</ins> : 
<img width="387" alt="image" src="https://github.com/SabLeCam/OUTILS_MOL/assets/20643860/d7f80ce5-89c9-4a9f-99bf-10dce47fe80c">

http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0001913 

<ins>Objectif</ins>: 
Faire un réseau d’haplotype. Un réseau d’haplotypes permet d’explorer les liens qui existent entre les différents haplotypes présents dans une même espèce. 
Il existe différentes méthodes de construction des réseaux.

<ins>Installation</ins> :
Aller à l’adresse suivante : http://popart.otago.ac.nz/downloads.shtml 
Télécharger *__popart__* selon votre système d’exploitation et installer le dans votre répertoire.

<ins>Glossaire</ins> :
Un haplotype= combinaison donnée de nucléotides/allèles polymorphes.
Diversité haplotypique= probabilité que deux allèles tirés au hasard soient différents.

## Part I: Faire les haplotypes  (déjà fait pour vous)

Pour faire un réseau d’haplotypes, il faut d’abord produire un fichier qui contient l’ensemble de nos séquences dans un format particulier (format fasta voir exemple). Ensuite nous devons aligner les séquences, c’est-à-dire les mettre en ordre selon les différentes positions sur l’ADN. Plusieurs logiciels sont utilisés pour les alignements (*__bioedit__* et *__mega__* qui sont déjà installés sur les ordinateurs). Une fois que les séquences sont alignées, il faut produire un fichier (format *_nexus_*) qui servira à produire les réseaux d’haplotypes. Ces fichiers sont produits par le logiciel *__DNAsp__* (http://www.ub.edu/dnasp/) qui regroupe les séquences identiques ensemble afin d’établir la fréquence de chaque haplotype. Ce fichier doit également inclure des coordonnées GPS pour chaque séquence afin de pouvoir placer la localisation de chaque groupe de séquence sur une carte. 



## Partie II : Réseaux.

- Aller à : http://www.popart.otago.ac.nz/index.shtml puis Download
- Choisir votre système d’exploitation puis télécharger le fichier.
- Installer le fichier dans votre répertoire. Cliquer sur **popart.exe** pour lancer l’installation (sous windows).
- Lorsque vous ouvrez le fichier, il y aura un message d'erreurs à propos de caractères inconnus, cliquez ok.
- Aller dans Network puis choisissez une méthode de construction (TCS network).
- Ils demandent d'inscrire un chiffre epsilon, entrez 33
- Le réseau apparait dans la fenêtre à droite. 


## Partie III : Réseaux en fonction du site d’échantillonnage.

**Pour** *__Gammarus oceanicus__* 
<img width="387" alt="image" src="https://github.com/SabLeCam/OUTILS_MOL/assets/20643860/02ccbf69-e8a4-4f3d-8be8-42e23ede57e5">


- Ouvrir le fichier  **Go_haplotype.nex** dans *__PopArt__*.
- Aller dans **Network** puis choisissez le type de construction.
- Le réseau apparait dans la fenêtre à droite. 
- Vous pouvez modifier les couleurs ainsi que la police en cliquant sur les points de couleurs ou en allant dans **Edit**.
- Exporter le graphique en allant dans **File** puis **Export graphics**. Vous pouvez choisir le format svg, png ou pdf. Par commodité préférez l’option pdf.

<p align="center">
  <img width="587" alt="image" src="https://github.com/SabLeCam/OUTILS_MOL/assets/20643860/97b140bb-abcf-4574-b7a0-559fd2ab9844">
  </p>

- Vous pouvez aussi visualiser la répartition de vos haplotypes en fonction du site d’échantillonnage en allant dans **View** puis switch to **map view**. Vous pouvez ajouter des légendes (échelle, rose des vents à votre carte en faisant clic droit puis **Info boxes**.
- Vous pouvez exporter de la même manière que le réseau.

 <p align="center">
  <img width="587" alt="image" src="https://github.com/SabLeCam/OUTILS_MOL/assets/20643860/da2aad71-908e-4fe7-8ea3-fe5bdd61e091">
  </p> 

> Note : Comme il y a plusieurs sites provenant de la même région, vous pourrez recolorer les haplotypes de ces sites avec la même couleur.

**Liste pour les couleurs**:
1. Labrador, Northern Peninsula, Bonavista, Fogo Island, Avalon Peninsula, Bonne Bay, Saint_John
2. Westmorland, Albert, St. John, Restigouche, Westmorland, Gloucester, Cumberland, Richmond, Pictou, Guysborough, Antigonish, Inverness, Halifax, queens, Inverness, twillingate, Northumberland, Richmond
3. Cote-Nord, Gaspesie, Charlevoix, Bas-Saint-Laurent, Madeleine
4. Hudson, Churchill, Kangiqsualujjuag
9. Norvège, Pologne, Suède, Baltique, Écosse, Mer Nord, Mer_Blanche,
10. Faroe, Groenland, Svalbard, Islande

## Partie IV : Analyses statistiques
Allez dans *__Statistics__* et cochez *__All stats__*

Rapporter les résultats sous forme de tableau.

Allez dans Statistics et cochez *__Amova__*
- Vous devez définir des groupes (Amérique du Nord vs Europe pour *_Gammarus oceanicus_*). 
- Rapporter les valeurs de `Fst`, `Fsc` et `Fct` pour les deux jeux de données.


> Question Gammarus:
> 1. Que nous indique ces réseaux d’haplotype sur l’impact des glaciations sur cette espèce des deux côtés de l’Atlantique ? Ou étaient les refuges glaciaires ?
> 2. Pouvez-vous formuler une hypothèse sur l’origine des haplotypes présents à Savalbard, au Groenland et à Churchill nord du Manitoba ?


### Poursuivez le TP pour les deux espèces arctiques: *_Themisto libellula_* (gauche) et *_Themisto abyssorum_* (droite)
<p float="left">
  <img width="387" alt="image" src="https://github.com/SabLeCam/OUTILS_MOL/assets/20643860/804fb46f-ac1d-44ec-8459-822633eee03e">
  <img width="387" alt="image" src="https://github.com/SabLeCam/OUTILS_MOL/assets/20643860/4acad8f2-89d4-4a82-8b66-0d9fa8c3bbdc">
  </p>


Pour *_Themisto libellula_*, ouvrez le fichier tlmt597.nex et pour *_Themisto abyssorum_*, ouvrez le fichier **ta498.nex**.
Allez dans *__Statistics_** et cochez *__All stats__*
Rapporter les résultats sous forme de tableau pour chaque espèce

Allez dans *__Statistics__* et cochez *__AMOVA__*
-Explorer les différents groupements qui maximisent la séparation des groupes pour chacune des espèces). 
-Rapporter les valeurs de  `Fst`, `Fsc` et `Fct` pour les deux jeux de données.
Cochez `D de Tajima` et rapporter les valeurs

Allez dans *__Statistics__* et cochez *__All stats__*
Rapporter les résultats sous forme de tableau.
Cochez `D de Tajima` et rapporter les valeurs


> Question Themisto:
> 1. Quelle espèce possède la plus grande diversité génétique ?
> 2. Est-ce-que les glaciations ont affecté les deux espèces de la même façon ? Quel scénario de recolonisation postglaciaire peut être envisagé pour ces espèces ?



