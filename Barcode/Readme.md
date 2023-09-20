# Exercice Code-barre

Objectif: Identifier l’espèce à laquelle appartient un individu en utilisant un fragment du gène mitochondrial du cytochrome oxydase I; en comparant la séquence de l’individu à une base de données. 

## A)  Identification d’une espèce par codebarre

<img width="1440" alt="image" src="https://user-images.githubusercontent.com/20643860/216438975-07312084-eda5-4fa7-b6b9-d33f7712458b.png">


- allez à l’adresse suivante : http://www.barcodinglife.org/
- puis dans l’onglet **Identification**
- copiez votre séquence et collez la dans le carré en dessous de **Enter** sequences in fasta format.
ATTENTION il faut que votre séquence soit d’au moins 500bp.
- cliquez sur **Submit**

Le résultat apparait en plusieurs sections:
  - <ins>**Search Result**</ins>: Dans cette section, vous trouverez le nom de l’espèce à laquelle votre séquence correspond, ainsi que des informations relatives à cette espèce si elles existent dans ce site. 
  - <ins>**Identification Summary**</ins>: Ici vous trouverez la probabilité d’identification de votre séquence en fonction du niveau taxonomique. 
  - <ins>**Similarity Scores of Top 99 Matches**</ins>: La similarité des 99 séquences les plus proches de votre séquence.
  - <ins>**TOP 20 Matches**</ins>: Les 20 entrées les plus proches de votre séquence. Vous trouverez leur identification complète ainsi que le pourcentage de similarité avec votre séquence et leur statut (publique ou privé). Toutes les séquences ne sont pas accessibles au grand public, car certaines font parties d’études qui ne sont pas encore publiées. Ces séquences privées existent mais vous ne pouvez accéder à d’autres informations.
  - <ins>**Sampling Sites For Top Hits (>98% Match)**</ins>: Ici les points sur la carte représentent les sites où les séquences ayant plus de 98% de similarité avec la vôtre ont été échantillonnées.

En absence de référence :

BOLD est une base de données mais qui est *incomplète* et il se peut que la référence n’existe pas. Dans ce cas, un message *no hit match* va s’afficher.
BOLD vous suggère d’aller chercher dans une autre base de données BLAST.

<img width="1189" alt="image" src="https://user-images.githubusercontent.com/20643860/216441493-58834143-55c1-4e59-9e8f-822a4d7e4307.png">


  - Cliquez sur **Blast Sequence On Genbank**. Votre séquence sera automatiquement rentrée sur le site de BLAST
  - Laissez les paramètres par défaut. Ici la recherche va s’effectuer dans la base de données nucléotidique. 
  - Cliquez sur **Highly similar sequences (megablast)**. Cette option va vous permettre de rechercher les séquences les plus similaires.
  - Cliquez sur **BLAST**

Après un certain temps (de l’ordre de quelques minutes maximum), votre résultat va s’afficher en plusieurs parties :
  - <ins>**Nucleotide Sequence**</ins>: Dans cette section, vous trouverez un résumé de votre requête : longueur de votre séquence, à quelle base de données vous l’avez comparé, etc...
  - <ins>**Graphic Summary**</ins>: Ici les résultats les plus proches sont résumés sous forme de traits et de couleurs. La longueur des traits représente la longueur de la séquence trouvée dans BLAST. Les couleurs allant du noir au rouge symbolisent le score d’alignement entre la séquence trouvée et votre séquence. Quand c’est noir c’est faible (<40), quand c’est rouge c’est très similaire (>200). Si vous cliquez sur une ligne vous allez été redirigé à la description détaillée de la séquence.
  - <ins>**Descriptions**</ins>: La description détaillée des résultats précédents. Vous trouverez pour chaque séquence son score d’alignement, sa couverture par rapport à votre séquence, sa e-value, le pourcentage d’identité et son numéro d’accession.
score d’alignement : quand deux séquences sont comparées, une valeur d’alignement (score) va être générée. Ce score va dépendre du nombre de mismatchs, de gaps et de matchs.
query-cover= couverture: vous pouvez comparée des séquences de différentes longueur, par exemple une plus longue et une plus courte. Votre séquence plus courte ne couvrira qu’un certain pourcentage de votre séquence plus longue.
  - <ins>**e-value**</ins>: la probabilité d’observer ce résultat par hasard. Plus la e-value est petite et moins vous avez de chance d’obtenir ce résultat par hasard.
  - <ins>**Identity**</ins> : pourcentage d’identité. Nombre de bases identiques entre votre séquence et celle comparée.
  - <ins>**Accession**</ins> : Identifiant de Genbank.

```
Pour chacune des séquences inconnues, vérifiez s’il y a un match dans BOLD.
Si vous devez aller dans genbank, quelle est l’espèce la plus proche ? 
Donnez les valeurs de couverture, la e-value, l’identité ainsi que le numéro d’accession dans genbank
```


## B) Détermination du barcode gap chez Aurelia, littorina, Zostera,Salvelinus (crédit photo wikipédia)

<p float="left">
<img width="200" alt="image" src="https://user-images.githubusercontent.com/20643860/216440139-bd7a70a9-aadc-412e-b431-9d9a56e1deff.png">
<img width="200" alt="image" src="https://user-images.githubusercontent.com/20643860/216439922-9999cbf2-6b6a-472a-aa2b-7696468522a1.png">
<img width="200" alt="image" src="https://user-images.githubusercontent.com/20643860/216440032-1f48bbd9-b05b-4f17-b176-341151bb1f0b.png">
<img width="200" alt="image" src="https://user-images.githubusercontent.com/20643860/216440886-cbe1e785-8042-4154-9224-de0e8b6e712f.png">
</p>




Sur le site de BOLD SYSTEMS (https://www.boldsystems.org), 
  - allez dans Taxonomy, et rentrer le nom du genre recherché
  - cliquez sur **public data** puis **combined** <img width="256" alt="image" src="https://user-images.githubusercontent.com/20643860/216434548-314a77a7-8d6e-461a-b055-2414f2979fbf.png">
  - allez rechercher le fichier txt téléchargé, copier le et coller le dans un fichier excel.
  - copier les **process id** dans la 1ère colonne
  
  - retourner sur Bold, cliquez sur **workbench** (username: biodivcons, password:FAU72816-05)
  - puis **record search** 
  - coller les données **process id** (cochez *include public records*) et **search**

<img width="1080" alt="image" src="https://user-images.githubusercontent.com/20643860/217006962-b457f55a-3d96-4fb9-8e74-9ae3e0305363.png">


  - selectionnez **all data**
  - allez dans l'onglet **_Sequence analyses_** et **taxon ID tree**
  - puis cochez **Muscle** dans alignement et dans **Colorize tree based on Barcode** : Bin 
  (Pour Zostera, faite un arbre avec rbcl et ensuite avec matk)
  - puis cliquez sur ‘Build tree’ et visualiser le ficher pdf
  - inspectez l’arbre coloré selon les Bin (individus regroupés par cluster de 2.2 %)
```
2)	Est-ce-que l’arbre parvient à résoudre les différentes espèces ?
  - Puis cliquez sur ‘Build tree’ et visualiser le ficher pdf
Est-ce-que globalement les individus portant les mêmes noms d’espèces se retrouvent dans des groupes (clusters) séparés ?  Indiquez si en général c’est bien regroupé ou s’il y a des problèmes d’identification pour certains genres..
3)	Indiquez les espèces qui causent problème 
```

Allez dans barcode gap analysis

```
4)	Quelle est la distance moyenne intraspécifique ? Quelles sont les valeurs minimum et maximum ?
La valeur de la distance moyenne est indiquée sous l’histogramme Mean intraspecific distance. Vous trouverez les valeurs min et max de la divergence intraspécifique sur l’axe des X du premier graphique à gauche (max intraspécifique divergence)
5)	Quelle est la distance à l’espèce voisine la plus proche ?
La moyenne se trouve sous l’histogramme ‘distance to nearest neighbor
6)	Est-ce-qu’il y a un barcode gap ?
Pour répondre à cette question, regarder s’il y a beaucoup de valeurs sous la diagonale rouge du graphique ‘max intraspecific distance to nearest neighbor’ et aussi dans l’analyse ABGD, si vous voyez un ‘gap’ entre la distribution des fréquences.

7)	Y-a-t’il des problèmes de classification ?
Vous pouvez voir ces problèmes quand vous cliquez sous ‘show warnings only’ dans barcode analyses dans bold.
``

