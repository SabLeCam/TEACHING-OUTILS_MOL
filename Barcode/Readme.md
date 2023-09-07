
TP   Code-barre
Objectif: Identifier l’espèce à laquelle appartient un individu en utilisant un fragment du gène mitochondrial du cytochrome oxydase I; en comparant la séquence de l’individu à une base de données. 
A)  Identification d’une espèce par codebarre
-Allez à l’adresse suivante : http://www.barcodinglife.org/
-puis dans l’onglet Identification
-copiez votre séquence 
	 
-cliquez sur Submit
-le résultat apparait en plusieurs sections :
	- Search Result:
Dans cette section, vous trouverez le nom de l’espèce à laquelle votre séquence correspond, ainsi que des informations relatives à cette espèce si elles existent dans ce site. 
	- Identification Summary:
Ici vous trouverez la probabilité d’identification de votre séquence en fonction du niveau taxonomique. 
-Similarity Scores of Top 99 Matches:
La similarité des 99 séquences les plus proches de votre séquence.
- TOP 20 Matches :
Les 20 entrées les plus proches de votre séquence. Vous trouverez leur identification complète ainsi que le pourcentage de similarité avec votre séquence et leur statut (publique ou privé). Toutes les séquences ne sont pas accessibles au grand public, car certaines font parties d’études qui ne sont pas encore publiées. Ces séquences privées existent mais vous ne pouvez accéder à d’autres informations.
- Sampling Sites For Top Hits (>98% Match):
Ici les points sur la carte représentent les sites où les séquences ayant plus de 98% de similarité avec la vôtre ont été échantillonnées.

En absence de référence
BOLD est une base de données mais qui est incomplète et il se peut que a référence n’existe pas. Dans ce cas un message no hit match va s’afficher.
BOLD vous suggère d’aller chercher dans une autre base de données BLAST.
- Cliquez sur Blast Sequence On Genbank. Votre séquence sera automatiquement rentrée sur le site de BLAST
-Laissez les paramètres par défaut. Ici la recherche va s’effectuer dans la base de données nucléotidique. 
-Cliquez sur Highly similar sequences (megablast). Cette option va vous permettre de rechercher les séquences les plus similaires.
-Cliquez sur BLAST
Après un certain temps (de l’ordre de quelques minutes maximum), votre résultat va s’afficher en plusieurs parties :
-Nucleotide Sequence 
Dans cette section, vous trouverez un résumé de votre requête : longueur de votre séquence, à quelle base de données vous l’avez comparé, ect..
- Graphic Summary
Ici les résultats les plus proches sont résumés sous forme de traits et de couleurs. La longueur des traits représente la longueur de la séquence trouvée dans BLAST. Les couleurs allant du noir au rouge symbolisent le score d’alignement entre la séquence trouvée et votre séquence. Quand c’est noir c’est faible (<40), quand c’est rouge c’est très similaire (>200). Si vous cliquez sur une ligne vous allez été redirigé à la description détaillée de la séquence.
- Descriptions
La description détaillée des résultats précédents. Vous trouverez pour chaque séquence son score d’alignement, sa couverture par rapport à votre séquence, sa e-value, le pourcentage d’identité et son numéro d’accession.
score d’alignement : quand deux séquences sont comparées, une valeur d’alignement (score) va être générée. Ce score va dépendre du nombre de mismatchs, de gaps et de matchs.
query-cover= couverture: vous pouvez comparée des séquences de différentes longueur, par exemple une plus longue et une plus courte. Votre séquence plus courte ne couvrira qu’un certain pourcentage de votre séquence plus longue.
	e-value: la probabilité d’observer ce résultat par hasard. Plus la e-value est petite et moins vous avez de chance d’obtenir ce résultat par hasard.
	Identity : pourcentage d’identité. Nombre de bases identiques entre votre séquence et celle comparée.
	Accession : Identifiant de Genbank.


1)	Pour chacune des séquences inconnues, indiquez s’il y avait un match dans BOLD.  Si vous avez dû aller dans genbank, quelle était l’espèce la plus proche ? Donnez les valeurs de couverture, la e-value.


B) Détermination du barcode gap chez Aurelia, Poecile, Salvelinus
Sur le site de BOLD SYSTEMS, allez dans Taxonomy, et rentrer le nom du genre recherché.
Cliquez sur ‘public data’ puis ‘combined’ ‘TSV’ (en rouge)
Allez rechercher le fichier txt téléchargé, copier le et coller le dans un fichier excel.
Copier les ‘process id’ dans la 1ère colonne
Retourner sur Bold, cliquez sur ‘workbench’, puis ‘record search’
Coller les données ‘process id’ et search
Cochez ‘public record’
Deroulez la page et cocher ‘all’ data
Allez dans taxon ID tree, puis cochez ‘Muscle’ dans alignement et dans ‘Colorize tree based on ‘Barcode ‘ : Bin
Puis cliquez sur ‘Build tree’ et visualiser le ficher pdf
Inspectez l’arbre coloré selon les Bin (individus regroupés par cluster de 2.2 %
2)	Est-ce-que l’arbre parvient à résoudre les différentes espèces ?
3)	Indiquez les espèces qui causent problème 
Allez dans barcode gap analysis

4)	Quelle est la distance moyenne intraspécifique ? Quelles sont les valeurs minimum et maximum ?
5)	Quelle est la distance à l’espèce voisine la plus proche ?
6)	Est-ce-qu’il y a un barcode gap ?
7)	Y-a-t’il des problèmes de classification ?






8)	Discutez des résultats que vous avez obtenus pour l’identification des espèce des trois différents genres. Pour Salvelinus, discutez de la séparation entre les malma et les alpinus. Quelles sont les principales limites associées à cette méthode pour la découverte de la biodiversité ? 




Neighbor-Joining: méthode de reconstruction phylogénétique se basant sur les distances évolutives. À partir des distances mesurées entre chaque séquence on part d’un arbre en étoile (i.e les séquences sont placées arbitrairement les unes par rapport aux autres). Puis, on calcule pour chaque paire de séquence la longueur des branches entre celles-ci. On retient la paire de séquence pour laquelle la longueur minimale des branches est obtenue; puis on les regroupe. On recommence avec les N-2 séquences restantes et jusqu’à avoir placé toutes les séquences.

