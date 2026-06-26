# TD — Bio-informatique Marine : Exploration de BOLD (stage seconde)
## Identifier les espèces bretonnes par leur code-barres ADN



---

## Contexte



Grâce au **code-barres ADN**, et à la base de données mondiale **BOLD** (*Barcode of Life Data System*), identifier les espèces détéctée dans un échantillon d'eau de mer.

> **BOLD** : [https://www.boldsystems.org](https://www.boldsystems.org)  
> Contient plus de **10 millions de séquences** pour plus de **300 000 espèces**.

---

## Principe : qu'est-ce qu'un code-barres ADN ?

Tout comme un code-barres de supermarché identifie un produit, un fragment du gène **COI** (*Cytochrome Oxydase sous-unité I*, ~650 paires de bases) identifie une espèce animale.

```
Espèce A : AACTTTATATTTTATTTTTGGAGCCTGAGCC...
Espèce B : AACTTTATATTTCATTTTTGGAGCATGAGCC...
              ↑ différence ici → espèce différente !
```

La ressemblance entre séquences s'appelle **similarité** (exprimée en %). Deux individus de la même espèce ont en général **> 97 % de similarité** sur le COI.

---

## Le jeu de données

Vous disposez du fichier **`bold_bretagne.tsv`**, un extrait fictif de séquences rétrouvées dans not échantillon d'eau.



---

## Partie 1 — Explorer le fichier (20 min)

### 1.1 Ouvrir le fichier TSV

1. Ouvrez **LibreOffice Calc** (ou Excel)
2. Faites **Fichier → Ouvrir** et sélectionnez `bold_bretagne.tsv`
3. Choisissez le séparateur **Tabulation**


### 1.2 Questions d'exploration

**Q1.** Combien y a-t-il de colonnes dans ce fichier ? Combien de lignes (hors en-tête) ?

**Q2.** Repérez les colonnes suivantes et notez leur position (lettre de colonne) :

- `nucleotides` → colonne ___
- `markercode` → colonne ___

**Q3.** Quel marqueur génétique est utilisé pour les poissons ? Pour le phytoplancton ?


## Partie 2 — BOLD en ligne (25 min)

### 3.1 Identifier une séquence inconnue

Voici une séquence mystère trouvée dans un échantillon d'eau de la rade de Brest :

```
AACTTTATATTTCATTTTTGGAGCTTGAGCTGGAATAGTAGGTACAGCTTTAAGACTCCTAATTCGAGCCGAATTAGGACAACCCGGTGCACTAATTGGAGAC
```

**Q10.** Rendez-vous sur BOLD : [https://www.boldsystems.org/index.php/IDS_OpenIdEngine](https://www.boldsystems.org/index.php/IDS_OpenIdEngine)

- Copiez-collez la séquence ci-dessus
- Sélectionnez le moteur **COI**
- Cliquez sur **Identify**

**Q11.** Quelle espèce BOLD identifie-t-il ? Quel est le pourcentage de similarité ?

**Q12.** recommencez avec toutes les séquences de `bold_bretagne.tsv` ? Quelles espèces retrouve-t-on dans notre échantillon.

### 3.2 Explorer une espèce sur BOLD

Choisissez **une espèce** du fichier résultat (au choix) et cherchez-la sur BOLD :
[https://www.boldsystems.org](https://www.boldsystems.org)

**Q13.** Dans combien de pays cette espèce a-t-elle été barcodée ?

**Q14.** Existe-t-il plusieurs « BIN » (Barcode Index Numbers) pour cette espèce ? Qu'est-ce que cela signifie ?

---

## Partie 4 — Réflexion (20 min)

### 4.1 Chaîne alimentaire

En vous aidant du fichier et de vos connaissances :

**Q15.** Construisez une chaîne alimentaire simplifiée  en utilisant uniquement des espèces présentes dans le fichier.

```
Producteur primaire → ... → ... → Prédateur supérieur
```

### 4.2 Espèce invasive

**Q16.** Une des espèces du fichier est considérée comme **invasive** en Bretagne. Laquelle ? Comment le savez-vous en regardant le fichier ?

### 4.3 Conservation

**Q17.** L'espèce *Lamna nasus* (requin taupe commun) est classée **Vulnérable** sur la liste rouge IUCN. Pourtant, elle figure dans notre fichier comme "bycatch" (prise accessoire). 

- Qu'est-ce que le bycatch ?
- Comment la bio-informatique (BOLD) peut-elle aider à surveiller les populations de cette espèce ?

---

## Bonus — Pour aller plus loin

### Visualiser les données sur une carte

Exportez les colonnes `species_name`, `lat`, `lon` du fichier TSV complete en CSV, puis importez-les dans [Google My Maps](https://mymaps.google.com) ou [uMap](https://umap.openstreetmap.fr).

Vous obtiendrez une **carte de distribution** des 20 espèces bretonnes !

### Comparaison de séquences

Sur le site [MUSCLE Alignment](https://www.ebi.ac.uk/jdispatcher/msa/muscle), alignez les séquences COI de *Mytilus edulis* (SEQ011), *Crassostrea gigas* (SEQ012) et *Ostrea edulis* (SEQ013) du fichier.

- Où sont les différences entre les trois séquences ?
- Ces trois espèces sont-elles proches ou éloignées génétiquement ?
