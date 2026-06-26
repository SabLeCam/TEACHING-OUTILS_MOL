# TD — Bio-informatique Marine : Exploration de BOLD
## Identifier les espèces bretonnes par leur code-barres ADN

**Niveau :** Seconde — Sciences Numériques et Technologie  
**Durée :** 1h30  
**Outils :** navigateur web, tableur (LibreOffice Calc ou Excel)

---

## Contexte

L'océan breton est l'un des plus riches d'Europe. Dans la rade de Brest, la mer d'Iroise ou le golfe de Gascogne se croisent des milliers d'espèces : du phytoplancton microscopique aux requins pélagiques.

Comment les scientifiques de l'IFREMER (basé à Brest) savent-ils *exactement* quelle espèce est présente dans un échantillon d'eau ? Grâce au **code-barres ADN**, et à la base de données mondiale **BOLD** (*Barcode of Life Data System*).

> **BOLD** : [https://www.boldsystems.org](https://www.boldsystems.org)  
> Contient plus de **10 millions de séquences** pour plus de **300 000 espèces**.

---

## Principe : qu'est-ce qu'un code-barres ADN ?

Tout comme un code-barres de supermarché identifie un produit, un fragment du gène **COI** (*Cytochrome Oxydase sous-unité I*, ~650 paires de bases) identifie une espèce animale. Pour les algues et bactéries, on utilise les gènes **18S rRNA** ou **16S rRNA**.

```
Espèce A : AACTTTATATTTTATTTTTGGAGCCTGAGCC...
Espèce B : AACTTTATATTTCATTTTTGGAGCATGAGCC...
              ↑ différence ici → espèce différente !
```

La ressemblance entre séquences s'appelle **similarité** (exprimée en %). Deux individus de la même espèce ont en général **> 97 % de similarité** sur le COI.

---

## Le jeu de données

Vous disposez du fichier **`bold_bretagne.tsv`**, un extrait fictif mais réaliste d'une base de données BOLD avec **20 spécimens** prélevés lors de campagnes IFREMER en Bretagne (2023).

Le fichier contient cinq groupes d'organismes :

| Groupe | Espèces | Code groupe |
|--------|---------|-------------|
| Phytoplancton | *Pseudo-nitzschia pungens*, *Alexandrium minutum*, *Micromonas pusilla*, *Synechococcus elongatus*, *Emiliania huxleyi* | PHYTO |
| Zooplancton | *Meganyctiphanes norvegica*, *Calanus helgolandicus*, *Sagitta elegans*, *Aequorea victoria*, *Liocarcinus holsatus* | ZOO |
| Bivalves | *Mytilus edulis*, *Crassostrea gigas*, *Ostrea edulis* | BIVALVE |
| Poissons | *Sardina pilchardus*, *Merluccius merluccius*, *Labrus bergylta*, *Solea solea* | FISH |
| Requins | *Lamna nasus*, *Prionace glauca*, *Scyliorhinus canicula* | ELASMOBRANCH |

---

## Partie 1 — Explorer le fichier (20 min)

### 1.1 Ouvrir le fichier TSV

1. Ouvrez **LibreOffice Calc** (ou Excel)
2. Faites **Fichier → Ouvrir** et sélectionnez `bold_bretagne.tsv`
3. Choisissez le séparateur **Tabulation**

> **Conseil :** Le fichier a beaucoup de colonnes. Utilisez **Ctrl+Fin** pour voir jusqu'où il s'étend.

### 1.2 Questions d'exploration

**Q1.** Combien y a-t-il de colonnes dans ce fichier ? Combien de lignes (hors en-tête) ?

**Q2.** Repérez les colonnes suivantes et notez leur position (lettre de colonne) :
- `species_name` → colonne ___
- `lat` / `lon` → colonnes ___ / ___
- `depth` → colonne ___
- `nucleotides` → colonne ___
- `markercode` → colonne ___

**Q3.** Quel marqueur génétique est utilisé pour les poissons ? Pour le phytoplancton ?

**Q4.** Quelle est la profondeur de prélèvement la plus grande ? À quelle espèce correspond-elle ?

---

## Partie 2 — Requêtes dans le tableur (25 min)

### 2.1 Filtrage par groupe

En utilisant le **filtre automatique** (Données → AutoFiltre) :

**Q5.** Combien d'espèces appartiennent au groupe `ZOO` (zooplancton) ?

**Q6.** Listez toutes les espèces prélevées dans la **Rade de Brest** (colonne `region` ou `exactsite`).

**Q7.** Quelles espèces ont été collectées en **profondeur > 100 m** ?

### 2.2 Formule de localisation

Sélectionnez les 20 spécimens et regardez les colonnes `lat` et `lon`.

**Q8.** Toutes les coordonnées GPS sont-elles en Bretagne ? Vérifiez sur [Google Maps](https://maps.google.com) en tapant `47.29, -4.78` (format: latitude, longitude).

**Q9.** Quelle est l'espèce prélevée le plus au large (longitude la plus négative) ?

### 2.3 Tableau de synthèse

Complétez ce tableau à partir du fichier :

| Espèce | Groupe | Profondeur (m) | Lieu exact | Mois |
|--------|--------|---------------|-----------|------|
| *Pseudo-nitzschia pungens* | | | | |
| *Calanus helgolandicus* | | | | |
| *Mytilus edulis* | | | | |
| *Sardina pilchardus* | | | | |
| *Lamna nasus* | | | | |

---

## Partie 3 — BOLD en ligne (25 min)

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

**Q12.** Cette espèce figure-t-elle dans notre fichier `bold_bretagne.tsv` ? Quel est son processid ?

### 3.2 Explorer une espèce sur BOLD

Choisissez **une espèce** du fichier TSV (au choix) et cherchez-la sur BOLD :
[https://www.boldsystems.org](https://www.boldsystems.org)

**Q13.** Dans combien de pays cette espèce a-t-elle été barcodée ?

**Q14.** Existe-t-il plusieurs « BIN » (Barcode Index Numbers) pour cette espèce ? Qu'est-ce que cela signifie ?

---

## Partie 4 — Réflexion (20 min)

### 4.1 Chaîne alimentaire

En vous aidant du fichier et de vos connaissances :

**Q15.** Construisez une chaîne alimentaire simplifiée avec au moins 4 niveaux en utilisant uniquement des espèces présentes dans le fichier.

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

Exportez les colonnes `species_name`, `lat`, `lon` du fichier TSV en CSV, puis importez-les dans [Google My Maps](https://mymaps.google.com) ou [uMap](https://umap.openstreetmap.fr).

Vous obtiendrez une **carte de distribution** des 20 espèces bretonnes !

### Comparaison de séquences

Sur le site [MUSCLE Alignment](https://www.ebi.ac.uk/jdispatcher/msa/muscle), alignez les séquences COI de *Mytilus edulis* (SEQ011), *Crassostrea gigas* (SEQ012) et *Ostrea edulis* (SEQ013) du fichier.

- Où sont les différences entre les trois séquences ?
- Ces trois espèces sont-elles proches ou éloignées génétiquement ?

---

## Corrigé rapide (pour l'enseignant)

| Question | Réponse |
|----------|---------|
| Q3 | Poissons : COI-5P ; Phyto : 18S-rRNA, 16S-rRNA |
| Q4 | 400 m — *Meganyctiphanes norvegica* (krill, migrateur vertical) |
| Q5 | 5 espèces zooplancton |
| Q9 | *Lamna nasus* et *Prionace glauca* — lon = -8.0000 (Atlantique large) |
| Q10–11 | La séquence correspond à *Mytilus edulis* (~99 % similarité) |
| Q12 | Oui — BOLD_BRE011 |
| Q16 | *Crassostrea gigas* — colonne `notes` : "invasive in Brittany" |

---

## Ressources

| Ressource | Lien |
|-----------|------|
| BOLD Systems | https://www.boldsystems.org |
| IFREMER | https://www.ifremer.fr |
| Tara Oceans | https://fondation.tara.org |
| NCBI BLAST | https://blast.ncbi.nlm.nih.gov |
| GBIF (cartes de distribution) | https://www.gbif.org |
| Liste rouge IUCN | https://www.iucnredlist.org |

---

*Document produit dans le cadre du cours de Sciences Numériques et Technologie, Seconde — IFREMER / Académie de Rennes*
