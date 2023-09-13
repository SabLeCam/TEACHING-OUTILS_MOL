# TP DIVERSITE GENETIQUE ET ECART A HARDY WEINBERG

(Sabrina Le Cam et France Dufresne, d'après un document produit par Tom Jenkins source https://tomjenkins.netlify.app/tutorials/r-popgen-getting-started/)

Dans ce TP, nous allons apprendre les bases de l'analyses de données de génétique des populations. Tout d'abord, explorer les données et calculer les différents indices de la diverisité génétique abordés en cours.

Le jeu de données est issu des articles suivants
[Jenkins TL, Ellis CD, Triantafyllidis A, Stevens JR (2019). Single nucleotide polymorphisms reveal a genetic cline across the north‐east Atlantic and enable powerful population assignment in the European lobster. Evolutionary Applications 12, 1881–1899.](https://onlinelibrary.wiley.com/doi/10.1111/eva.12849)
[Holland LP, Jenkins TL, Stevens JR (2017). Contrasting patterns of population structure and gene flow facilitate exploration of connectivity in two widely distributed temperate octocorals. Heredity 119, 35–48.](https://www.nature.com/articles/hdy201714)

Calculer les indices de diversité génétique de bases:

  - Polymorphisme/fréquence allèliques
  - Diversité génetique (hétérozygotie  attendue et observée)
  - Test d'écart à l'équilibre de Hardy Weinberg


 ## Ressources et packages nécessaires

Nous allons effectuer ces analyses sous R. 
Installez les packages suivants:

```r
install.packages("adegenet")
install.packages("dplyr")
install.packages("hierfstat")
install.packages("(reshape2")
install.packages("ggplot2")
install.packages("(RColorBrewer")
install.packages("scales")
install.packages("poppr")
install.packages("dartR")
```


Ouverture des librairies
```r
library("adegenet")
library("poppr")
library("dplyr")
library("hierfstat")
library("reshape2")
library("ggplot2")
library("RColorBrewer")
library("scales")
library("poppr")
library("dartR")
```
## Données

Nous allons analyser une jeu de données de SNP des 1278 individus (homards) issus de 38 populations et génotypés à 79 microsatellites.
Nous allons importer une fichier de format 'csv' en tant qu'objet "genind"

Vue du fichier csv (NB: dans ce fichier il y a plus que 79 SNP):

```
Site,ID,Locus,Genotype
Ale,Ale04,3441,GG
Ale,Ale04,4173,NA
Ale,Ale04,6157,NA
Ale,Ale04,7502,NA
Ale,Ale04,7892,TT
Ale,Ale04,8953,TT
Ale,Ale04,9441,NA
...
Ale,Ale16,29889,AA
Ale,Ale16,30339,NA
Ale,Ale16,31462,AA
Ale,Ale16,31618,GG
Ale,Ale16,31967,NA
Ale,Ale16,31979,NA
...

Sar13,Sr5,3441,GG
Sar13,Sr5,4173,CT
Sar13,Sr5,6157,GC
Sar13,Sr5,7502,TT
Sar13,Sr5,7892,TT
Sar13,Sr5,8953,GT
```



## ouverture de fichier

Vérifier que les données sont bien représentées dans l'objet ```data.frame```
```r
#setwd("C:/Users/marie/Documents/Deuxieme_cycleUQAR/")
lobster = read.csv("Lobster_SNP_Genotypes.csv")
str(lobster)
## 'data.frame':    125280 obs. of  4 variables:
##  $ Site    : chr  "Ale" "Ale" "Ale" "Ale" ...
##  $ ID      : chr  "Ale04" "Ale04" "Ale04" "Ale04" ...
##  $ Locus   : int  3441 4173 6157 7502 7892 8953 9441 11071 11183 11291 ...
##  $ Genotype: chr  "GG" NA NA NA ...
```

Convertir le ```data.frame```pour qu'un rang corresponde à un individu et une colonne à une locus, puis une colonne pour les ID et les noms de sites.
```r
lobster_wide = reshape(lobster, idvar = c("ID","Site"), timevar = "Locus", direction = "wide", sep = "")
## Warning in reshapeWide(data, idvar = idvar, timevar = timevar, varying =
## varying, : multiple rows match for Locus=3441: first taken

# Remove "Genotype" from column names
colnames(lobster_wide) = gsub("Genotype", "", colnames(lobster_wide))
```

Utilisez seulement une partie des données de Jenkins et al. 2019
```r
# Subset genotypes
snpgeno = lobster_wide[ , 3:ncol(lobster_wide)]

# Keep only SNP loci used in Jenkins et al. 2019
snps_to_remove = c("25580","32362","41521","53889","65376","8953","21197","15531","22740","28357","33066","51507","53052","53263","21880","22323","22365")
snpgeno = snpgeno[ , !colnames(snpgeno) %in% snps_to_remove]
```

Création de vecteurs pour les individus et les sites
```r
ind = as.character(lobster_wide$ID) # individual ID
site = as.character(lobster_wide$Site) # site ID
```

Convertir ```data.frame``` en objet genind object. Vérifiez que les génotypes pour les cinq premiers individus et loci sont bons.
```r
lobster_gen = df2genind(snpgeno, ploidy = 2, ind.names = ind, pop = site, sep = "")
lobster_gen$tab[1:5, 1:10]
##       3441.G 3441.A 4173.C 4173.T 6157.G 6157.C 7502.T 7502.C 7892.T 7892.A
## Ale04      2      0     NA     NA     NA     NA     NA     NA      2      0
## Ale05      1      1      2      0      1      1      2      0      1      1
## Ale06     NA     NA      2      0      2      0     NA     NA      2      0
## Ale08     NA     NA      0      2      2      0      2      0     NA     NA
## Ale13      2      0     NA     NA      2      0     NA     NA      2      0
```
Imprimer information de bases pour les objets genind

```r
lobster_gen
```
```r

## /// GENIND OBJECT /////////
## 
##  // 1,305 individuals; 79 loci; 158 alleles; size: 945.1 Kb
## 
##  // Basic content
##    @tab:  1305 x 158 matrix of allele counts
##    @loc.n.all: number of alleles per locus (range: 2-2)
##    @loc.fac: locus factor for the 158 columns of @tab
##    @all.names: list of allele names for each locus
##    @ploidy: ploidy of each individual  (range: 2-2)
##    @type:  codom
##    @call: df2genind(X = snpgeno, sep = "", ind.names = ind, pop = site, 
##     ploidy = 2)
```
Explorer les différentes infos disponibles:
```r
lobster_gen@loc.n.all
lobster_gen@pop
```
Importer le fichier de g?notypes Microsatellites du corail
Importer le fichier genepop et convertir en objet genind. V?rifiez que les g?notypes au locus  Ever002 pour trois individus choisis au hasard sont conformes.

seafan_gen = import2genind("Pinkseafan_13MicrosatLoci.gen", ncode = 3, quiet = TRUE)
set.seed(1)
tab(seafan_gen[loc = "Ever002"])[runif(3, 1, nInd(seafan_gen)), ]
##       Ever002.114 Ever002.117 Ever002.109 Ever002.105 Ever002.121
## Far10           2           0           0           0           0
## Han36           2           0           0           0           0
## Moh5            1           1           0           0           0

# Imprimer l'info de base de cet objet genind

seafan_gen
## /// GENIND OBJECT /////////
## 
##  // 877 individuals; 13 loci; 114 alleles; size: 478.2 Kb
## 
##  // Basic content
##    @tab:  877 x 114 matrix of allele counts
##    @loc.n.all: number of alleles per locus (range: 2-18)
##    @loc.fac: locus factor for the 114 columns of @tab
##    @all.names: list of allele names for each locus
##    @ploidy: ploidy of each individual  (range: 2-2)
##    @type:  codom
##    @call: read.genepop(file = file, ncode = 3, quiet = quiet)
## 




Mettre ? jour les identifiants des sites pour que les codes des sites soient utilis?s plut?t que les identifiants des individus

# Use gsub to extract only letters from a vector
popNames(seafan_gen) = gsub("[^a-zA-Z]", "", popNames(seafan_gen))
popNames(seafan_gen)
##  [1] "ArmI"   "ArmII"  "ArmIII" "Bla"    "Bov"    "Bre"    "Far"    "Fla"   
##  [9] "Han"    "Lao"    "Lio"    "Lun"    "Men"    "Mew"    "Moh"    "PorI"  
## [17] "PorII"  "Rag"    "RosI"   "RosII"  "Sko"    "Thu"    "Vol"    "Wtn"




Imprimer le nombrer d'all?les par locus
table(lobster_gen$loc.fac)
## 
##  3441  4173  6157  7502  7892  9441 11071 11183 11291 12971 14047 14742 15109 
##     2     2     2     2     2     2     2     2     2     2     2     2     2 
## 15128 15435 15581 18512 18652 19266 19460 20354 23146 23447 23481 23677 23787 
##     2     2     2     2     2     2     2     2     2     2     2     2     2 
## 24020 25229 25608 27329 29410 29801 29889 30339 31462 31618 31967 31979 32358 
##     2     2     2     2     2     2     2     2     2     2     2     2     2 
## 32435 33784 34443 34818 35584 36910 39107 39876 42395 42529 42821 44670 45154 
##     2     2     2     2     2     2     2     2     2     2     2     2     2 
## 45217 51159 53314 53720 53935 54240 54762 55111 55142 55564 56423 56785 57131 
##     2     2     2     2     2     2     2     2     2     2     2     2     2 
## 57989 58053 59503 59586 59967 60546 63140 63267 63581 63605 63771 63798 65064 
##     2     2     2     2     2     2     2     2     2     2     2     2     2 
## 65576 
##     2

table(seafan_gen$loc.fac)
## 
## Ever001 Ever002 Ever003 Ever004 Ever006 Ever007 Ever008 Ever010 Ever011 Ever013 
##      15       5       5      13      18      11       2       9       5      14 
## Ever014 
##       9

Imprimer la taille d'?chantillon
summary(lobster_gen$pop)
##   Ale   Ber   Brd   Cor   Cro   Eye   Flo   Gul   Heb   Hel   Hoo Idr16 Idr17 
##    28    33    36    32    35    26    36    35    36    35    36    32    29 
##   Iom   Ios   Jer   Kav   Kil   Laz   Loo   Lyn   Lys   Mul   Oos   Ork   Pad 
##    35    36    36    36    35     5    36    34    36    36    40    36    36 
##   Pem Sar13 Sar17   Sbs   She   Sin   Sky   Sul   Tar   The   Tor   Tro   Ven 
##    36     7    15    36    36    35    37    36     5    36    37    17    36 
##   Vig 
##    36

summary(seafan_gen$pop)
##   ArmI  ArmII ArmIII    Bla    Bov    Bre    Far    Fla    Han    Lao    Lio 
##     26     43     39     29     38     43     44     23     36     36     22 
##    Lun    Men    Mew    Moh   PorI  PorII    Rag   RosI  RosII    Sko    Thu 
##     21     43     44     28     39     34     42     39     35     39     48 
##    Vol    Wtn 
##     23     43

Imprimer le nombre d'all?les priv?s par site par locus

private_alleles(seafan_gen) %>% apply(MARGIN = 1, FUN = sum)
##   ArmI  ArmII ArmIII    Bla    Bov    Bre    Far    Fla    Han    Lao    Lio 
##      1      1      0      0      0      0      1      0      0      1      0 
##    Lun    Men    Mew    Moh   PorI  PorII    Rag   RosI  RosII    Sko    Thu 
##      1      1      1      0      1      4      0      2      1      0      2 
##    Vol    Wtn 
##      0      0


Imprimer la richesse all?lique moyenne par sites pour l'ensemble des loci

allelic.richness(seafan_gen,min.n=NULL,diploid=TRUE)
##   ArmI  ArmII ArmIII    Bla    Bov    Bre    Far    Fla    Han    Lao    Lio 
##  2.771  2.720  2.748  2.635  2.784  2.837  2.807  2.698  3.030  2.809  2.957 
##    Lun    Men    Mew    Moh   PorI  PorII    Rag   RosI  RosII    Sko    Thu 
##  2.915  2.824  2.895  2.791  2.900  2.833  2.895  2.831  2.966  2.905  2.650 
##    Vol    Wtn 
##  2.767  3.032

Calculer l'h?t?rozygote par site

# Calculer les stats de base en utilisant hierfstat
basic_lobster = basic.stats(lobster_gen, diploid = TRUE)
basic_seafan = basic.stats(seafan_gen, diploid = TRUE)


# H?t?rozygotie moyenne observ?e par site 
Ho_lobster = apply(basic_lobster$Ho, MARGIN = 2, FUN = mean, na.rm = TRUE) %>%
  round(digits = 2)
Ho_lobster
##   Ale   Ber   Brd   Cor   Cro   Eye   Flo   Gul   Heb   Hel   Hoo Idr16 Idr17 
##  0.32  0.36  0.37  0.38  0.37  0.37  0.35  0.38  0.39  0.35  0.39  0.39  0.39 
##   Iom   Ios   Jer   Kav   Kil   Laz   Loo   Lyn   Lys   Mul   Oos   Ork   Pad 
##  0.39  0.39  0.38  0.37  0.38  0.38  0.39  0.40  0.34  0.37  0.32  0.36  0.37 
##   Pem Sar13 Sar17   Sbs   She   Sin   Sky   Sul   Tar   The   Tor   Tro   Ven 
##  0.38  0.32  0.34  0.38  0.37  0.35  0.33  0.36  0.42  0.33  0.33  0.33  0.39 
##   Vig 
##  0.39

Ho_seafan = apply(basic_seafan$Ho, MARGIN = 2, FUN = mean, na.rm = TRUE) %>%
  round(digits = 2)
Ho_seafan
##   ArmI  ArmII ArmIII    Bla    Bov    Bre    Far    Fla    Han    Lao    Lio 
##   0.41   0.45   0.44   0.44   0.45   0.46   0.43   0.47   0.50   0.45   0.50 
##    Lun    Men    Mew    Moh   PorI  PorII    Rag   RosI  RosII    Sko    Thu 
##   0.49   0.47   0.51   0.41   0.44   0.45   0.51   0.49   0.49   0.53   0.40 
##    Vol    Wtn 
##   0.47   0.50


# H?t?rozygotie moyenne observ?e par site 
He_lobster = apply(basic_lobster$Hs, MARGIN = 2, FUN = mean, na.rm = TRUE) %>%
  round(digits = 2)
He_lobster
##   Ale   Ber   Brd   Cor   Cro   Eye   Flo   Gul   Heb   Hel   Hoo Idr16 Idr17 
##  0.34  0.36  0.37  0.39  0.37  0.37  0.35  0.36  0.38  0.35  0.39  0.39  0.39 
##   Iom   Ios   Jer   Kav   Kil   Laz   Loo   Lyn   Lys   Mul   Oos   Ork   Pad 
##  0.39  0.39  0.38  0.36  0.38  0.34  0.37  0.39  0.35  0.38  0.33  0.37  0.37 
##   Pem Sar13 Sar17   Sbs   She   Sin   Sky   Sul   Tar   The   Tor   Tro   Ven 
##  0.38  0.32  0.35  0.37  0.37  0.35  0.33  0.37  0.36  0.34  0.33  0.36  0.38 
##   Vig 
##  0.39

He_seafan = apply(basic_seafan$Hs, MARGIN = 2, FUN = mean, na.rm = TRUE) %>%

#Visualiser l'h?t?rozygote par site
# Create a data.frame of site names, Ho and He and then convert to long format
Het_lobster_df = data.frame(Site = names(Ho_lobster), Ho = Ho_lobster, He = He_lobster) %>%
  melt(id.vars = "Site")
Het_seafan_df = data.frame(Site = names(Ho_seafan), Ho = Ho_seafan, He = He_seafan) %>%
  melt(id.vars = "Site")

# Custom theme for ggplot2
custom_theme = theme(
  axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, face = "bold"),
  axis.text.y = element_text(size = 10),
  axis.title.y = element_text(size = 12),
  axis.title.x = element_blank(),
  axis.line.y = element_line(size = 0.5),
  legend.title = element_blank(),
  legend.text = element_text(size = 12),
  panel.grid = element_blank(),
  panel.background = element_blank(),
  plot.title = element_text(hjust = 0.5, size = 15, face="bold")
  )
# Italic label
hetlab.o = expression(italic("H")[o])
hetlab.e = expression(italic("H")[e])

# Barplot pour l'h?t?rozygotie du homard
ggplot(data = Het_lobster_df, aes(x = Site, y = value, fill = variable))+
  geom_bar(stat = "identity", position = position_dodge(width = 0.6), colour = "black")+
  scale_y_continuous(expand = c(0,0), limits = c(0,0.50))+
  scale_fill_manual(values = c("royalblue", "#bdbdbd"), labels = c(hetlab.o, hetlab.e))+
  ylab("Heterozygosity")+
  ggtitle("European lobster")+
  custom_theme

# Barplot pour l'h?t?rozygotie du corail
ggplot(data = Het_seafan_df, aes(x = Site, y = value, fill = variable))+
  geom_bar(stat = "identity", position = "dodge", colour = "black")+
  scale_y_continuous(expand = c(0,0), limits = c(0,0.60), breaks = c(0, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60))+
  scale_fill_manual(values = c("pink", "#bdbdbd"), labels = c(hetlab.o, hetlab.e))+
  ylab("Heterozygosity")+
  ggtitle("Pink sea fan")+
  custom_theme

Calculer le coefficient de consanguinit?
Calculerr le FIS moyen par site

# European lobster
apply(basic_lobster$Fis, MARGIN = 2, FUN = mean, na.rm = TRUE) %>%
  round(digits = 3)
##    Ale    Ber    Brd    Cor    Cro    Eye    Flo    Gul    Heb    Hel    Hoo 
##  0.057  0.003  0.003  0.021 -0.006 -0.004  0.005 -0.044 -0.034  0.013 -0.016 
##  Idr16  Idr17    Iom    Ios    Jer    Kav    Kil    Laz    Loo    Lyn    Lys 
## -0.007 -0.001 -0.024 -0.007  0.004 -0.024 -0.016 -0.115 -0.043 -0.032  0.018 
##    Mul    Oos    Ork    Pad    Pem  Sar13  Sar17    Sbs    She    Sin    Sky 
##  0.040  0.023  0.017 -0.010 -0.004 -0.009  0.018 -0.017 -0.006 -0.013  0.006 
##    Sul    Tar    The    Tor    Tro    Ven    Vig 
##  0.033 -0.153  0.029  0.010  0.066 -0.024  0.013

# Pink sea fan
apply(basic_seafan$Fis, MARGIN = 2, FUN = mean, na.rm = TRUE) %>%
  round(digits = 3)
##   ArmI  ArmII ArmIII    Bla    Bov    Bre    Far    Fla    Han    Lao    Lio 
##  0.166  0.085  0.076 -0.006  0.075  0.039  0.116  0.014  0.042  0.064  0.029 
##    Lun    Men    Mew    Moh   PorI  PorII    Rag   RosI  RosII    Sko    Thu 
##  0.057  0.067  0.030  0.153  0.137  0.089  0.010  0.048  0.077  0.013  0.056 
##    Vol    Wtn 
##  0.057  0.058


        )