# 🏭 Distillation Multi-étagée

Programme Fortran 90 simulant en **régime permanent** une colonne de distillation biphasique à NP plateaux théoriques et NC constituants, avec condenseur partiel et bouilleur. Le système MESH complet est résolu par **Newton-Raphson amorti**, au choix par un solveur matriciel plein (**MRSL01**) ou par un solveur tridiagonal par blocs (**MRSL21**). Le comportement thermodynamique est au choix **idéal** (γᵢ = 1) ou **non idéal** (modèle **NRTL**) les quatre combinaisons sont opérationnelles.

![Fortran](https://img.shields.io/badge/fortran-90-green)
![Code::Blocks](https://img.shields.io/badge/IDE-Code%3A%3ABlocks-blue)
![gfortran](https://img.shields.io/badge/compiler-gfortran-orange)
![Projet](https://img.shields.io/badge/projet-acad%C3%A9mique-lightgrey)

## 🎯 Fonctionnalités

- **Système MESH rigoureux** : bilans de matière (total et partiels), équilibre liquide-vapeur, sommation des fractions, bilan enthalpique — soit **2·NC + 3 équations par plateau**
- **Deux modèles thermodynamiques** : solution idéale (γᵢ = 1) et solution non idéale (**NRTL** à coefficients d'interaction binaire)
- **Deux solveurs linéaires au choix** : MRSL01 (matrice pleine, pivot partiel) et MRSL21 (tridiagonal par blocs, exploitant la structure de la colonne)
- **Jacobienne numérique** par perturbation relative (δ = 10⁻³), assemblée en matrice pleine ou en blocs A/B/C selon le solveur
- **Initialisation physique** : profil de température interpolé entre température de bulle et température de rosée de l'alimentation, compositions et débits initialisés par bilan global
- **Colonne entièrement paramétrable** : nombre de plateaux, nombre de constituants, alimentations multiples, soutirages liquides et vapeur, puissances du bouilleur et du condenseur
- **Données externalisées** : aucune recompilation nécessaire pour changer de colonne ou de mélange, tout passe par des fichiers `.txt`
- **Résultats exportés en CSV** : profils de débits, compositions et températures, plus le vecteur des résidus finaux
- **Architecture modulaire** : modules de données (dimensions, thermo, opératoires, indices, messages) et sous-programmes de calcul séparés

## 📐 Modèle

Chaque plateau *j* porte **2·NC + 3 inconnues** : Vⱼ, yⱼ,ᵢ (NC), Tⱼ, xⱼ,ᵢ (NC), Lⱼ et autant d'équations :

```
Bilan matiere total     : L(j-1) - (L(j)+Dliq(j)) - (V(j)+Dvap(j)) + V(j+1) + F(j) = 0

Bilan matiere partiel   : L(j-1)·x(j-1,i) - (L(j)+Dliq(j))·x(j,i)
                          - (V(j)+Dvap(j))·y(j,i) + V(j+1)·y(j+1,i) + F(j)·z(j,i) = 0

Bilan enthalpique       : L(j-1)·Hliq(j-1) - (L(j)+Dliq(j))·Hliq(j)
                          - (V(j)+Dvap(j))·Hvap(j) + V(j+1)·Hvap(j+1) + F(j)·Halim(j) ± Q(j) = 0

Equilibre L/V           : y(j,i) - K(j,i)·x(j,i) = 0

Sommation               : Sum x(j,i) - Sum y(j,i) = 0
```

avec

```
K(i) = gamma(i) · Psat(i) / P                  (coefficient de partage)

Psat(i) = exp[ A + B/T + C·ln(T) + D·T^E ]     (Antoine etendue, en Pa)

gamma(i) = 1                                   (mode 1 : melange ideal)
gamma(i) = NRTL(x, T, Aij, Cij)                (mode 2 : melange non ideal)
```

Le modèle **NRTL** utilise τᵢⱼ = Cᵢⱼ / (R·T) et Gᵢⱼ = exp(−αᵢⱼ·τᵢⱼ), où `Cij` porte les énergies d'interaction binaire et `Aij` les paramètres de non-aléatoricité αᵢⱼ. En mode non idéal, l'enthalpie liquide inclut l'**enthalpie d'excès**, estimée par dérivation numérique de ln γᵢ par rapport à T.

**Résolution.** Le système F(Y) = 0, de taille NP·(2·NC+3), est résolu par **Newton-Raphson amorti** :

```
J(Y_k) · dY = F(Y_k)        puis        Y_(k+1) = Y_k - alpha · dY,   alpha = 0,1
```

L'amortissement α = 0,1 stabilise la marche au prix de la vitesse : la norme du résidu décroît d'un facteur (1 − α) par itération, d'où une convergence en quelques centaines d'itérations. Le critère d'arrêt est ‖F‖₂ < 10⁻⁶.

La jacobienne est calculée **numériquement**, par perturbation relative de chaque inconnue (Y → Y·(1+δ), δ = 10⁻³). Selon le solveur choisi, elle est assemblée soit en matrice pleine (MRSL01), soit sous forme de trois tableaux de blocs A, B, C exploitant le fait qu'un plateau n'échange qu'avec ses voisins immédiats (MRSL21).

**Système fourni en exemple** : mélange ternaire **Acétone – Benzène – Chloroforme** sous 101 325 Pa, sur 30 plateaux, alimenté au plateau 14. Le couple acétone/chloroforme forme un **azéotrope**, ce qui rend le mode NRTL indispensable à une description réaliste.

## 🚀 Démarrage rapide

### Prérequis

- **Code::Blocks avec MinGW** à télécharger depuis <https://www.codeblocks.org/downloads/binaries/> en choisissant impérativement une version dont le nom contient `mingw` (elle embarque **gfortran**). La version utilisée pour le projet est `codeblocks-25.03mingw-setup.exe`.

### Installation de Code::Blocks et réglage du compilateur

Ces étapes sont détaillées, captures à l'appui, dans [doc/installation.pdf](doc/installation.pdf).

```
1. Telecharger et lancer codeblocks-25.03mingw-setup.exe
   Installation standard, chemin par defaut : C:\Program Files\CodeBlocks

2. Ouvrir Code::Blocks puis aller dans :
   Settings > Compiler... > Toolchain executables

3. Selected compiler          : GNU Fortran Compiler (default)

4. Compiler's installation directory :
   C:\Program Files\CodeBlocks\MinGW        (bouton "Auto-detect" possible)

5. Onglet "Program Files", verifier que l'on a bien :
   C compiler              -> gfortran.exe
   C++ compiler            -> gfortran.exe
   Linker for dynamic libs -> gfortran.exe
   Linker for static libs  -> ar.exe
   Debugger                -> GDB/CDB debugger : Default
   Resource compiler       -> windres.exe
   Make program            -> gfortran.exe

6. Valider par OK
```

### Ouverture et exécution du projet

```
1. Cloner ou telecharger le depot
   git clone https://github.com/asbgra341/distillation-multi-etagee.git

2. File > Open...  puis selectionner :
   DistillationMultiEtagee/Distillation.cbp

3. Compiler et executer :
   Build > Build and run   (ou F9)
```

> ⚠️ Le programme lit ses fichiers de données avec des **chemins relatifs**. Le répertoire d'exécution doit donc être `DistllationMultiEtagee/` (c'est le comportement par défaut de Code::Blocks). Si vous lancez l'exécutable `bin/Debug/Distillation.exe` à la main, placez-vous dans le dossier `DistllationMultiEtagee/` ou copiez les fichiers `Donnees*.txt` à côté de lui.

> ⚠️ Après avoir modifié un fichier source, utilisez **Build > Rebuild (Ctrl + F11)** plutôt qu'un simple *Build*. Une reconstruction partielle peut relier d'anciens fichiers objets et produire un exécutable qui ne correspond plus aux sources.

### Compilation en ligne de commande (alternative)

L'ordre de compilation est imposé : les modules d'abord, le programme principal en dernier.

```bash
cd DistllationMultiEtagee
gfortran -Wall -O2 -o Distillation \
  mod_DIMENSIONS.f90 mod_THERMOS.f90 mod_OPERATOIRES.f90 mod_INDICES.f90 mod_MESSAGES.f90 \
  Affecte_Dimensions.f90 Affecte_Thermos.f90 Affecte_Operatoires.f90 Affecte_Indices.f90 \
  Thermo_Psaturation.f90 Thermo_Calculgama.f90 Thermo_CalculmKi.f90 \
  Thermo_CalculTbulle.f90 Thermo_CalculTrosee.f90 \
  Thermo_EnthalpieHliq.f90 Thermo_EnthalpieHvap.f90 \
  Programme_Initialisation.f90 Programme_Residus.f90 Programme_Norme.f90 \
  Programme_Jacobien01.f90 Programme_Jacobien21.f90 \
  MRSL01.f MRSL21.f Programme_Principale.f90
```

### Utilisation

Le programme pose deux questions : le **modèle thermodynamique**, puis le **solveur linéaire**.

```
Etape 1. Choix du comportement decrit par le modele Thermodynamique NRTL

tapez 1 pour une solution ideale et tapez 2 pour une solution non ideale
> 2

Etape 2. Choix du Solveur

tapez 1 pour le Solveur MRSL01 , tapez 2 pour le Solveur MRSL21
> 2

Etape 3. Debut du iterations

Iterations :    1 Norme :   74131.327847739245      IER : 0   SIMULATION EN COURS
Iterations :    2 Norme :   66720.613756641949      IER : 0   SIMULATION EN COURS
...
Iterations :  245 Norme :   8.3073038326122183E-007 IER : 0   SIMULATION EN COURS

Etape 4. Messages du programme

Temps de simulation du solveur en seconde :    7.1875000000000000
```

Les résultats sont écrits dans `ResultatsYsolMRSL21.csv` (profils) et `ResultatsResMRSL21.csv` (résidus finaux), dans le répertoire d'exécution. Toutes les autres données : nombre de plateaux, nombre de constituants, alimentations, soutirages, puissances, coefficients d'Antoine et d'interaction — sont lues dans les fichiers `Donnees*.txt`.

## 📊 Résultats de référence

Mélange Acétone – Benzène – Chloroforme, 30 plateaux, alimentation au plateau 14, convergence à ‖F‖₂ < 10⁻⁶ (exécutable Debug, `gfortran -g -Wall`) :

| Modèle | Solveur | Méthode | Itérations | Temps CPU |
|---|---|---|---|---|
| Idéal | MRSL01 | Matrice pleine | 239 | 11,0 s |
| Idéal | MRSL21 | Tridiagonal par blocs | 239 | **1,3 s** |
| Non idéal (NRTL) | MRSL01 | Matrice pleine | 247 | 17,8 s |
| Non idéal (NRTL) | MRSL21 | Tridiagonal par blocs | 245 | **7,2 s** |

Les deux solveurs convergent vers la même solution à environ 12 chiffres significatifs — ce qui vaut validation croisée. **MRSL21 est 2,5 à 8 fois plus rapide** et nettement plus économe en mémoire : il ne stocke que les blocs non nuls (3 × 9 × 9 × 30 réels) là où MRSL01 alloue la matrice pleine (270 × 270).

### Format des fichiers de sortie

Les deux fichiers CSV contiennent une ligne par plateau, séparateur `;`, avec **2·NC + 3 colonnes** de valeurs :

| Fichier | Colonnes |
|---|---|
| `ResultatsYsol*.csv` | Nº plateau ; **débit vapeur V** ; **compositions vapeur y (NC)** ; **température T** ; **compositions liquide x (NC)** ; **débit liquide L** |
| `ResultatsRes*.csv` | Nº plateau ; bilan matière total ; bilans matière partiels (NC) ; bilan enthalpique ; équilibres L/V (NC) ; sommation |

> ℹ️ Le message affiché en fin d'exécution annonce l'ordre « Température ; Débit Liquide ; Composition Liquide » ; l'ordre réellement écrit est celui du tableau ci-dessus, fixé par `Affecte_Indices.f90` (T, puis x, puis L).

## 📁 Structure du projet

```
distillation-multi-etagee/
├── README.md
├── doc/
│   ├── installation.pdf                              # Procedure d'installation Code::Blocks + MinGW
│   └── Rapport MOU2. Bangoura.Aboubacar.Sidiki.pdf   # Rapport du projet
└── DistllationMultiEtagee/
    ├── Distillation.cbp              # Projet Code::Blocks
    ├── Programme_Principale.f90      # Programme principal : choix mode/solveur, boucle de Newton
    │
    ├── mod_DIMENSIONS.f90            # Module : nombre de plateaux NP, de constituants NC
    ├── mod_THERMOS.f90               # Module : Antoine, Aij, Cij, Teb, Cp, deltaH, R, P, Tref
    ├── mod_OPERATOIRES.f90           # Module : alimentations, soutirages, puissances, D et W
    ├── mod_INDICES.f90               # Module : position des variables et equations dans le vecteur
    ├── mod_MESSAGES.f90              # Module : messages d'accueil, de fin et d'erreur
    │
    ├── Affecte_Dimensions.f90        # Lecture de NP et NC
    ├── Affecte_Thermos.f90           # Lecture Antoine, thermo et coefficients d'interaction
    ├── Affecte_Operatoires.f90       # Lecture alimentations, soutirages, puissances
    ├── Affecte_Indices.f90           # Calcul des index IV, IYvap, ITemp, IXliq, IL, BMT, BMP...
    │
    ├── Thermo_Psaturation.f90        # Pressions de vapeur saturante (Antoine etendue)
    ├── Thermo_Calculgama.f90         # Coefficients d'activite : mode 1 ideal / mode 2 NRTL
    ├── Thermo_CalculmKi.f90          # Coefficients de partage Ki = gamma·Psat/P
    ├── Thermo_CalculTbulle.f90       # Temperature de bulle par dichotomie
    ├── Thermo_CalculTrosee.f90       # Temperature de rosee par Newton-Raphson
    ├── Thermo_EnthalpieHliq.f90      # Enthalpie liquide (+ enthalpie d'exces en mode NRTL)
    ├── Thermo_EnthalpieHvap.f90      # Enthalpie vapeur
    │
    ├── Programme_Initialisation.f90  # Profil initial de T, x, y, L, V
    ├── Programme_Residus.f90         # Vecteur des residus MESH
    ├── Programme_Norme.f90           # Norme euclidienne du residu
    ├── Programme_Jacobien01.f90      # Jacobienne pleine (pour MRSL01)
    ├── Programme_Jacobien21.f90      # Jacobienne en blocs A, B, C (pour MRSL21)
    ├── MRSL01.f                      # Solveur : systeme lineaire plein, pivot partiel
    ├── MRSL21.f                      # Solveur : systeme tridiagonal par blocs
    │
    ├── DonneesDimensions.txt         # Fichiers de donnees d'entree
    ├── DonneesThermos.txt
    ├── DonneesCoefAntoine.txt
    ├── DonneesCoefInteraction.txt
    ├── DonneesOperatoires.txt
    ├── bin/Debug/                    # Executable genere (Distillation.exe)
    └── obj/Debug/                    # Fichiers objets et .mod
```

## ⚙️ Configuration

Les fichiers de données se trouvent dans `DistllationMultiEtagee/`. Le texte placé après `!` est un commentaire, ignoré à la lecture.

**`DonneesDimensions.txt`** dimensions de la colonne

```
30      ! Nombre de plateaux (condenseur = plateau 1, bouilleur = plateau NP)
3       ! Nombre de constituants
```

**`DonneesCoefAntoine.txt`** une ligne par constituant : A, B, C, D, E (Psat en Pa, T en K)

```
69.006 -5599.6 -7.0985 6.2237E-06 2    ! Acetone
83.107 -6486.2 -9.2194 6.9844E-06 2    ! Benzene
146.43 -7792.3 -20.614 0.024578  1     ! Chloroforme
```

**`DonneesThermos.txt`** une ligne par grandeur, une colonne par constituant

```
329.44 353.24 334.33         ! Temperature d'ebullition [K]
2.96E+04 3.08E+04 2.95E+04   ! Enthalpie de vaporisation [J/mol]
8.05E+01 9.60E+01 6.89E+01   ! Cp vapeur [J/mol/K]
1.34E+02 1.47E+02 1.17E+02   ! Cp liquide [J/mol/K]
```

**`DonneesCoefInteraction.txt`** NC lignes pour αᵢⱼ, puis NC lignes pour Cᵢⱼ (utilisé en mode 2 uniquement)

```
0       0.3007  0.3034        ! Aij : parametres de non-aleatoricite alpha
0.3007  0       0
0.3034  0       0
0            -808.935  -2691.471   ! Cij : energies d'interaction binaire
2384.591304  0         0
955.864088   0         0
```

**`DonneesOperatoires.txt`** conditions opératoires, une valeur (ou une liste) par ligne

```
1          ! Nombre de plateaux alimentes
14         ! Numero des plateaux alimentes
1          ! Debit des alimentations [mol/s]
332        ! Temperature des alimentations [K]
0.6        ! Fraction molaire du constituant 1 (Acetone) dans chaque alimentation
0.3        ! Fraction molaire du constituant 2 (Benzene)
0.1        ! Fraction molaire du constituant 3 (Chloroforme)
7.942E+04  ! Puissance du bouilleur [J/s]
6.27E+04   ! Puissance du condenseur [J/s]
0          ! Nombre de soutirages liquide   (0 si aucun)
0          ! Numero des plateaux soutires   (0 si aucun)
0          ! Debit des soutirages liquide   (0 si aucun)
0          ! Nombre de soutirages vapeur    (0 si aucun)
0          ! Numero des plateaux soutires   (0 si aucun)
0          ! Debit des soutirages vapeur    (0 si aucun)
```

Pour changer de colonne ou de mélange : mettre à jour `NP` et `NC`, ajuster les lignes correspondantes dans chaque fichier, puis relancer —> aucune recompilation n'est nécessaire. Les blocs de soutirage doivent rester présents même en l'absence de soutirage : conserver les zéros par défaut.

## 🔧 Développement

```
# Recompiler entierement le projet
Build > Rebuild            (Ctrl + F11)

# Nettoyer les fichiers objets
Build > Clean

# Passer en configuration optimisee
Build > Select target > Release
```

Le projet compile avec `-Wall` ; la cible *Debug* ajoute `-g`, la cible *Release* `-O2`.

Pour traquer une divergence ou un `NaN`, recompiler en ligne de commande avec des gardes supplémentaires :

```bash
gfortran -g -Wall -fcheck=all -finit-real=nan -ffpe-trap=invalid,zero,overflow ...
```

`-finit-real=nan` est particulièrement utile : il fait échouer immédiatement tout calcul qui s'appuierait sur une variable non initialisée, au lieu de laisser le résultat dépendre du contenu de la pile.

### État actuel

| Fonctionnalité | État |
|---|---|
| Mode 1 : solution idéale (γᵢ = 1) | ✅ opérationnel |
| Mode 2 : solution non idéale (NRTL) | ✅ opérationnel |
| Solveur MRSL01 — matrice pleine | ✅ opérationnel |
| Solveur MRSL21 — tridiagonal par blocs | ✅ opérationnel |
| Alimentations multiples | ✅ prévu par l'architecture, peu testé |
| Soutirages liquide et vapeur | ✅ prévu par l'architecture, peu testé |
| Condenseur total | 🔜 perspective (condenseur partiel uniquement) |
| Spécification par taux de reflux | 🔜 perspective (puissances imposées uniquement) |

### Limites connues

- **Division par une inconnue nulle.** La jacobienne numérique divise par `Ysol(j)·δ` (`Programme_Jacobien01.f90`, `Programme_Jacobien21.f90`). Si une inconnue passe par zéro, la fraction molaire d'un constituant très dilué, par exemple, le calcul produit un `NaN`. Une perturbation absolue, ou mixte (`max(|Y|, Y_ref)·δ`), lèverait cette fragilité.
- **Divergence silencieuse.** La boucle principale s'arrête sur `ERREUR > 1E-6` ; or `NaN > 1E-6` est **faux**. Une simulation qui diverge sort donc de la boucle comme si elle avait convergé, et écrit des fichiers CSV entièrement remplis de `NaN`. Un nombre maximal d'itérations et un test explicite de `NaN` restent à ajouter.
- **Amortissement fixe.** α = 0,1 est codé en dur dans `Programme_Principale.f90`. Il assure la robustesse mais impose plusieurs centaines d'itérations ; un amortissement adaptatif (α croissant à mesure que le résidu décroît) réduirait fortement le temps de calcul.
- **Le choix idéal / non idéal transite par un `common MODE`** partagé entre le programme principal, `CalculGama` et `EnthalpieHliq` ; une variable de module serait plus propre et plus sûre.
- **Intervalle de recherche de la température de bulle.** `CalculTbulle` procède par dichotomie sur [Tref, 1000 K] et s'arrête sur « Pas de solution dans l'intervalle de T donnee » si la solution n'y est pas encadrée.

## 📚 Contexte

Ce projet a été réalisé dans le cadre du **projet de modélisation des opérations unitaires (MOU2)** de **3ᵉ année**, parcours **CPAO** (Conception des Procédés Assistée par Ordinateur), à l'**ENSGTI** (École Nationale Supérieure en Génie des Technologies Industrielles), Université de Pau et des Pays de l'Adour, année universitaire **2022-2023**.

**Réalisé par** : Aboubacar Sidiki BANGOURA  
**Encadrant** : Pr Frédéric MARIAS

Le rapport complet (mise en équations, logigrammes, tests et analyse des résultats) est disponible dans [doc/Rapport MOU2. Bangoura.Aboubacar.Sidiki.pdf](doc/Rapport%20MOU2.%20Bangoura.Aboubacar.Sidiki.pdf).

Ce travail prolonge un projet antérieur de 1ʳᵉ année consacré au [calcul de la température de bulle](https://github.com/asbgra341/temperature-de-bulle), dont il reprend l'architecture modulaire et le principe des données externalisées, en y ajoutant le modèle NRTL, le système MESH complet et la résolution multi-plateaux.

## 📄 Licence

Projet académique, mis à disposition à des fins pédagogiques.

## 🙏 Remerciements

- Pr Frédéric MARIAS pour l'encadrement scientifique et la mise en équations du problème
- L'ENSGTI pour le cadre pédagogique
