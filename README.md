# Stage M2 – Projet GazFrom
Ce dépôt contient les scripts utilisés pour l’analyse métatranscriptomique de communautés microbiennes de fromages modèles à pâte pressée cuite dans le cadre du projet GazFrom.
Les analyses incluent :
- le contrôle qualité des lectures RNA-seq avec FastQC et MultiQC ;
- le trimming des lectures avec Cutadapt ;
- la filtration des ARN ribosomiques avec SortMeRNA ;
- l’alignement des lectures sur les génomes bactériens de référence avec Bowtie2 ;
- la quantification des lectures associées aux gènes avec featureCounts ;
- l’analyse différentielle d’expression avec DESeq2 ;
- la visualisation des profils transcriptomiques par ACP ;
- la représentation des gènes différentiellement exprimés par volcano plots ;
- l’estimation de l’abondance transcriptionnelle relative des espèces du consortium ;
- l’annotation fonctionnelle et l’analyse d’enrichissement KEGG.

Les données brutes de séquençage ne sont pas incluses dans ce dépôt pour des raisons de volume.

## Organisation générale du pipeline
Le pipeline suit les principales étapes suivantes :
1. contrôle qualité des données brutes ;
2. trimming et nettoyage des lectures ;
3. contrôle qualité post-trimming ;
4. filtration des ARN ribosomiques ;
5. indexation des génomes de référence ;
6. alignement des lectures non ribosomiques ;
7. quantification des lectures par gène ;
8. analyses statistiques et graphiques sous R.

## Structure du dépôt
```text
GazFrom/
└── scripts/
    ├── 01_pretraitement/
    ├── 02_mapping_quantification/
    ├── 03_expression_differentielle/
    └── 04_annotation_fonctionnelle/
```
## Organisation des scripts
Les scripts sont organisés selon les principales étapes du pipeline bioinformatique :

- 01_pretraitement/ : contrôle qualité, trimming et filtration des ARN ribosomiques ;
- 02_mapping_quantification/ : indexation des génomes, alignement et quantification ;
- 03_expression_differentielle/ : analyses statistiques et visualisations transcriptomiques ;
- 04_annotation_fonctionnelle/ : annotation eggNOG et enrichissement fonctionnel KEGG.

## Environnement logiciel
Les analyses ont été réalisées sous Linux sur l’infrastructure Migale (INRAE).
Les principaux outils utilisés incluent :
- FastQC ;
- MultiQC ;
- Cutadapt ;
- SortMeRNA ;
- Bowtie2 ;
- samtools ;
- featureCounts ;
- DESeq2 ;
- eggNOG-mapper ;
- clusterProfiler.

Les dépendances logicielles étaient gérées avec Conda.

## Principales sorties du pipeline
Le pipeline génère notamment :
- des rapports de contrôle qualité FastQC/MultiQC ;
- des fichiers FASTQ nettoyés ;
- des fichiers BAM alignés ;
- des matrices de comptages featureCounts ;
- des résultats d’expression différentielle DESeq2 ;
- des ACP transcriptomiques ;
- des volcano plots ;
- des tableaux d’annotation fonctionnelle ;
- des analyses d’enrichissement KEGG.

## Exécution des scripts sur Migale
Les analyses ont été exécutées sur l’infrastructure Migale (INRAE) à l’aide de scripts de soumission SGE (qsub) adaptés à chaque étape du pipeline.

Pour des raisons de lisibilité, seuls les principaux scripts d’analyse bioinformatique sont inclus dans ce dépôt. Les scripts de lancement (run_*.sh) utilisés pour soumettre les jobs sur le cluster suivent tous la même structure générale : définition des chemins, vérification des entrées, chargement des environnements Conda et exécution des scripts Bash ou R correspondants.

## Important
**Avant d’exécuter les scripts, vérifiez et adaptez les chemins d’accès aux répertoires de travail selon votre installation locale ou votre environnement de calcul.
En particulier, les variables correspondant aux répertoires d’entrée, de sortie, aux annotations et aux environnements logiciels doivent être modifiées si nécessaire.**

## Utilisation de l’intelligence artificielle

Certaines reformulations rédactionnelles, aides à la structuration de la documentation ainsi que l’harmonisation des commentaires des scripts ont été réalisées avec l’assistance d’un modèle d’intelligence artificielle générative (OpenAI ChatGPT).

L’auteur a systématiquement relu, vérifié, corrigé et validé l’ensemble du contenu scientifique, des analyses bioinformatiques, des scripts, des résultats et des références bibliographiques, et assume l’entière responsabilité du contenu final du dépôt GitHub et du mémoire de stage.
