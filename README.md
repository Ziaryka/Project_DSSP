# Project_DSSP

# DSSP
dssp.py est un programme dédié à la détermination du repliement des protéines (feuillet, hélice, boucle, tour...etc) grâce aux liaisons hydrogènes.
Ce script peut être appliqué à n'importe quel pdb (vous trouverez deux pdb à utiliser comme exemples).
Ce script peut également être trouvé sur github à l'adresse: https://github.com/YFeriel/dssp.git

Pour trouver la liaison hydrogène entre les résidus constituant la protéine:

- Il faudra calculer la distance entre les atomes N,O,H,C de chaque résidu.
- Calculer l'énergie d'interaction.

# Installation :

Avant d’exécuter dssp.py, il est impératif de suivre les lignes de commande présentent ci dessous : 

# Activation de l’environnement Conda : (lignes séparées)

conda config --add channels conda-forge

conda update --yes conda


# Installation des packages nécessaires : 
math
conda install -c conda-forge -c schrodinger pymol-bundle

# Les fichiers d’entrée (input) : 

L'input doit être au format .pdb.

argument 2 = dssp.py
argument 3 = fichier pdb
argument 4 = présence de liaisons H ou non dans le pdb

# Vous trouverez dans ce script :

- Une interface utilisateur et de vérification permettant de déterminer si tous les arguments ont été entrés et répondent aux conditions. Le but de cette fonction est de s’assurer que le fichier d’entrée a bien été importé, n'est pas vide et présente la bonne extension.

- Une fonction "calcul_dist" qui permet de calculer la distance entre les atomes O,N,C,H des résidus.

- une fonction calcul_energy" qui permet à partir des distances enregistrées dans la précédente fonction, les énergies de liaison.

- Une fonction main() qui permet de démarrer et d’exécuter le programme en tant que programme principal. Elle renvoie les sorties de chaque fonction définie précédemment. Et pour ce faire, il faudra rentrer les lignes suivantes : 
import dssp
dssp.main()
		
# Comprendre les lignes de commande
- ligne 70 à 90 : ici, le code extrait chaque atome (N,O,H,C) du fichier pdb et les conserve dans des listes. La liste "elements" contient le nom et la position de chaque résidu. Cette dernière nous sera utile ultérieurement.
- ligne 92 à 110 : ces lignes de commande servent à enregistrer chaque distance calculée dans des fichiers.
- ligne 112 à 136 : ces lignes permettent de conserver dans des listes, uniquement la valeur des distances afin de les utiliser dans le calcul des énergies.
- ligne 147 à 186 : ici, on commence d'abord par déterminer un cutoff (<-0.5kcal/mol) pour définir une liaison hydrogène établie entre deux résidus. Ensuite, nous créons une matrice appelée "final" contenant des indices (0 pour absence de liaison hydrogène, et 1 pour liaison hydrogène).

# Détermination des repliements : 

- ligne 188 à 237 : ces lignes permettent de détecter les hélices (une hélice est une succession de deux tours minimum) et brins.
- ligne 239 à 271: ces lignes génèrent les output en stockant la détermination des structures dans un fichier, avec la première colonne indiquant le code du résidu, et la deuxième colonne indiquant si ce dernier est impliqué ou non dans une hélice ou un brin.

# NB 1 : 
Nous avons volontairement laissé le calcul entre deux paires de résidus se répeter (exemple: Distance entre SER1-GLY2 et GLY2-SER1) et ce, dans le but de nous faciliter l'implémentation de la matrice et plus précisémment la diagonale.
# NB 2 :
Nous avons également choisi de remplacer les atomes H par les atomes CA, afin d'implémenter notre algorithme DSSP, car dans les fichiers pdb, les structures résolues par cristallographie ne contiennent pas de H. Les résultats dans les output ainsi obtenus peuvent être biaisés et faussés.
Néanmoins, nous avons essayé d'implémenter un code pour remédier à cela. Ainsi, pour lancer le script avec cette option :
- Si vous n'avez pas ajouté en amont les atomes H, veuillez mettre comme 4ème argument no_h .
- Sinon, veuillez mettre yes_h .
