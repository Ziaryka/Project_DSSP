#SCRIPT POUR ANALYSE DSSP

#importation des librairies nécessaires
import sys, os, math


###INTERFACE UTILISATEUR

#Importation des fichiers

try:
    pdb_file = sys.argv[1]
    h = sys.argv[2]

#Verification des fichiers:

#Si biens specifies:
except:
    sys.exit("Erreur: le fichier choisi pour l'analyse DSSP ou la présence ou non d'H dans le fichier pdb n'ont pas été spécifié.\
\nVeuillez relancer la commande en précisant en 3ème argument le fichier PDB et en 4ème argument: no_h (si le pdb ne contient pas d'H) ou yes_h (si le pdb en contient).")

#Si bon format:
if not pdb_file.endswith("pdb"):
    sys.exit("Erreur: Mauvais format. Le fichier choisi doit être au format \"pdb\".")

#Si existe et bien rempli:
if not os.path.exists(pdb_file) or os.stat(pdb_file).st_size == 0:
    sys.exit("Le fichier choisi est vide ou n'existe pas.")

#PREPARATION FICHIER PDB
#Ajout des H
if h =="no_h":
    import pymol2
    with pymol2.PyMOL() as pymol:
        pymol.cmd.load(pdb_file, 'prot')
        pymol.cmd.h_add()
        pymol.cmd.save('prot_H.pdb')
    pdb_file = "prot_H.pdb"

elif h =="yes_h":
    pass

else:
    sys.exit("Erreur: La présence d'hydrogène ou non dans le PDB a été spécifié de manière incorrecte.\
\nVeuillez relancer la commande en précisant en 4ème argument: no_h (si le pdb ne contient pas d'H) ou yes_h (si le pdb en contient).")

###CALCULS

def calcul_dist(line1, line2):  #fonction pour calculer la distance
    cor1 = line1.split()
    cor2 = line2.split()
    x1, y1, z1 = float(cor1[6]), float(cor1[7]), float(cor1[8])
    x2, y2, z2 = float(cor2[6]), float(cor2[7]), float(cor2[8])
    distance_on = math.sqrt(math.pow(x1-x2, 2) + math.pow(y1-y2, 2) + math.pow(z1-z2, 2))
    return cor1[3]+" "+ cor1[5]+"-"+cor2[3]+" "+cor2[5]+ " " " : {0:.3f}\n".format(distance_on)

#Fonction energie :
def calcul_energy(distance_on,distance_oh,distance_ch,distance_cn):
    e = 1.06217662*(10^-19)
    q1= 0.42*e
    q2= 0.2*e
    f= 332
    energy=[]
    for i in range(len(distance_on)):
        energy.append((q1*q2)*((1/distance_on[i])+(1/distance_ch[i])-(1/distance_oh[i])-(1/distance_cn[i]))*f)
    return energy

def main():
    print("dssp__name__==%s" % __name__)
    o_atome=[]
    n_atome=[]
    alpha=[]
    c_atome=[]

    elements=[]
    elt = []
    with open("chaine1.pdb") as chain:
        for line in chain:

            if line[13:15].startswith('O') and line[13:15].endswith(' '):
                o_atome.append(line)
                elt=line.split()
                chaine=elt[3]+elt[5]
                elements.append(chaine)
            if line[13:15].startswith('N') and line[13:15].endswith(' '):
                n_atome.append(line)
            if line[13:15].startswith('C') and line[13:15].endswith(' '):
                c_atome.append(line)
            if line[13:15]=='CA':
                alpha.append(line)

    with open("Distance_O_N","w") as calculator:
        for i in range(len(o_atome)):
            for j in range(len(n_atome)):
                calculator.write(str(calcul_dist(o_atome[i],n_atome[j])))

    with open("Distance_O_H","w") as calculator:
        for i in range(len(o_atome)):
            for j in range(len(alpha)):
                calculator.write(str(calcul_dist(o_atome[i],alpha[j])))

    with open("Distance_C_N","w") as calculator:
        for i in range(len(c_atome)):
            for j in range(len(n_atome)):
                calculator.write(str(calcul_dist(c_atome[i],n_atome[j])))

    with open("Distance_C_H","w") as calculator:
        for i in range(len(c_atome)):
            for j in range(len(n_atome)):
                calculator.write(str(calcul_dist(c_atome[i],alpha[j])))

    distance_on = []
    distance_cn = []
    distance_oh = []
    distance_ch = []

    temp = []

    with open("Distance_O_N","r") as f:
        for lines in f:
            line = lines.split()
            distance_on.append(float(line[-1]))
    with open("Distance_C_N","r") as f:
        for lines in f:
            line = lines.split()
            distance_cn.append(float(line[-1]))
    with open("Distance_O_H","r") as f:
        for lines in f:
            line = lines.split()
            distance_oh.append(float(line[-1]))
    with open("Distance_C_H","r") as f:
        for lines in f:
            line = lines.split()
            distance_ch.append(float(line[-1]))

    energy = calcul_energy(distance_on,distance_cn,distance_ch,distance_oh)

    i=0
    with open("Distance_O_N","r") as f:
        with open("Energy","w") as e :
            for lines in f:
                line = lines.split()
                e.write(line[0]+line[1]+line[2]+line[3]+" "+str(energy[i])+"\n")
                i+=1


    liaison_h = []
    pas_liaison = []
    results=[]
    with open("Energy","r") as f:
        for lines in f :
            line = lines.split()
            line[-1] = float(line[-1])
            if float(line[-1]) < -0.5 :
                liaison_h.append(line)
                results.append(1)
            else:
                pas_liaison.append(line)
                results.append(0)
    cpt=1
    header = [""]
    for i in range(len(elements)):
        header.append(elements[i])
        cpt+=1

    final = []
    final.append(header)

    temp=[]

    j=0
    k=1
    var = 0
    temp.append(header[1])
    while j < len(results):
        temp.append(results[j])
        var+=1
        if var == len(o_atome):
            k += 1
            final.append(temp)
            temp=[]
            if k == len(o_atome)+1:
                break
            temp.append(header[k])
            var=0
        j += 1

#Assignement des structures secondaires

    #Helices
    resid_helix=[]
    x=0
    for i in range(len(final)) :
        for j in range(len(final[i])-3):
            if final[i][j]==1:
                x+=1
            else:
                x=0
            if(x==3 and final[i][j+1]==0):
                #print(final[0][j-2],final[i][0],'+',final[0][j-1],final[i][0],'+',final[0][j],final[i][0])
                resid_helix.extend([final[0][j-2],final[i][0],final[0][j-1],final[i][0],final[0][j],final[i][0]])

            if(x==3  and final[i][j+1] ==1) :
                #print(final[0][j - 2], final[i][0], '+', final[0][j - 1], final[i][0], '+', final[0][j], final[i][0],'+',final[0][j+1],final[i][0])
                resid_helix.extend([final[0][j - 2],final[i][0],final[0][j - 1],final[i][0],final[0][j], final[i][0],final[0][j+1],final[i][0]])

            if(x==3  and final[i][j+1] ==1 and final[i][j+2]==1 and final[i][j+3]==0 ) :
                #print(final[0][j - 2], final[i][0], '+', final[0][j - 1], final[i][0], '+', final[0][j], final[i][0],'+',final[0][j+1],final[i][0],'+',final[0][j+2],final[i][0])
                resid_helix.extend([final[0][j - 2], final[i][0],final[0][j - 1], final[i][0], final[0][j], final[i][0],final[0][j+1],final[i][0],final[0][j+2],final[i][0]])

    # Bridges
    resid_bridges = []
    x = 0
    for i in range(len(final) - 1):
        for j in range(len(final[i]) - 2):
            if final[i][j] == 1:
                x += 1
            else:
                x = 0
            # parallel bridges
            if x == 2 and final[i - 1][j] == 1 and final[j][i + 1] == 1:
                # print(final[0][j - 1], final[i][0],'+',final[0][j], final[i][0])
                resid_bridges.extend([final[0][j - 1], final[i][0], final[0][j], final[i][0]])
            if x == 2 and final[j - 1][i] == 1 and final[i][j + 1] == 1:
                # print(final[0][j - 1], final[i][0],'+',final[0][j], final[i][0])
                resid_bridges.extend([final[0][j - 1], final[i][0], final[0][j], final[i][0]])

            # antiparallel bridges
            if x == 2 and final[i][j] == 1 and final[j][i] == 1:
                # print(final[0][j - 1], final[i][0],'+',final[0][j], final[i][0])
                resid_bridges.extend([final[0][j - 1], final[i][0], final[0][j], final[i][0]])
            if x == 2 and final[i - 1][j + 1] == 1 and final[j - 1][i + 1] == 1:
                # print(final[0][j - 1], final[i][0],'+',final[0][j], final[i][0])
                resid_bridges.extend([final[0][j - 1], final[i][0], final[0][j], final[i][0]])

    resid_helix = list(dict.fromkeys(resid_helix))
    resid_bridges = list(dict.fromkeys(resid_helix))

    # Output helices
    allresid = []
    for i in range(len(final)):
        allresid.append(final[i][0])

    allresid_set = set(allresid)
    intersection = allresid_set.intersection(resid_helix)
    intersection_list = list(intersection)

    helix = ["not involved in an helix"] * (len(allresid))
    for resid1 in allresid:
        for i in range(len(intersection_list)):
            if resid1 == intersection_list[i]:
                helix[i] = "H"

    with open("DSSP_results_helix", "w") as output:
        for i in range(len(allresid)):
            output.write(str(allresid[i]) + "-" + str(helix[i]) + "\n")

    # Output bridges
    allresid_set = set(allresid)
    intersection2 = allresid_set.intersection(resid_bridges)
    intersection2_list = list(intersection2)

    bridge = ["not involved in a bridge"] * (len(allresid))
    for resid1 in allresid:
        for i in range(len(intersection2_list)):
            if resid1 == intersection2_list[i]:
                bridge[i] = "B"

    with open("DSSP_results_bridge", "w") as output:
        for i in range(len(allresid)):
            output.write(str(allresid[i]) + "-" + str(bridge[i]) + "\n")
if __name__ == "__main__" :
    main()
