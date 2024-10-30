import sys
import stat
import math
import pandas as pd
import os
import subprocess

RADII = {
        "O": 0.48,
        "C": 0.67,
        "N": 0.56,
        "P": 0.98
        }

def main():
    # cl args:
    # 1) 4-alphanumeric id for a structure in the PDB
    # 2) a Ligand ID
    # 3) one of O, C, N, and P
    structID = sys.argv[1]
    ligandID = sys.argv[2]
    element = sys.argv[3]

    # Check args
    if (len(structID)!=4 and structID not in ['5con', '5coo', '5cop', '4hla', '5bry', '5bs4', '4hdp', '5dgw', '3cyw']):
        print("Invalid structID")
        return

    # Run bash script if pdbFiles folder does not exist
    if os.path.isdir('pdbFiles') is False:
        os.chmod('./fetch_pdb.sh', os.stat('./fetch_pdb.sh').st_mode | stat.S_IEXEC)
        subprocess.call("./fetch_pdb.sh")

    # Open PDB File
    pdbFileChars = open(f"pdbFiles/{structID}.cif", "r")
    pdbFileChars = pdbFileChars.read()

    # Check if ligand in pdb file
    ligandStr = f"{ligandID}{(4-len(ligandID)) * ' '}non-polymer"
    if (ligandStr not in pdbFileChars):
        print("Ligand not found in pdbFile")
        return

    # Check element
    if (element not in 'ONCP'):
        print("Invalid element")
        return


    '''
    Calculate the 10 distinct atoms in the ligand (HETATM) that are closest to the atoms (ATOM) in the PDB file
    and print:
        a) atom ID and name in the ligand
        b) atom ID and name in protein ordered by sep distance
        c) their sep distance
        d) the O,C,N, or P atom specified as the 3rd cla fits in each of the sep distances output by the prog
    '''
    pdbFile = open(f"pdbFiles/{structID}.cif", "r")
    pdbFile = pdbFile.readlines()
    df = pd.DataFrame(columns=["type", "atomNum", "elementId", "ligandId", "x", "y", "z"])

    hetatms = []
    atoms = []

    for newline in pdbFile:
        line = newline.split()
        new_row = {"type":None,
                   "atomNum":None,
                   "elementId":None,
                   "ligandId":None,
                   "x":None,
                   "y":None,
                   "z":None}

        if ("HETATM" in newline):
            new_row["type"] = "HETATM"
            new_row["ligandId"] = line[5]
        elif (line[0]=="ATOM"):
            new_row["type"]="ATOM"
        else:
            continue

        new_row["atomNum"] = line[1]
        new_row["elementId"] = line[2][0]
        new_row["x"] = float(line[10])
        new_row["y"] = float(line[11])
        new_row["z"] = float(line[12])

        # Store HETATMs with specified ligandID
        if (new_row['ligandId'] == ligandID):
            hetatms.append(new_row)

        # Store ATOMs
        if (new_row['type'] == 'ATOM'):
            atoms.append(new_row)


    # Calculate distances between all HETATMs and ATOMS
    distances = []
    for hetatm in hetatms:
        for atom in atoms:
            distance = round(math.sqrt(((hetatm['y'] - atom['y'])**2 + (hetatm['y'] - atom['y'])**2 + (hetatm['z'] - atom['z'])**2)), 2)
            fitness = f"{element} doesn't fit" if RADII[element] > distance else f"{element} fits"
            string = f"HETATM {hetatm['atomNum']} {hetatm['elementId']} {hetatm['ligandId']}, ATOM {atom['atomNum']} {atom['elementId']} {structID} separated by {distance}A, {fitness}"
            distances.append((distance, string))

    # Take top 10
    distances.sort(key=lambda x: x[0])
    for dist in distances[:10]:
        print(dist[1])


if __name__=="__main__":
    main()
