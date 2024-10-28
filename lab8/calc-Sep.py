import sys
import pandas as pd

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

    # Open PDB File
    pdbFileChars = open(f"pdbFiles/{structID}.cif", "r")
    pdbFileChars = pdbFileChars.read()

    # Check if ligand in pdb file
    ligandStr = f"{ligandID}{4%len(ligandID) * ' '}non-polymer"
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
    pdbFile = pdbFile.readline()
    df = pd.DataFrame(columns=["type", "atomNum", "elementId", "ligandId", "x", "y", "z"])
    for newline in pdbFile:
        line = newline.split(" ")
        new_row = {"type",
                   "atomNum",
                   "elementId",
                   "ligandId",
                   "x",
                   "y",
                   "z"}
        if ("HETATM" in newline):
            breakpoint()
            new_row["type"] = "HETATM"
            new_row["ligandId"] = line[21:24]
        elif (line[:4]=="ATOM"):
            new_row["type"]="ATOM"
        else:
            continue
        new_row["atomNum"] = line[7:11]
        new_row["elementId"] = line[12]
        new_row["x"] = line[34:41]
        new_row["y"] = line[42:48]
        new_row["z"] = line[49:56]











if __name__=="__main__":
    main()
