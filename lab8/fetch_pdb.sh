#!/bin/bash

# declare arrays of PDBs
pdbIDs=('5con' '5coo' '5cop' '4hla' '5bry' '5bs4' '4hdp' '5dgw', '3cyw')
target="5dgu"
numPDBs=${#pdbIDs[@]}
echo $target
echo ${pdbIDs[@]}
echo $numPDBs

if [-d "pdbFiles"]; then
    exit

else
    mkdir "pdbFiles"

fi

# fetch all of the drug candidate structures
for pdbID in "${pdbIDs[@]}"
do
    echo "fetch $pdbID, async=0" >> fetchCmnds.pml
done

# save pdf files
for pdbID in "${pdbIDs[@]}"
do
    echo "save ${pdbID}.pdb, $pdbID" >> fetchCmnds.pml
done

# invoke pymol
pymol fetchCmnds.pml -qc

mv *.cif ./pdbFiles
mv *.pdb ./pdbFiles
