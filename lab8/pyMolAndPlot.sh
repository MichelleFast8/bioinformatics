#!/bin/bash
echo "hello world"

# declare arrays of PDBs
pdbIDs=('5con' '5coo' '5cop' '4hla' '5bry' '5bs4' '4hdp' '5dgw')
target="5dgu"
numPDBs=${#pdbIDs[@]}
echo $target
echo ${pdbIDs[@]}
echo $numPDBs

# fetch the target pdbID
echo "fetch 5dgu, async=0" > pymolCmnds.pml

# fetch all of the drug candidate structures
for pdbID in "${pdbIDs[@]}"
do
    echo "fetch $pdbID, async=0" >> pymolCmnds.pml
done

# stucturally align each candidate structure against the target
for pdbID in "${pdbIDs[@]}"
do
    echo "align $pdbID & res 1-40+60-90 & name n+ca+c+o, $target & res 1-40+60-90 & name n+ca+c+o" >> pymolCmnds.pml
    echo "save ${pdbID}_aligned.pdb, $pdbID" >> pymolCmnds.pml
done

# invoke pymol
pymol pymolCmnds.pml -qc

# retrieve the CA atom IDs of the residues 25 and 51 of chain A of all structures
for i in `seq 0 7`;
do
    aa25x[$i]=`grep "A  25" ${pdbIDs[$i]}_aligned.pdb | grep "ATOM" | grep "CA" | head -1 | awk '{print $7}'`
    aa25y[$i]=`grep "A  25" ${pdbIDs[$i]}_aligned.pdb | grep "ATOM" | grep "CA" | head -1 | awk '{print $8}'`
    aa25z[$i]=`grep "A  25" ${pdbIDs[$i]}_aligned.pdb | grep "ATOM" | grep "CA" | head -1 | awk '{print $9}'`

    aa51x[$i]=`grep "A  51" ${pdbIDs[$i]}_aligned.pdb | grep "ATOM" | grep "CA" | head -1 | awk '{print $7}'`
    aa51y[$i]=`grep "A  51" ${pdbIDs[$i]}_aligned.pdb | grep "ATOM" | grep "CA" | head -1 | awk '{print $8}'`
    aa51z[$i]=`grep "A  51" ${pdbIDs[$i]}_aligned.pdb | grep "ATOM" | grep "CA" | head -1 | awk '{print $9}'`
done

# perform flap/active site distance calculations
for i in `seq 0 7`;
do
    d=`echo "(${aa51x[$i]} - ${aa25x[$i]})^2+(${aa51y[$i]} - ${aa25y[$i]})^2+(${aa51z[$i]} - ${aa25z[$i]})^2" | bc`
    distance[$i]=`echo "sqrt($d)" | bc`
done

echo ${distance[@]}

# make gnuplot data file
rm gnuplotData.txt
for i in `seq 0 7`;
do
    echo "$((i+1)) ${distance[$i]}" >> gnuplotData.txt
done

# gnuplot commands for terminal, x and y ranges, labels, title, etc.
echo "set term png" > gnuplot.gpl
echo "set output 'flapCatSitePlot.png'" >> gnuplot.gpl
echo "set xrange[0:9]" >> gnuplot.gpl

# to-do : Add an X label, Y label, and Title for the plot
echo "set xlabel 'Structure'" >> gnuplot.gpl
echo "set ylabel 'Separation Distance (Angstrom)'" >> gnuplot.gpl
echo "set title 'Flap-Catalytic Site Separation Distance, 8 HIV-1 Protease Variants" >> gnuplot.gpl

# make xtic labels and generate the plot, and invoke gnuplot
xticsString="("
for i in `seq 0 7`;
do
    xticsString+="'${pdbIDs[$i]}' $((i+1)),"
done
xticsString+="'${pdbIDs[7]}' $((i+1)))"
echo "set xtics $xticsString" >> gnuplot.gpl
echo "plot 'gnuplotData.txt' using 1:2 with lines title 'sep distance'" >> gnuplot.gpl
gnuplot gnuplot.gpl
