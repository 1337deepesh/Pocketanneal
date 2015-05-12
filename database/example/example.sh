#This script will guide you through an example 'Pocketanneal' usage:

#Inputs:
cat protein_1dhi_A.pdb
#This file contains all ATOMS of the selected protein. HETATM inclusion is not necessary.
cat pocket_1dhi_A.pdb
#This file must contain only those residues constituting the ligand-binding pocket.
#Ligand HETATM coordinates are needed here. They must be placed after all ATOM coordinates.
#Note that ligand hydrogens may or may not be included, depending on user choice.
#If ligand hydrogens are not included, no attempt will be made to include them.

#Step-1: format the protein:
#formatting is required to check/repair syntax irregularities and structural anomalies.
#Furthermore, existing hydrogen atoms (if any) will be stripped and replaced.
Pocketanneal_directory=<insert location here>
Pocketanneal_database=<insert location here>
bash $Pocketanneal_directory/Format.sh -i protein_1dhi_A.pdb -o Format_output -d $Pocketanneal_database

#Note that the pocket does not need to be formatted.

#The formatted protein will be available in this file:
cat Format_output/output.pdb

#Once the protein is formatted, Autanneal.py will redesign the protein-ligand interface:
python $Pocketanneal_directory/Pocketanneal.py -P format_OP/output.pdb -p pocket_1dhi_A.pdb -d $Pocketanneal_database -n 1000 -a 2 -o output

#A typical run takes 10-20 minutes.

#for usage explanations, run the command without flags:
python $Pocketanneal_directory/Pocketanneal.py

#The output folder will contain the file 'BEST_protein.pdb'. This file contains the redesigned binding interface. The file 'BEST_score.csv' contains the corresponding AutoDock score of the best designed interface.
