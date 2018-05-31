Pocketanneal: Creating protein-ligand interfaces using AutoDock and Simulated annealing algorithms.
==============

Copyright (C) 2015 Deepesh Nagarajan (1337deepesh@gmail.com)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, version 3 of the License.

This program is distributed in the hope that it will be useful, but 
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License <http://www.gnu.org/licenses/>
for more details.

This README file will guide you through an example 'Pocketanneal' usage:

1) Enter the Pocketanneal directory through terminal
    **cd Pocketanneal**

2) Run 'Format.sh' to add hydrogens & clean your input PROTEIN.pdb file:
    **bash Format.sh -i example/protein_1dhi_A.pdb -o ~/Desktop/Pocketanneal_OUTPUT_01 -d database/**

3) After pre-processing, run Pocketanneal with these arguments:
    **python Pocketanneal.py -P ~/Desktop/Pocketanneal_OUTPUT_01/output.pdb -p example/pocket_1dhi_A.pdb -d database -n 1000 -a 2 -o ~/Desktop/Pocketanneal_OUTPUT_02**

A more detailed guide explaining every step can be found in:
    **gedit Pocketanneal/example/example.sh**

DEPENDENCIES:
--------------
    Install autodocktools:
    **sudo apt-get install autodocktools**
    On Ubuntu 16 onwards, you have to downgrade numpy:
    **sudo pip install numpy==1.8**
    
    This script runs on Python 2.7.12 (Ubuntu 16.04). Don't run on Python 3.
    I can't make any guesses about how Pocketanneal will run on future Ubuntu versions.

BUGS:
--------------
Ignore this output on terminal (generated by autodock, makes no difference to the score output):
**swig/python detected a memory leak of type 'BHtree *', no destructor found.**

Python image libraries are in a state of flux. Simply comment out line 20 and line 995 in Pocketanneal.py if the Image library is giving you trouble.
**import Image, ImageDraw, ImageFont** (Pocketanneal.py, line 26)
**grapher(XY_data, str(working_directory)+"/score_log.png", str(IP_database))** (Pocketanneal.py, line 995)


Contact
--------------
If you have any queries/suggestions, contact me:
Deepesh Nagarajan: 1337deepesh@gmail.com

I will eventually leave my current position, so I can't guarantee indefinite support.









