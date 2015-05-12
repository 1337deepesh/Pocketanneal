#!/usr/bin/python
#Script started on: Sunday 13 April 2014 10:25:55 PM IST  
#LDI modification started on: Tue 24 Mar 2015 04:16:47 PM IST 
import os
import csv
import sys
import copy
import math
import numpy
import string
import random
import shutil
import getopt
import subprocess

import AutoDockTools
import MolKit.protein
import MolKit.molecule
from MolKit import Read
from PyAutoDock.InternalEnergy import InternalEnergy
from PyAutoDock.AutoDockScorer import AutoDock41Scorer
from PyAutoDock.MolecularSystem import MolecularSystem
from AutoDockTools.MoleculePreparation import AD4LigandPreparation
from AutoDockTools.MoleculePreparation import AD4ReceptorPreparation

import Image, ImageDraw, ImageFont

########################################
#       FUNCTION: GRAPH-PLOTTING       #
########################################

#create graph of score_log (recycle code from dinoplot):
def grapher(XY_data, OP_image, IP_database):
    
    #create a 'diamond' that will serve the standard graph-point:
    diamond_blue = [[255, 255, 255, 255, 255, 255],
                    [255, 255, 255, 255, 255, 255],
                    [255, 255, 255, 255, 255, 255],
                    [255, 255, 255, 255, 255, 255],
                    [255, 255, 255, 255, 255, 255],
                    [255, 255, 255, 255, 255, 255]]
    diamond_null = [[255, 255,   0,   0, 255, 255],
                    [255,   0,   0,   0,   0, 255],
                    [  0,   0,   0,   0,   0,   0],
                    [  0,   0,   0,   0,   0,   0],
                    [255,   0,   0,   0,   0, 255],
                    [255, 255,   0,   0, 255, 255]]
    
    diamond = [diamond_null, diamond_null, diamond_blue]
    
    W = 1000
    H = 1000
    
    img = Image.new("RGB", (W, H), "white")
    draw = ImageDraw.Draw(img)
    
    #create output image margin:
    for i in range(99,901):
        color = (0,0,0)
        draw.point((i,99), fill=color)
        draw.point((i,100), fill=color)
        draw.point((i,901), fill=color)
        draw.point((i,902), fill=color)
        draw.point((99,i), fill=color)
        draw.point((98,i), fill=color)
        draw.point((901,i), fill=color)
        draw.point((902,i), fill=color)

    #calculate axes range:
    X_min = float("inf")
    X_max = float("inf")*-1
    Y_min = float("inf")
    Y_max = float("inf")*-1

    for i in range(0, len(XY_data)):
        if XY_data[i][0] > X_max:
            X_max = float(XY_data[i][0])
        if XY_data[i][0] < X_min:
            X_min = float(XY_data[i][0])
        if XY_data[i][1] > Y_max:
            Y_max = float(XY_data[i][1])
        if XY_data[i][1] < Y_min:
            Y_min = float(XY_data[i][1])
    
    #calculate all points as fractions of min/max, convert to pixel location:
    #also plot on graph:
    for i in range(0, len(XY_data)):
        Y_point = (((XY_data[i][0]-X_min)/(X_max-X_min)*800)+100)
        X_point = 1000-(((XY_data[i][1]-Y_min)/(Y_max-Y_min)*800)+100)
        for j in range(0,len(diamond[0])):
            for k in range(0,len(diamond[0][0])):
                if diamond[0][j][k] == 0:
                    K = X_point+j-math.ceil(len(diamond[0])/2)
                    J = Y_point+k-math.ceil(len(diamond[0][0])/2)
                    color = (diamond[0][j][k], diamond[1][j][k], diamond[2][j][k])
                    draw.point((J,K), fill=color)

    #draw axes labels, ranges:
    font = ImageFont.truetype(IP_database+"/formatting/arial.ttf", 30)
    draw.text((445,950), "iteration", fill=(0,0,0), font=font)
    draw.text((45,420), "s", fill=(0,0,0), font=font)
    draw.text((45,450), "c", fill=(0,0,0), font=font)
    draw.text((45,480), "o", fill=(0,0,0), font=font)
    draw.text((45,510), "r", fill=(0,0,0), font=font)
    draw.text((45,540), "e", fill=(0,0,0), font=font)
    draw.text((310,50), "Simulated annealing (score log)", fill=(0,0,0), font=font)
    draw.text((100,930), str(X_min), fill=(128,128,128), font=font)
    draw.text((900,930), str(X_max), fill=(128,128,128), font=font)
    draw.text((10, 890), '% 4.2f' % Y_min, fill=(128,128,128), font=font)
    draw.text((10,90), '% 4.2f' % Y_max, fill=(128,128,128), font=font)

    #output image:
    img.save(OP_image, "png")

########################################
#       FUNCTION: CLASH-DETECTOR       #
########################################

#check for steric clashes between any 2 residues:
#PRE_function variables are generated for faster execution
#source: Seeliger, Daniel, and Bert L. de Groot. "Atomic contacts in protein structures. A detailed analysis of atomic radii, packing, and overlaps." Proteins: Structure, Function, and Bioinformatics 68.3 (2007): 595-601.
SF=0.6
R = {'H0':1.19*SF,   'HAR':1.14*SF,  'HA':1.03*SF,   'H':1.05*SF,
     'HC':0.58*SF,   'HDR':0.67*SF,  'C':1.43*SF,    'CA':1.48*SF,
     'CH1E':1.92*SF, 'CH2E':1.89*SF, 'CH3E':1.81*SF, 'CR1E':1.81*SF,
     'CR1W':1.76*SF, 'C5':1.76*SF,   'C5W':1.86*SF,  'CW':1.74*SF,
     'CH2G':1.76*SF, 'CH2P':1.47*SF, 'CY':1.87*SF,   'CY2':1.63*SF,
     'CF':1.83*SF,   'CDR':1.69*SF,  'CR1H':1.75*SF, 'CRHH':1.63*SF,
     'O':1.41*SF,    'OC':1.33*SF,   'OH1':1.31*SF,  'NH1':1.37*SF,
     'NH2':1.45*SF,  'NH3':1.35*SF,  'NC2':1.45*SF,  'NHS':1.40*SF,
     'SM':1.79*SF,   'S':1.83*SF}

#determine radii of all atoms of all amino acids:
GLY_radii = {' N  ':R['NH1'],  ' CA ':R['CH2G'],
             ' C  ':R['C'],    ' O  ':R['O'], 
             ' H  ':R['H'],    ' HA2':R['HA'],
             ' HA3':R['HA']}

ALA_radii = {' N  ':R['NH1'],  ' CA ':R['CA'],
             ' C  ':R['C'],    ' O  ':R['O'],
             ' CB ':R['CH3E'], ' H  ':R['H'],
             ' HA ':R['HA'],   ' HB1':R['H0'],
             ' HB2':R['H0'],   ' HB3':R['H0']}

VAL_radii = {' N  ':R['NH1'],  ' CA ':R['CA'], 
             ' C  ':R['C'],    ' O  ':R['O'], 
             ' CB ':R['CH1E'], ' CG1':R['CH3E'], 
             ' CG2':R['CH3E'], ' H  ':R['H'], 
             ' HA ':R['HA'],   ' HB ':R['H0'], 
             'HG11':R['H0'],   'HG12':R['H0'], 
             'HG13':R['H0'],   'HG21':R['H0'], 
             'HG22':R['H0'],   'HG23':R['H0']}

LEU_radii = {' N  ':R['NH1'],  ' CA ':R['CA'], 
             ' C  ':R['C'],    ' O  ':R['O'], 
             ' CB ':R['CH2E'], ' CG ':R['CH1E'], 
             ' CD1':R['CH3E'], ' CD2':R['CH3E'], 
             ' H  ':R['H'],    ' HA ':R['HA'], 
             ' HB2':R['H0'],   ' HB3':R['H0'], 
             ' HG ':R['H0'],   'HD11':R['H0'], 
             'HD12':R['H0'],   'HD13':R['H0'], 
             'HD21':R['H0'],   'HD22':R['H0'], 
             'HD23':R['H0']}

ILE_radii = {' N  ':R['NH1'],  ' CA ':R['CA'], 
             ' C  ':R['C'],    ' O  ':R['O'], 
             ' CB ':R['CH1E'], ' CG1':R['CH2E'], 
             ' CG2':R['CH3E'], ' CD1':R['CH3E'], 
             ' H  ':R['H'],    ' HA ':R['HA'], 
             ' HB ':R['H0'],   'HG12':R['H0'], 
             'HG13':R['H0'],   'HG21':R['H0'], 
             'HG22':R['H0'],   'HG23':R['H0'], 
             'HD11':R['H0'],   'HD12':R['H0'], 
             'HD13':R['H0']}

MET_radii = {' N  ':R['NH1'],  ' CA ':R['CA'], 
             ' C  ':R['C'],    ' O  ':R['O'], 
             ' CB ':R['CH2E'], ' CG ':R['CH2E'], 
             ' SD ':R['SM'],   ' CE ':R['CH3E'], 
             ' H  ':R['H'],    ' HA ':R['HA'],  
             ' HB2':R['H0'],   ' HB3':R['H0'], 
             ' HG2':R['H0'],   ' HG3':R['H0'], 
             ' HE1':R['H0'],   ' HE2':R['H0'], 
             ' HE3':R['H0']}

PHE_radii = {' N  ':R['NH1'],  ' CA ':R['CA'], 
             ' C  ':R['C'],    ' O  ':R['O'], 
             ' CB ':R['CH2E'], ' CG ':R['CF'], 
             ' CD1':R['CR1E'], ' CD2':R['CR1E'], 
             ' CE1':R['CR1E'], ' CE2':R['CR1E'], 
             ' CZ ':R['CR1E'], ' H  ':R['H'], 
             ' HA ':R['HA'],   ' HB2':R['H0'], 
             ' HB3':R['H0'],   ' HD1':R['HAR'], 
             ' HD2':R['HAR'],  ' HE1':R['HAR'], 
             ' HE2':R['HAR'],  ' HZ ':R['HAR']}

TYR_radii = {' N  ':R['NH1'],  ' CA ':R['CA'], 
             ' C  ':R['C'],    ' O  ':R['O'], 
             ' CB ':R['CH2E'], ' CG ':R['CY'],
             ' CD1':R['CR1E'], ' CD2':R['CR1E'], 
             ' CE1':R['CR1E'], ' CE2':R['CR1E'],
             ' CZ ':R['CY2'],  ' OH ':R['OH1'], 
             ' H  ':R['H'],    ' HA ':R['HA'], 
             ' HB2':R['H0'],   ' HB3':R['H0'], 
             ' HD1':R['HAR'],  ' HD2':R['HAR'], 
             ' HE1':R['HAR'],  ' HE2':R['HAR'], 
             ' HH ':R['H']}

TRP_radii = {' N  ':R['NH1'],  ' CA ':R['CA'],
             ' C  ':R['C'],    ' O  ':R['O'],
             ' CB ':R['CH2E'], ' CG ':R['C5W'], 
             ' CD1':R['CR1E'], ' CD2':R['CR1E'], 
             ' NE1':R['NH1'],  ' CE2':R['CR1W'], 
             ' CE3':R['CR1W'], ' CZ2':R['CR1W'], 
             ' CZ3':R['CR1W'], ' CH2':R['CR1E'], 
             ' H  ':R['H'],    ' HA ':R['HA'], 
             ' HB2':R['H0'],   ' HB3':R['H0'], 
             ' HD1':R['HAR'],  ' HE1':R['H'], 
             ' HE3':R['HAR'],  ' HZ2':R['HAR'], 
             ' HZ3':R['HAR'],  ' HH2':R['HAR']}

PRO_radii = {' N  ':R['NH1'],  ' CA ':R['CA'], 
             ' C  ':R['C'],    ' O  ':R['O'], 
             ' CB ':R['CH2P'], ' CG ':R['CH2P'], 
             ' CD ':R['CH2P'], ' HA ':R['HA'], 
             ' HB2':R['H0'],   ' HB3':R['H0'], 
             ' HG2':R['H0'],   ' HG3':R['H0'], 
             ' HD2':R['H0'],   ' HD3':R['H0']}

SER_radii = {' N  ':R['NH1'],  ' CA ':R['CA'], 
             ' C  ':R['C'],    ' O  ':R['O'], 
             ' CB ':R['CH2E'], ' OG ':R['OH1'], 
             ' H  ':R['H'],    ' HA ':R['HA'], 
             ' HB2':R['H0'],   ' HB3':R['H0'], 
             ' HG ':R['H']}

CYS_radii = {' N  ':R['NH1'],  ' CA ':R['CA'], 
             ' C  ':R['C'],    ' O  ':R['O'], 
             ' CB ':R['CH2E'], ' SG ':R['S'], 
             ' H  ':R['H'],    ' HA ':R['HA'], 
             ' HB2':R['H0'],   ' HB3':R['H0'], 
             ' HG ':R['H0']}

THR_radii = {' N  ':R['NH1'],  ' CA ':R['CA'], 
             ' C  ':R['C'],    ' O  ':R['O'], 
             ' CB ':R['CH1E'], ' OG1':R['OH1'], 
             ' CG2':R['CH3E'], ' H  ':R['H'], 
             ' HA ':R['HA'],   ' HB ':R['H0'], 
             ' HG1':R['H'],    'HG21':R['H0'], 
             'HG22':R['H0'],   'HG23':R['H0']}
    
ASN_radii = {' N  ':R['NH1'],  ' CA ':R['CA'], 
             ' C  ':R['C'],    ' O  ':R['O'], 
             ' CB ':R['CH2E'], ' CG ':R['C'], 
             ' OD1':R['O'],    ' ND2':R['NH2'], 
             ' H  ':R['H'],    ' HA ':R['HA'], 
             ' HB2':R['H0'],   ' HB3':R['H0'], 
             'HD21':R['H'],    'HD22':R['H']}

GLN_radii = {' N  ':R['NH1'],  ' CA ':R['CA'], 
             ' C  ':R['C'],    ' O  ':R['O'], 
             ' CB ':R['CH2E'], ' CG ':R['CH2E'], 
             ' CD ':R['C'],    ' OE1':R['O'], 
             ' NE2':R['NH2'],  ' H  ':R['H'], 
             ' HA ':R['HA'],   ' HB2':R['H0'], 
             ' HB3':R['H0'],   ' HG2':R['H0'], 
             ' HG3':R['H0'],   'HE21':R['H'], 
             'HE22':R['H']}

ASP_radii = {' N  ':R['NH1'],  ' CA ':R['CA'], 
             ' C  ':R['C'],    ' O  ':R['O'], 
             ' CB ':R['CH2E'], ' CG ':R['C'], 
             ' OD1':R['O'],    ' OD2':R['O'], 
             ' H  ':R['H'],    ' HA ':R['HA'], 
             ' HB2':R['H0'],   ' HB3':R['H0']}

GLU_radii = {' N  ':R['NH1'],  ' CA ':R['CA'], 
             ' C  ':R['C'],    ' O  ':R['O'], 
             ' CB ':R['CH2E'], ' CG ':R['CH2E'], 
             ' CD ':R['C'],    ' OE1':R['O'], 
             ' OE2':R['O'],    ' H  ':R['H'], 
             ' HA ':R['HA'],   ' HB2':R['H0'], 
             ' HB3':R['H0'],   ' HG2':R['H0'], 
             ' HG3':R['H0']}

LYS_radii = {' N  ':R['NH1'],  ' CA ':R['CA'], 
             ' C  ':R['C'],    ' O  ':R['O'],
             ' CB ':R['CH2E'], ' CG ':R['CH2E'], 
             ' CD ':R['CH2E'], ' CE ':R['CH2E'], 
             ' NZ ':R['NH3'],  ' H  ':R['H'], 
             ' HA ':R['HA'],   ' HB2':R['H0'], 
             ' HB3':R['H0'],   ' HG2':R['H0'], 
             ' HG3':R['H0'],   ' HD2':R['H0'], 
             ' HD3':R['H0'],   ' HE2':R['H0'], 
             ' HE3':R['H0'],   ' HZ1':R['HC'], 
             ' HZ2':R['HC'],   ' HZ3':R['HC']}

ARG_radii = {' N  ':R['NH1'],  ' CA ':R['CA'], 
             ' C  ':R['C'],    ' O  ':R['O'],
             ' CB ':R['CH2E'], ' CG ':R['CH2E'], 
             ' CD ':R['CDR'],  ' NE ':R['NC2'], 
             ' CZ ':R['CR1E'], ' NH1':R['NH2'], 
             ' NH2':R['NH2'],  ' H  ':R['H'], 
             ' HA ':R['HA'],   ' HB2':R['H0'], 
             ' HB3':R['H0'],   ' HG2':R['H0'], 
             ' HG3':R['H0'],   ' HD2':R['HDR'], 
             ' HD3':R['HDR'],  ' HE ':R['HC'], 
             'HH11':R['HC'],   'HH12':R['HC'], 
             'HH21':R['HC'],   'HH22':R['HC']}

HIS_radii = {' N  ':R['NH1'],  ' CA ':R['CA'], 
             ' C  ':R['C'],    ' O  ':R['O'],
             ' CB ':R['CH2E'], ' CG ':R['CR1E'], 
             ' ND1':R['NHS'],  ' CD2':R['CR1H'], 
             ' CE1':R['CRHH'], ' NE2':R['NHS'], 
             ' H  ':R['H'],    ' HA ':R['HA'], 
             ' HB2':R['H0'],   ' HB3':R['H0'], 
             ' HD2':R['H0'],   ' HE1':R['H0']}

PTN_radii = {'GLY':GLY_radii, 'ALA':ALA_radii, 'VAL':VAL_radii, 'LEU':LEU_radii, 
             'ILE':ILE_radii, 'MET':MET_radii, 'PHE':PHE_radii, 'TYR':TYR_radii, 
             'TRP':TRP_radii, 'PRO':PRO_radii, 'SER':SER_radii, 'CYS':CYS_radii, 
             'THR':THR_radii, 'ASN':ASN_radii, 'GLN':GLN_radii, 'ASP':ASP_radii, 
             'GLU':GLU_radii, 'LYS':LYS_radii, 'ARG':ARG_radii, 'HIS':HIS_radii}

#ligand radii are based on approximations/assumptions taken from protein radii:
#whenever a choice is presnted, the smallest possible radius is used.
#it is assumed that AutoDock will penalise residues too close to the ligand:
LIG_radii = {'H':R['HC'], 'C':R['C'],  'N':R['NH3'],  'O':R['OH1'],
             'S':R['SM'], 'P':R['SM'], 'X':R['C']}

def clasher(residue_1, residue_2):
    
    #set atomic radii to a global variable.
    #This imports the radii given right before this function:
    global PTN_radii 
    global LIG_radii
    clash_flag = 0
    residue_list = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
    ligand_list = ['H', 'C', 'N', 'O', 'S', 'P']

    #identify residue names:
    residue_1_name = residue_1[0][17:20]
    residue_2_name = residue_2[0][17:20]

    for i in range(0,len(residue_1)):
        if len(str(residue_1[i])) != 81:
            print "warning: <clasher> has encountered an incomplete residue (1)."
            print "now exiting <clasher> function..."
            return "WARNING"
        residue_1_coords = []
        residue_1_coords.append(float(residue_1[i][30:38]))
        residue_1_coords.append(float(residue_1[i][38:46]))
        residue_1_coords.append(float(residue_1[i][46:54]))
        #determine whether selected atom (i) belongs to backbone:
        residue_1_atom = copy.copy(residue_1[i][12:16])
        residue_1_backbone_flag = 0
        if residue_1_atom == " N  " or residue_1_atom == " CA " or residue_1_atom == " C  " or residue_1_atom == " O  ":
            residue_1_backbone_flag = 1
        
        for j in range(0,len(residue_2)):
            if len(str(residue_2[j])) != 81:
                print "warning: <clasher> has encountered an incomplete residue (2)."
                print "now exiting <clasher> function..."
                return "WARNING"
        #determine whether selected atom (j) belongs to backbone:
            residue_2_atom = copy.copy(residue_2[j][12:16])
            residue_2_backbone_flag = 0
            if residue_2_atom == " N  " or residue_2_atom == " CA " or residue_2_atom == " C  " or residue_2_atom == " O  ":
                residue_2_backbone_flag = 1
            #IGNORE backbone-backbone clashes:
            if residue_1_backbone_flag != 1 or residue_2_backbone_flag != 1:  
                #calculate the minimum acceptable inter-atomic radius:
                #read atomic X,Y,Z co-ordinates:
                residue_2_coords = []
                residue_2_coords.append(float(residue_2[j][30:38]))
                residue_2_coords.append(float(residue_2[j][38:46]))
                residue_2_coords.append(float(residue_2[j][46:54]))
                #calculate interatomic distance:
                distance = math.sqrt(pow(residue_1_coords[0]-residue_2_coords[0],2) +pow(residue_1_coords[1]-residue_2_coords[1],2) +pow(residue_1_coords[2]-residue_2_coords[2],2))

                residue_2_atom = residue_2[j][12:16]
                #extract atomic radii from double dictionary (radii):
                #check if inputted residue is protein or ligand:
                if residue_1_name in residue_list:
                    #identify atom types:
                    residue_1_atom = residue_1[i][12:16]
                    residue_1_radius = PTN_radii[residue_1_name][residue_1_atom]
                else:
                    residue_1_atom = residue_1[i][77]
                    if residue_1_atom not in ligand_list:
                        residue_1_atom = 'X'
                    residue_1_radius = LIG_radii[residue_1_atom]
                if residue_2_name in residue_list:
                    #identify atom types:
                    residue_2_atom = residue_2[j][12:16]
                    residue_2_radius = PTN_radii[residue_2_name][residue_2_atom]
                else:
                    residue_2_atom = residue_2[j][77]
                    if residue_2_atom not in ligand_list:
                        residue_2_atom = 'X'
                    residue_2_radius = LIG_radii[residue_2_atom]
                #calculate minimum permissible interatomic distance:
                min_distance = residue_1_radius +residue_2_radius
                if distance <= min_distance:
                    clash_flag = 1
    #print "CLASH: "+str(clash_flag)+" | "+residue_1_name+":"+residue_1_atom+"| "+residue_2_name+":"+residue_2_atom
    return clash_flag

########################################
#      FUNCTION: COOLING SCHEDULE      #
########################################

#Select cooling schedule based on user input:
def cooling_scheduler(i, cycles, cooling_schedule, start_temperature, end_temperature):
    i = float(i)
    T0 = float(start_temperature)
    Tn = float(end_temperature)
    N = float(cycles)
    
    if cooling_schedule == 0:
        Ti = T0 -i*(T0-Tn)/N
    
    if cooling_schedule == 1:
        Ti = T0*(Tn/T0)**(i/N)
    
    if cooling_schedule == 2:
        A = (T0-Tn)*float(N+1)/N
        B = T0 -A
        Ti = A/(i+1) +B
    
    if cooling_schedule == 3:
        print "warning: cooling_schedule '3' does not work as described"
        print "switching to cooling_schedule '5'"
        cooling_schedule = 5
    
    if cooling_schedule == 4:
        Ti = (T0-Tn)/(1+math.exp(0.01*(i-N/2))) +Tn;
    
    if cooling_schedule == 5:
        Ti = 0.5*(T0 -Tn)*(1+math.cos(i*math.pi/N)) +Tn
    
    if cooling_schedule == 6:
        Ti = 0.5*(T0-Tn)*(1-math.tanh(i*10/N-5)) +Tn;
    
    if cooling_schedule == 7:
        Ti = (T0-Tn)/math.cosh(i*10/N) +Tn;
    
    if cooling_schedule == 8:
        A = (1/N)*math.log(T0/Tn)
        Ti = T0*math.exp(-A*i)

    if cooling_schedule == 9:
        A = (1/N**2)*math.log(T0/Tn)
        Ti = T0*math.exp(-A*i**2);
    
    return Ti

########################################
#          AUTODOCK FUNCTIONS          #
########################################

#prepare ligand for AutoDock4.1 scoring:
def ligand_preparation(ligand, output):
    ligand_filename =  None
    verbose = None
    add_bonds = False
    repairs = ""
    charges_to_add = 'gasteiger'
    preserve_charge_types=''
    cleanup  = "nphs_lps"
    allowed_bonds = "backbone"
    root = 'auto'
    outputfilename = output
    check_for_fragments = False
    bonds_to_inactivate = ""
    inactivate_all_torsions = False
    attach_nonbonded_fragments = False
    attach_singletons = False
    mode = 'automatic'
    dict = None

    mols = Read(ligand)
    mol = mols[0]
    mol.buildBondsByDistance()
        
    LPO = AD4LigandPreparation(mol, mode, repairs, charges_to_add, cleanup,
          allowed_bonds, root, outputfilename=output, dict=dict,
          check_for_fragments=check_for_fragments,
          bonds_to_inactivate=bonds_to_inactivate,
          inactivate_all_torsions=inactivate_all_torsions,
          attach_nonbonded_fragments=attach_nonbonded_fragments,
          attach_singletons=attach_singletons)

#prepare protein for AutoDock4.1 scoring:
def protein_preparation(receptor, output):
    receptor_filename =  receptor
    verbose = None
    repairs = 'checkhydrogens'
    charges_to_add = 'gasteiger'
    preserve_charge_types=None
    cleanup  = "nphs_lps_waters_nonstdres"
    mode = 'automatic'
    delete_single_nonstd_residues = None
    dictionary = None

    mols = Read(receptor)
    mol = mols[0]
    RPO = AD4ReceptorPreparation(mol, mode, repairs, charges_to_add, cleanup,
          delete_single_nonstd_residues=delete_single_nonstd_residues,
          outputfilename=output,dict=dictionary)

#score protein-ligand complex using AutoDock4.1 (pdbqt files only):
def autodock_scoring(receptor, ligand):
    receptorfilename =  receptor
    ligandfilename =  ligand
    write_file_mode = False
    parameter_library_filename = None
    exclude_torsFreeEnergy = False
    verbose = None
    ad_scorer = AutoDock41Scorer(exclude_torsFreeEnergy=exclude_torsFreeEnergy)
    supported_types = ad_scorer.supported_types
    receptor = Read(receptorfilename)[0]
    receptor.buildBondsByDistance()
    ligand = Read(ligandfilename)[0]
    ligand.buildBondsByDistance()

    ms = MolecularSystem()
    ms.add_entities(receptor.allAtoms)
    ms.add_entities(ligand.allAtoms)
    ad_scorer.set_molecular_system(ms)
    #get the scores, score per term:
    [estat, hb, vdw ,dsolv] = ad_scorer.get_score_per_term()
    torsEnrg = ligand.TORSDOF * ad_scorer.tors_weight
    score = estat +hb +vdw +dsolv +torsEnrg
    output_score = {'score':score, 'estat':estat, 'hb':hb, 'vdw':vdw, 'dsolv,':dsolv, 'torsEnrg':torsEnrg}
    
    return output_score


#score protein-ligand complex using AutoDock4.1:
def autoscore(protein, ligand_pdbqt, autoscore_directory):
    
    #definitions:
    anchor = os.getcwd()
    
    #cleanup/make autoscore_directory:
    if os.path.isdir(autoscore_directory) == True:
        shutil.rmtree(autoscore_directory)
    os.makedirs(autoscore_directory)
    
    #output protein:
    OP_protein = open(autoscore_directory+'/OP_protein.pdb', 'w+')
    for residue in range(0, len(protein)):
        for line in range(0, len(protein[residue].IP_residue)):
            OP_protein.write("%s" % protein[residue].IP_residue[line])
    
    OP_protein.close()
    
    #create protein.pdbqt file:
    protein_preparation(autoscore_directory+'/OP_protein.pdb', autoscore_directory+'/OP_protein.pdbqt')
       
    #score protein/ligand complex (as pdbqt files):
    score = autodock_scoring(autoscore_directory+'/OP_protein.pdbqt', ligand_pdbqt)
    
    #cleanup:
    shutil.rmtree(autoscore_directory)
    return score

########################################
#   OUTPUT POCKET AT EACH ITERATION    #
########################################

def iterpocket_maker(protein, i, OLD_score, working_directory):
    iterpocket_directory = str(working_directory)+"/iterpocket_output"
    if not os.path.exists(iterpocket_directory):
        os.makedirs(iterpocket_directory)
    number_string = "%06d" % i
    ITER_protein = open(str(iterpocket_directory)+'/ITER_protein_'+number_string+'.pdb', 'w+')
    ITER_protein.write("REMARK: iteration: "+number_string+", autodsock score: "+str(OLD_score)+"\n")
    for residue in range(0, len(protein)):
        for line in range(0, len(protein[residue].IP_residue)):
            ITER_protein.write("%s" % protein[residue].IP_residue[line])

    ITER_protein.close()
    return

########################################
#         USAGE/ERROR MESSAGE          #
########################################

#ensure the correct number of arguments are provided:
if len(sys.argv)<13 or len(sys.argv)>17:
    print "usage: Pocketanneal.py -P <protein.pdb> -p <pocket.pdb> -d <database> -n 1000 -a 2 -o <output_directory>"
    print "MANDATORY FLAGS:"
    print "  -P/--protein: input entire protein structure (after formatting)."
    print "  -p/--pocket: input ligand + surrounding pocket residues."
    print "  -d/--database: path to Autoanneal database."
    print "  -n/--iterations: number of simulated annealing iterations."
    print "  -a/--annealer: choose an annealing schedule (0 to 9, recommended 2)."
    print "  -o/--output: choose an output directory."
    print "OPTIONAL FLAGS:"
    print "  -r/--report: output continuous progress reports (for GUI/webserver only)."
    print "  -q/--iterpocket: output the interface design after every annealing iteration (debugging only)."
    print "  -s/--seed: define a random-number seed for reproducible runs."
    sys.exit()

########################################
#              DATA INPUT              #
########################################

#definitions:
opts, args = getopt.getopt(sys.argv[1:], "P:p:d:n:a:o:rqs:", ["protein=", "pocket=", "database=", "iterations=", "annealer=", "output=", "report", "iterpocket", "seed="])

#'report_flag' works as both a flag and counter:
report_flag = 0
iterpocket_flag = 0
seed = random.randint(0, sys.maxint)
for o,a in opts:
    if o == "-P" or o == "--protein":
        IP_protein = a
    if o == "-p" or o == "--pocket":
        IP_pocket = a
    if o == "-d" or o == "--database":
        IP_database = a
    if o == "-n" or o == "--iterations":
        cycles = int(a)
        if cycles < 10:
            print "warning: too few iterations selected. Defaulting to 10."
            cycles = 10
    if o == "-a" or o == "--annealer":
        cooling_schedule = int(a)
    if o == "-o" or o == "--output":
        working_directory = a
    if o == "-r" or o == "--report":
        print "progress-report generation initiated..."
        report_flag = 1
    if o == "-q" or o == "--iterpocket":
        print "outputting interface after each iteration..."
        iterpocket_flag = 1
    if o == "-s" or o == "--seed":
        seed = a

#Use the random-number seed to initialize random-number generation:
random.seed(seed)

#Import class 'PDB_residue'. Useful for rotamer generation.
execfile(IP_database+'/classes/PDB_residue.py')

#Import class 'Dunbrack'. Useful for easy access to Dunbrack's rotamer library.
execfile(IP_database+'/classes/Dunbrack.py')

#make working_directory (if absent):
if os.path.isdir(working_directory) == False:
    os.makedirs(working_directory)
#check contents of working_directory. Terminate if not completely empty:
if len(os.listdir(working_directory)) != 0:
    print "error: --output directory is not empty."
    print "Autoanneal.py will now terminate..."
    sys.exit()

#use a dictionary for quick access to naive residues based on residue name:
naive_residue = {}
residue_list = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]

for i in range(0, len(residue_list)):
    naive_file = open(IP_database+"/PDBnaive_files/"+residue_list[i]+"_naiveH.pdb")
    naive_residue[residue_list[i]] = naive_file.readlines()

#extract inputPDB_file, store as list of objects:
PDB_file = open(IP_protein)
PDB_data = PDB_file.readlines()
old_residue_ID = PDB_data[0][21:27]
new_residue_ID = []
residue_store = []
protein = []
for i in range(0, len(PDB_data)):
    if PDB_data[i][0:6] == "ATOM  ":
        #extract individual amino acid residues, store as 2D lists:
        new_residue_ID = copy.copy(PDB_data[i][21:27])
        if new_residue_ID != old_residue_ID or i == len(PDB_data)-1:
            if i == len(PDB_data)-1:
                residue_store.append(copy.copy(PDB_data[i]))
            protein.append(PDB_residue(residue_store, naive_residue[residue_store[0][17:20]]))
            residue_store = []
        residue_store.append(copy.copy(PDB_data[i]))
        old_residue_ID = new_residue_ID

#extract pocket:
PDB_pocket_file = open(IP_pocket)
PDB_pocket_data = PDB_pocket_file.readlines()
old_residue_ID = PDB_pocket_data[0][21:27]
new_residue_ID = []
residue_store = []
pocket = []
for i in range(0, len(PDB_pocket_data)):
    if PDB_pocket_data[i][0:6] == "ATOM  ":
        #extract individual amino acid residues, store as 2D lists:
        new_residue_ID = copy.copy(PDB_pocket_data[i][21:27])
        if new_residue_ID != old_residue_ID or i == len(PDB_pocket_data)-1 or (i != len(PDB_pocket_data)-1 and PDB_pocket_data[i+1][0:6] != "ATOM  "):
            pocket.append(residue_store)
            residue_store = []
        residue_store.append(copy.copy(PDB_pocket_data[i]))
        old_residue_ID = new_residue_ID

#mark pocket residues on PDB_residue object:
for i in range(0, len(protein)):
    for j in range(0, len(pocket)):
        if protein[i].IP_residue[0][21:27] == pocket[j][0][21:27]:
            protein[i].pocket_flag = "YES"
            protein[i].Rama()
            PDB_residue.pocket_list.append(i)

#extract ligand. Store as list (no special object):
#repair ligand. In some cases, ligands do not contain 80 characters & atom names:
ligand = []
spacer="                                "
counter = 0
for i in range(0, len(PDB_pocket_data)):
    if PDB_pocket_data[i][0:6] == "HETATM":
        ligand.append(PDB_pocket_data[i])
        if len(PDB_pocket_data[i]) != 81:
            ligand_atom = ligand[counter][13]
            ligand[counter] = ligand[counter]+spacer
            ligand[counter] = ligand[counter][0:77]+ligand_atom+"  \n"
        counter = counter+1

#output ligand:
OP_ligand = open(str(working_directory)+'/OP_ligand.pdb', 'w+')
for line in range(0, len(ligand)):
     OP_ligand.write("%s" % ligand[line])
OP_ligand.close()

#create ligand.pdbqt file:
ligand_preparation(working_directory+"/OP_ligand.pdb", working_directory+"/OP_ligand.pdbqt")
#calculate 'clash_zones' for all pocket residues.
#'clash zones' include all atoms within 7 Angstroms from a selected residue:
distance_threshold = 7
for i in range(0, len(protein)):
    if protein[i].pocket_flag == "YES":
        for j in range(0, len(protein)):
            clash_flag = 0
            if i != j:
                for a in range(0, len(protein[i].IP_residue)):
                    for b in range(0, len(protein[j].IP_residue)):
                        distance = math.sqrt(math.pow(float(protein[i].IP_residue[a][30:38]) -float(protein[j].IP_residue[b][30:38]),2) +math.pow(float(protein[i].IP_residue[a][38:46]) -float(protein[j].IP_residue[b][38:46]),2) +math.pow(float(protein[i].IP_residue[a][46:54]) -float(protein[j].IP_residue[b][46:54]),2))
                        if distance <= distance_threshold:
                            clash_flag = 1
            if clash_flag == 1:
                protein[i].clash_zone.append(j)

#extract Dunbrack's rotamer library (2010). Store in a dictionary-based class:
Dunbrack = Dunbrack_library(IP_database+"/dunbrack/dunbrack_rotamers_2010.csv")

#write first progress report
if report_flag >= 1:
    report_FIRST_file = open(str(working_directory)+"/report.txt", 'w+')
    report_FIRST_file.write("All inputs processed. Initiating simulated annealing.\n")
    report_FIRST_file.close()
    report_flag = report_flag+1

########################################
#         SIMULATED ANNEALING          #
########################################

#preserve initial [0] residue for possible future PRO to XXX mutations:
start_residue = copy.copy(protein[0])
if start_residue.name == "PRO":
    print "warning: N-terminal PRO residue has been selected as part of the pocket."
    print "ramachandaran angles are incalculable without a previous C=O group."
    print "deselecting N-terminal PRO from pocket..."
    PDB_residue.pocket_list.remove(0)

#design a customised protein-alteration function.
#Unlike the 'Alter' function of the 'PDB_residue' class, this function is
#especially designed for the simulated annealing code-block of 'Autoanneal.py'.
def SA_alter(res_number, alteration):
    global protein
    global naive_residue
    global Dunbrack
    chi = Dunbrack.Surprise_SD(alteration, protein[res_number].phi, protein[res_number].psi)
    if res_number == 0 and protein[res_number].name == "PRO" and alteration != "PRO":
        print "note: PRO to XXX mutation requested at N-terminal residue"
        protein[res_number] = copy.copy(start_residue)
        protein[res_number].Alter(naive_residue[alteration], protein[res_number-1].IP_residue, chi)
    else:
        protein[res_number].Alter(naive_residue[alteration], protein[res_number-1].IP_residue, chi)

OLD_score = {'score':float("inf"), 'estat':float("inf"), 'hb':float("inf"), 'vdw':float("inf"), 'dsolv,':float("inf"), 'torsEnrg':float("inf")}
NEW_score = {}
BEST_score = copy.copy(OLD_score)
score_log = []
protein_log = []
autoscore_directory = str(working_directory)+str("/autoscore")
start_temperature = 1
end_temperature = 0.01
temperature = start_temperature

#score given pocket, determine starting score:
START_score = autoscore(protein, working_directory+"/OP_ligand.pdbqt", autoscore_directory)
print "input pocket has an Autodock4.1 score of: "+str(START_score['score'])

#mutate all pocket residues to alanine, also insert glycine as backbone proxy:
for i in range(0, len(protein)):
    protein[i].Backup()
    SA_alter(i, "GLY")
    protein[i].backbone_residue = copy.copy(protein[i].IP_residue)
    protein[i].Restore()
    if protein[i].pocket_flag == "YES":
        SA_alter(i, "ALA")

#score alanine pocket. This score will act as a baseline:
ALA_score = autoscore(protein, working_directory+"/OP_ligand.pdbqt", autoscore_directory)
print "ALA-nated pocket has an Autodock4.1 score of: "+str(ALA_score['score'])

#iterate simulated annealing thorugh user-specified number of cycles:
#a 'while' loop rather than the traditional 'for' loop is used here.
#this is to account for conditional 'i' decrements:
i = 0
while i<cycles:
    i = i+1
    
    #write a progress report after every 10% job-completion:
    if report_flag >= 1 and i%(cycles/10) == 0:
        report_SA_file = open(str(working_directory)+"/report.txt", 'a')
        report_SA_file.write("Simulated annealing: "+str(i)+" iterations completed.\n")
        report_SA_file.close()
        report_flag = report_flag+1
    
    #I'm using cooling-schedule #2: http://www.btluke.com/simanf1.html
    #temperature ranges from 1:0.01
    temperature = cooling_scheduler(i, cycles, cooling_schedule, start_temperature, end_temperature)
    print ""
    print "iteration: "+str(i)+"/"+str(cycles)+". Temperature: "+str(temperature)+". Starting score: "+str(OLD_score['score'])

    #backup protein:
    for j in range(0, len(protein)):
        if protein[j].pocket_flag == "YES":
            protein[j].Backup()
    
    #randomly select a pocket residue:
    chosen_residue = random.choice(PDB_residue.pocket_list)
    #randomly select a mutation (change this line later if you want weighted mutations):
    random_mutant = random.choice(residue_list)
    #mutate selected residue:
    print "mutating \""+protein[chosen_residue].IP_residue[0][17:26]+"\" to \""+random_mutant+"\""
    SA_alter(chosen_residue, random_mutant)

    #(check whether any clashes exist with newly mutated residue)
    #(clashes may be protein-protein or protein-ligand clashes:)
    clash_flag = 0
    warning_flag = 0
    for j in range(0, len(protein[chosen_residue].clash_zone)):

        #Detect incomplete residues, reset protein if any:
        clasher_output = clasher(protein[chosen_residue].IP_residue, protein[protein[chosen_residue].clash_zone[j]].IP_residue)
        if clasher_output == "WARNING":
            warning_flag = 1
            break
        
        #detect clashes:
        if clasher_output == 1:
            clash_flag = 1
            break
    
    if warning_flag == 1:
        print "<main> will end iteration and reset protein."
        protein[chosen_residue].Restore()
        i = i-1
        continue
    
    if clasher(protein[chosen_residue].IP_residue, ligand) == 1:
        clash_flag = 1
    
    if clash_flag == 1:
        #(check whether residue clashes with backbone/ligand/non-pocket residues)
        backbone_clash_flag = 0
        for j in range(0, len(protein[chosen_residue].clash_zone)):
            if clasher(protein[chosen_residue].IP_residue, protein[protein[chosen_residue].clash_zone[j]].backbone_residue) == 1:
                backbone_clash_flag = 1
            if clasher(protein[chosen_residue].IP_residue, ligand) == 1:
                backbone_clash_flag = 1
            if protein[protein[chosen_residue].clash_zone[j]].pocket_flag == "NO":
                backbone_clash_flag = 1
        #IF the residue clashes with the backbone/ligand/non-pocket residues:
        if backbone_clash_flag == 1:
            #END iteration:
            protein[chosen_residue].Restore()
            score_log.append(OLD_score)
            if iterpocket_flag == 1:
                iterpocket_maker(protein, i, OLD_score, working_directory)
            continue
        #IF the residue does NOT clash with the backbone/ligand/non-pocket residues:
        else:
            #IF the 'annealing threshold' has NOT been met:
            if random.random() > temperature:
                #END iteration:
                protein[chosen_residue].Restore()
                score_log.append(OLD_score)
                if iterpocket_flag == 1:
                    iterpocket_maker(protein, i, OLD_score, working_directory)
                continue
            #IF the 'annealing threshold' has been met:
            else:
                #(mutate clashing residues to alanine):
                for k in range(0, len(protein[chosen_residue].clash_zone)):
                    if clasher(protein[chosen_residue].IP_residue, protein[protein[chosen_residue].clash_zone[k]].backbone_residue) == 1:
                        SA_alter(k, "ALA")
                #output new pocket:
                NEW_protein = open(str(working_directory)+'/NEW_protein.pdb', 'w+')
                for residue in range(0, len(protein)):
                    for line in range(0, len(protein[residue].IP_residue)):
                        NEW_protein.write("%s" % protein[residue].IP_residue[line])
                #score new pocket:
                NEW_score = autoscore(protein, working_directory+"/OP_ligand.pdbqt", autoscore_directory)
                score_log.append(copy.copy(NEW_score))
                if iterpocket_flag == 1:
                    iterpocket_maker(protein, i, OLD_score, working_directory)
                OLD_score = copy.copy(NEW_score)
                #END iteration
                continue
    #IF NO clashes of any sort exist:
    else:
        #output new pocket:
        NEW_protein = open(str(working_directory)+'/NEW_protein.pdb', 'w+')
        for residue in range(0, len(protein)):
            for line in range(0, len(protein[residue].IP_residue)):
                NEW_protein.write("%s" % protein[residue].IP_residue[line])
        #score new pocket:
        NEW_score = autoscore(protein, working_directory+"/OP_ligand.pdbqt", autoscore_directory)
        #IF the new score is WORSE than the old score:
        if float(NEW_score['score']) > float(OLD_score['score']):
            #IF the 'annealing threshold' has NOT been met:
            if random.random() > temperature:
                #END iteration:
                protein[chosen_residue].Restore()
                score_log.append(OLD_score)
                if iterpocket_flag == 1:
                    iterpocket_maker(protein, i, OLD_score, working_directory)
                continue
            #IF the annealing threshold has been met:
            else:
                #END iteration:
                OLD_score = copy.copy(NEW_score)
                score_log.append(OLD_score)
                if iterpocket_flag == 1:
                    iterpocket_maker(protein, i, OLD_score, working_directory)
                continue
        #IF the new score is BETTER than the old score:
        else:
            #IF the new score is BETTER than the best score:
            if NEW_score['score'] < BEST_score['score']:
                BEST_score = copy.copy(NEW_score)
                for j in range(0, len(protein)):
                    protein[j].best_residue = copy.copy(protein[j].IP_residue)
            #END iteration:
            OLD_score = copy.copy(NEW_score)
            score_log.append(OLD_score)
            if iterpocket_flag == 1:
                iterpocket_maker(protein, i, OLD_score, working_directory)
            continue

########################################
#              DATA OUTPUT             #
########################################

#output BEST_protein:
BEST_protein = open(str(working_directory)+'/BEST_protein.pdb', 'w+')
for residue in range(0, len(protein)):
    for line in range(0, len(protein[residue].best_residue)):
        BEST_protein.write("%s" % protein[residue].best_residue[line])

BEST_protein.close()

#output BEST_score:
BEST_score_file = open(working_directory+"/BEST_score.csv", "w+")
BEST_score_file.write("%s\n" % BEST_score)

BEST_score_file.close()

#output score log:
score_log_file = open(str(working_directory)+'/score_log.csv', 'w+')
for i in range(0, len(score_log)):
    score_log_file.write("%s\n" % score_log[i]['score'])

score_log_file.close()

#format XY_data:
XY_data = []
for i in range(0, len(score_log)):
    if score_log[i]['score'] < float('inf'):
        XY_data.append([i, score_log[i]['score']])

#create graph of score_log:
grapher(XY_data, str(working_directory)+"/score_log.png", str(IP_database))

#write last progress report
if report_flag >= 1:
    report_LAST_file = open(str(working_directory)+"/report.txt", 'a')
    report_LAST_file.write("This is Autoanneal.py signing off.\n")
    report_LAST_file.close()
    report_flag = report_flag+1

#end script:
sys.exit()















