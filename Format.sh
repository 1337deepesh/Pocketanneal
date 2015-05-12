#!/bin/bash

#residue-specific functions for checking completeness:
function GLY {
    function_input=$1
    counter=0

    val=`cat $function_input|grep '^.\{12\} N  '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CA '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} C  '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} O  '|wc -l`
    counter=`echo $counter+$val|bc`

    if [[ $counter -ne 4 ]]
    then
        #The residue is incomplete:
        echo "1"
    else
        #The residue is complete:
        echo "0"
    fi
    }

function ALA {
    function_input=$1
    counter=0

    val=`cat $function_input|grep '^.\{12\} N  '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CA '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} C  '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} O  '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CB '|wc -l`
    counter=`echo $counter+$val|bc`

    if [[ $counter -ne 5 ]]
    then
        #The residue is incomplete:
        echo "1"
    else
        #The residue is complete:
        echo "0"
    fi
    }

function VAL {
    function_input=$1
    counter=0

    val=`cat $function_input|grep '^.\{12\} N  '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CA '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} C  '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} O  '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CB '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CG1'|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CG2'|wc -l`
    counter=`echo $counter+$val|bc`

    if [[ $counter -ne 7 ]]
    then
        #The residue is incomplete:
        echo "1"
    else
        #The residue is complete:
        echo "0"
    fi
    }

function LEU {
    function_input=$1
    counter=0

    val=`cat $function_input|grep '^.\{12\} N  '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CA '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} C  '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} O  '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CB '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CG '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CD1'|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CD2'|wc -l`
    counter=`echo $counter+$val|bc`

    if [[ $counter -ne 8 ]]
    then
        #The residue is incomplete:
        echo "1"
    else
        #The residue is complete:
        echo "0"
    fi
    }

function ILE {
    function_input=$1
    counter=0

    val=`cat $function_input|grep '^.\{12\} N  '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CA '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} C  '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} O  '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CB '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CG1'|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CG2'|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CD1'|wc -l`
    counter=`echo $counter+$val|bc`

    if [[ $counter -ne 8 ]]
    then
        #The residue is incomplete:
        echo "1"
    else
        #The residue is complete:
        echo "0"
    fi
    }

function MET {
    function_input=$1
    counter=0

    val=`cat $function_input|grep '^.\{12\} N  '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CA '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} C  '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} O  '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CB '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CG '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} SD '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CE '|wc -l`
    counter=`echo $counter+$val|bc`

    if [[ $counter -ne 8 ]]
    then
        #The residue is incomplete:
        echo "1"
    else
        #The residue is complete:
        echo "0"
    fi
    }

function PHE {
    function_input=$1
    counter=0

    val=`cat $function_input|grep '^.\{12\} N  '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CA '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} C  '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} O  '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CB '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CG '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CD1'|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CD2'|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CE1'|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CE2'|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CZ '|wc -l`
    counter=`echo $counter+$val|bc`

    if [[ $counter -ne 11 ]]
    then
        #The residue is incomplete:
        echo "1"
    else
        #The residue is complete:
        echo "0"
    fi
    }

function TYR {
    function_input=$1
    counter=0

    val=`cat $function_input|grep '^.\{12\} N  '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CA '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} C  '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} O  '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CB '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CG '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CD1'|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CD2'|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CE1'|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CE2'|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CZ '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} OH '|wc -l`
    counter=`echo $counter+$val|bc`

    if [[ $counter -ne 12 ]]
    then
        #The residue is incomplete:
        echo "1"
    else
        #The residue is complete:
        echo "0"
    fi
    }

function TRP {
    function_input=$1
    counter=0

    val=`cat $function_input|grep '^.\{12\} N  '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CA '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} C  '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} O  '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CB '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CG '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CD1'|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CD2'|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} NE1'|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CE2'|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CE3'|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CZ2'|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CZ3'|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CH2'|wc -l`
    counter=`echo $counter+$val|bc`

    if [[ $counter -ne 14 ]]
    then
        #The residue is incomplete:
        echo "1"
    else
        #The residue is complete:
        echo "0"
    fi
    }

function PRO {
    function_input=$1
    counter=0

    val=`cat $function_input|grep '^.\{12\} N  '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CA '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} C  '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} O  '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CB '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CG '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CD '|wc -l`
    counter=`echo $counter+$val|bc`

    if [[ $counter -ne 7 ]]
    then
        #The residue is incomplete:
        echo "1"
    else
        #The residue is complete:
        echo "0"
    fi
    }

function SER {
    function_input=$1
    counter=0

    val=`cat $function_input|grep '^.\{12\} N  '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CA '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} C  '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} O  '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CB '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} OG '|wc -l`
    counter=`echo $counter+$val|bc`

    if [[ $counter -ne 6 ]]
    then
        #The residue is incomplete:
        echo "1"
    else
        #The residue is complete:
        echo "0"
    fi
    }

function CYS {
    function_input=$1
    counter=0

    val=`cat $function_input|grep '^.\{12\} N  '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CA '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} C  '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} O  '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CB '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} SG '|wc -l`
    counter=`echo $counter+$val|bc`

    if [[ $counter -ne 6 ]]
    then
        #The residue is incomplete:
        echo "1"
    else
        #The residue is complete:
        echo "0"
    fi
    }

function THR {
    function_input=$1
    counter=0

    val=`cat $function_input|grep '^.\{12\} N  '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CA '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} C  '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} O  '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CB '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} OG1'|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CG2 '|wc -l`
    counter=`echo $counter+$val|bc`

    if [[ $counter -ne 7 ]]
    then
        #The residue is incomplete:
        echo "1"
    else
        #The residue is complete:
        echo "0"
    fi
    }

function ASN {
    function_input=$1
    counter=0

    val=`cat $function_input|grep '^.\{12\} N  '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CA '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} C  '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} O  '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CB '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CG '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} OD1'|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} ND2'|wc -l`
    counter=`echo $counter+$val|bc`

    if [[ $counter -ne 8 ]]
    then
        #The residue is incomplete:
        echo "1"
    else
        #The residue is complete:
        echo "0"
    fi
    }

function GLN {
    function_input=$1
    counter=0

    val=`cat $function_input|grep '^.\{12\} N  '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CA '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} C  '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} O  '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CB '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CG '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CD '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} OE1'|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} NE2'|wc -l`
    counter=`echo $counter+$val|bc`

    if [[ $counter -ne 9 ]]
    then
        #The residue is incomplete:
        echo "1"
    else
        #The residue is complete:
        echo "0"
    fi
    }

function ASP {
    function_input=$1
    counter=0

    val=`cat $function_input|grep '^.\{12\} N  '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CA '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} C  '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} O  '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CB '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CG '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} OD1'|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} OD2'|wc -l`
    counter=`echo $counter+$val|bc`

    if [[ $counter -ne 8 ]]
    then
        #The residue is incomplete:
        echo "1"
    else
        #The residue is complete:
        echo "0"
    fi
    }

function GLU {
    function_input=$1
    counter=0

    val=`cat $function_input|grep '^.\{12\} N  '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CA '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} C  '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} O  '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CB '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CG '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CD '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} OE1'|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} OE2'|wc -l`
    counter=`echo $counter+$val|bc`

    if [[ $counter -ne 9 ]]
    then
        #The residue is incomplete:
        echo "1"
    else
        #The residue is complete:
        echo "0"
    fi
    }

function LYS {
    function_input=$1
    counter=0

    val=`cat $function_input|grep '^.\{12\} N  '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CA '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} C  '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} O  '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CB '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CG '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CD '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CE '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} NZ '|wc -l`
    counter=`echo $counter+$val|bc`

    if [[ $counter -ne 9 ]]
    then
        #The residue is incomplete:
        echo "1"
    else
        #The residue is complete:
        echo "0"
    fi
    }

function ARG {
    function_input=$1
    counter=0

    val=`cat $function_input|grep '^.\{12\} N  '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CA '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} C  '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} O  '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CB '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CG '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CD '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} NE '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CZ '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} NH1'|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} NH2'|wc -l`
    counter=`echo $counter+$val|bc`

    if [[ $counter -ne 11 ]]
    then
        #The residue is incomplete:
        echo "1"
    else
        #The residue is complete:
        echo "0"
    fi
    }

function HIS {
    function_input=$1
    counter=0

    val=`cat $function_input|grep '^.\{12\} N  '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CA '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} C  '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} O  '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CB '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CG '|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} ND1'|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CD2'|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} CE1'|wc -l`
    counter=`echo $counter+$val|bc`
    val=`cat $function_input|grep '^.\{12\} NE2'|wc -l`
    counter=`echo $counter+$val|bc`

    if [[ $counter -ne 10 ]]
    then
        #The residue is incomplete:
        echo "1"
    else
        #The residue is complete:
        echo "0"
    fi
    }

#master function for checking completeness:
function COMPLETE {
    
    #definitions:
    residue_name=$1
    residue_file=$2
    incomplete_flag=0
    
    if [ $residue_name == "GLY" ]
    then
        incomplete_flag=`GLY $residue_file`
    elif [ $residue_name == "ALA" ]
    then
        incomplete_flag=`ALA $residue_file`
    elif [ $residue_name == "VAL" ]
    then
        incomplete_flag=`VAL $residue_file`
    elif [ $residue_name == "LEU" ]
    then
        incomplete_flag=`LEU $residue_file`
    elif [ $residue_name == "ILE" ]
    then
        incomplete_flag=`ILE $residue_file`
    elif [ $residue_name == "MET" ]
    then
        incomplete_flag=`MET $residue_file`
    elif [ $residue_name == "PHE" ]
    then
        incomplete_flag=`PHE $residue_file`
    elif [ $residue_name == "TYR" ]
    then
        incomplete_flag=`TYR $residue_file`
    elif [ $residue_name == "TRP" ]
    then
        incomplete_flag=`TRP $residue_file`
    elif [ $residue_name == "PRO" ]
    then
        incomplete_flag=`PRO $residue_file`
    elif [ $residue_name == "SER" ]
    then
        incomplete_flag=`SER $residue_file`
    elif [ $residue_name == "CYS" ]
    then
        incomplete_flag=`CYS $residue_file`
    elif [ $residue_name == "THR" ]
    then
        incomplete_flag=`THR $residue_file`
    elif [ $residue_name == "ASN" ]
    then
        incomplete_flag=`ASN $residue_file`
    elif [ $residue_name == "GLN" ]
    then
        incomplete_flag=`GLN $residue_file`
    elif [ $residue_name == "ASP" ]
    then
        incomplete_flag=`ASP $residue_file`
    elif [ $residue_name == "GLU" ]
    then
        incomplete_flag=`GLU $residue_file`
    elif [ $residue_name == "LYS" ]
    then
        incomplete_flag=`LYS $residue_file`
    elif [ $residue_name == "ARG" ]
    then
        incomplete_flag=`ARG $residue_file`
    elif [ $residue_name == "HIS" ]
    then
        incomplete_flag=`HIS $residue_file`
    fi

    #output:
    echo $incomplete_flag
    }

function HYDROGENATOR {
    
    #definitions:
    protein=$1
    output=$2
    
    #remove altloc residues, if any:
    for i in `cat $protein|grep "^ATOM  "`
    do
        #check for altloc residues:
        val=`echo $i|cut -c 17|grep " "|wc -l`
        if [[ $val -ne 1 ]]
        then
            A_locator=`echo $i|cut -c 17|grep "A"|wc -l`
            if [[ $A_locator -eq 1 ]]
            then
                echo $i|sed 's/^\(.\{16\}\)A/\1 /'
            fi
        else
            echo $i
        fi
    done > $working_directory/Step1_NOaltloc.pdb
    
    #add 'END' to 'Step1_NOaltloc.pdb' to make it easier to detect incomplete residues:
    echo "END                                                                             " >> $working_directory/Step1_NOaltloc.pdb
    
    #check if 'residue.pdb' exists. If so, delete it:
    if [ -f $working_directory/residue.pdb ]
    then
        rm $working_directory/residue.pdb
    fi
    
    #check for completeness of residues:
    master_Incomplete_flag=0
    first_residue_flag=0
    
    residue_ID=`cat $working_directory/Step1_NOaltloc.pdb|head -1|awk '{print substr ($0,18,10)}'`
    file_length=`cat $working_directory/Step1_NOaltloc.pdb|wc -l`
    for i in `cat $working_directory/Step1_NOaltloc.pdb`
    do
        current_ID=`echo $i|awk '{print substr ($0,18,10)}'`
        #IF the same residue is being read:
        if [ $current_ID == $residue_ID ]
        then
            #continue adding atoms to complete the current residue:
            echo $i >> $working_directory/residue.pdb
        #IF a different residue has been read:
        else
            #pass completed residue through completness-detecting function:
            aa_name=`cat $working_directory/residue.pdb|head -1|awk '{print substr($0,18,3)}'`
            Incomplete_flag=`COMPLETE $aa_name $working_directory/residue.pdb`
            if [[ $Incomplete_flag -eq 1 ]]
            then
                echo "error: \"$residue_ID\" contains missing atoms"
                master_Incomplete_flag=1
            fi
            #store first residue for later computations. rename and remove all hydrogen atoms:
            if [[ $first_residue_flag == 0 ]]
            then
                first_residue_ID=`cat $working_directory/residue.pdb|awk '{print substr ($0,21,7)}'`
                cat $working_directory/residue.pdb|awk 'substr($0,13,1)!="H" && substr($0,14,1)!="H" {print substr($0,1,20)" A   1 "substr($0,28,53)}' > $working_directory/first_residue.pdb
                first_residue_flag=1
            fi
            #overwrite previous residue:
            echo $i > $working_directory/residue.pdb
            #reset residue_ID:
            residue_ID=$current_ID
        fi
    done
    rm $working_directory/residue.pdb
    
    if [[ $master_Incomplete_flag == 1 ]]
    then
        echo "Format.sh script will now terminate..."
        rm -r $working_directory
        exit
    fi
 
    #strip the protein of all its hydrogens, as well as the 'END' tag:
    awk 'substr($0,1,6)=="ATOM  " && substr($0,13,1)!="H" && substr($0,14,1)!="H" && substr($0,1,3)!="END" {print $0}' $working_directory/Step1_NOaltloc.pdb > $working_directory/Step2_PREreduce.pdb
    
    #reduce the first residue only, rename back to original:
    first_residue_ID=`cat $working_directory/Step2_PREreduce.pdb|head -1|awk '{print substr ($0,18,10)}'`
    
    $installation_directory/third_party_software/reduce $working_directory/first_residue.pdb|awk 'substr($0,1,6)=="ATOM  " && substr($0,13,4)!=" H2 " && substr($0,13,4)!=" H3 " {print substr($0,1,80)}'|sed 's/^\(.\{12\}\) H1 /\1 H  /'|awk '{print $0"                                "}'|cut -c 1-80|awk -v first_residue_ID=$first_residue_ID '{print substr($0,1,17)"|"first_residue_ID"|"substr($0,28,53)}'|sed s/"|"//g > $working_directory/reduced_first_residue.pdb
    
    #perform reduce, add new hydrogen atoms to protein chain, remove first residue:
    $installation_directory/third_party_software/reduce $working_directory/Step2_PREreduce.pdb|awk 'substr($0,1,6)=="ATOM  " && substr($0,13,4)!=" H2 " && substr($0,13,4)!=" H3 " {print substr($0,1,80)}'|sed 's/^\(.\{12\}\) H1 /\1 H  /'|awk '{print $0"                                "}'|cut -c 1-80|awk -v first_residue_ID=$first_residue_ID 'substr($0,18,10)!=first_residue_ID {print $0}' > $working_directory/Step3_POSTreduce.pdb
    
    #combine (first residue) + (rest of protein) to form the output:
    cat $working_directory/reduced_first_residue.pdb $working_directory/Step3_POSTreduce.pdb > $working_directory/Step4_PREatom.pdb
    
    #add atom names to char[77] of every line:
    for i in `cat $working_directory/Step4_PREatom.pdb`
    do
        atom_1=`echo $i|cut -c 13`
        atom_2=`echo $i|cut -c 14`
        if [ $atom_1 == "H" ]
        then
            atom_name="H"
        else
            atom_name=$atom_2
        fi
        echo $i|cut -c 1-77|awk -v atom_name=$atom_name '{print $0atom_name"  "}'
    done > $output
    
    #cleanup
    rm $working_directory/first_residue.pdb
    rm $working_directory/reduced_first_residue.pdb
    rm $working_directory/Step*
    }

#standard usage/error message: 
if [[ $# -ne 6 ]]
then
    echo "usage: Format.sh -i <protein.pdb> -o <output_directory> -d <database>"
    echo "FLAGS:"
    echo "  -i: input entire protein structure."
    echo "  -o: choose an output directory."
    echo "  -d: path to Autoanneal database."
    exit
fi

#definitions:
while getopts ":i::o::d:" opt; do
    case $opt in
    i)
        IP_protein=$OPTARG
        ;;
    o)
        working_directory=$OPTARG
        ;;
    d)
        installation_directory=$OPTARG
        ;;
    \?)
        echo "error: flag -$OPTARG does not exist." >&2
        exit
        ;;
    :)
        echo "error: flag -$OPTARG requires an argument." >&2
        exit
        ;;
    esac
done
IFS=$'\n' 

#ensure that the <output_directory> given can be worked with:
if [ ! -d $working_directory ]
then
    mkdir $working_directory
fi

contents=`ls $working_directory|wc -l`
if [[ $contents -ge 1 ]]
then
    echo "error: <output_directory> is not empty."
    echo "Format.sh script will now terminate..."
    exit
fi

cp $IP_protein $working_directory/input.pdb

#check completeness of input protein. It must show the following attributes:
#  -80 characters per line +'\n' character.
#  -chain name must be present.
#  -NO altloc residues. By default, the 'A' altloc will be taken.
#  -ALL heavy atoms in each residue (hydrogens optional).

line_counter=1
error_flag=0
for i in `grep "^ATOM  " $IP_protein|grep "GLY\|ALA\|VAL\|LEU\|ILE\|MET\|PHE\|TYR\|TRP\|PRO\|SER\|CYS\|THR\|ASN\|GLN\|ASP\|GLU\|LYS\|ARG\|HIS"`
do

    #check if chain name is present:
    val=`echo $i|cut -c 22|grep "[A-Z]"|wc -l`
    if [[ $val -ne 1 ]]
    then
        echo "error: line $line_counter does not contain a chain name"
        error_flag=1
    fi

    #check for altloc residues:
    val=`echo $i|cut -c 17|grep " "|wc -l`
    if [[ $val -ne 1 ]]
    then
        echo "warning: line $line_counter contains an altloc indicator"
        echo "only 'A' altloc atoms will be used for further calculations"
    fi
    
    line_counter=$((line_counter+1))
    
done

if [[ $error_flag == 1 ]]
then
    echo "Format.sh script will now terminate..."
    rm -r $working_directory
    exit
fi

#add hydrogen atoms to each chain individually:
for i in `grep "^ATOM  " $IP_protein|cut -c 22|sort -u`
do
    grep "^ATOM  " $IP_protein|egrep ^.{21}$i > $working_directory/chain.$i.pdb
    HYDROGENATOR $working_directory/chain.$i.pdb $working_directory/output.$i.pdb
done

#combine all output chains:
cat $working_directory/output.* > $working_directory/pre_output.pdb

#post-output formatting:
# -select only "ATOM  " entries:
# -add end-line atom & ensure 80 characters:
# -remove all atom/residue names with "X"

for i in `grep "^ATOM  " $working_directory/pre_output.pdb|grep -v "X"`
do
    pre_element=`echo "$i"|cut -c 13`
    if [ $pre_element == 'H' ]
    then
        element='H'
    else
        element=`echo "$i"|cut -c 14`
    fi
    echo "$i                                "|cut -c 1-77|awk -v element=$element '{print $0 element"  "}'
done > $working_directory/output.pdb

#cleanup:
rm $working_directory/chain.* $working_directory/pre_output.pdb
ls $working_directory|grep "output"|grep -v "output\.pdb"|awk -v working_directory=$working_directory '{print "rm "working_directory"/"$1}'|bash

#exit script:
echo "Format.sh script will now terminate..."
exit











