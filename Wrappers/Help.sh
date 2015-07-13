#!/bin/bash

set -e
set -u
set -o pipefail

#   A script to serve as an interactive help

#       A function to handle the intial
function helpMe() {
    echo "Welcome to the ANGSD-Wrapper interactive help!"
    echo "Please specify which topic you would like help on:"
    echo " Choose from SFS, 2DSFS, ABBA_BABA, ANC_SEQ, Genotypes, Thetas"
    read choice
    export choice
}

#   A function to figure out if user needs more help
function areDone() {
    echo
    echo "Are you done? [y/n]"
    read finished
    case "$finished" in
        "y" | "yes" | "Y" | "Yes" )
            exit 0
            ;;
        "n" | "no" | "N" | "No" )
            continue
            ;;
    esac
}

function variableDefinitions() {
    source Wrappers/variable_definitions
    while true
    do
        echo
        echo "Please select the argument you are interested in"
        echo "type 'fin' without quotes when finished"
        read argument
        case "$argument" in
            "bam" )
                echo "$BAM"
                continue
                ;;
            "out" )
                continue
                ;;
            "indF" )
                continue
                ;;
            "doSaf" )
                continue
                ;;
            "uniqueOnly" )
                continue
                ;;
            "anc" )
                continue
                ;;
            "minMapQ" )
                continue
                ;;
            "minQ" )
                continue
                ;;
            "nInd" )
                continue
                ;;
            "minInd" )
                continue
                ;;
            "baq" )
                continue
                ;;
            "ref" )
                continue
                ;;
            "GL" )
                continue
                ;;
            "doGLF" )
                continue
                ;;
            "P" )
                continue
                ;;
            "doMajorMinor" )
                continue
                ;;
            "doMaf" )
                continue
                ;;
            "doGeno" )
                continue
                ;;
            "rf" )
                continue
                ;;
            "r" )
                continue
                ;;
            "doPost" )
                continue
                ;;
            "all" )
                continue
                ;;
            "fin" )
                break
                ;;
            * )
                continue
                ;;
        esac
    done
}

while true
do
    helpMe
    echo `pwd`
    case "$choice" in
        "SFS" )
            echo "Site Frequency Spectrum"
            echo "The site frequency spectrum is the number of segregating sites"
            echo
            echo "The argument for SFS are: bam, out, indF, doSaf, uniqueOnly, anc, minMapQ, minQ, nInd, minInd, baq, ref, GL, doGLF, P, doMajorMinor, doMaf, doGeno, rf, r, and doPost"
            variableDefinitions
            areDone
            ;;
        "2DSFS" )
            echo "2D Site Frequency Specturm"
            echo "The 2D site frequency spectrum ist he number of segregating sites for two groups"
            areDone
            ;;
        "ABBA_BABA" )
            echo "ABBA_BABA Test"
            variableDefinitions
            areDone
            ;;
        "ANC_SEQ" )
            echo "Extracting the Ancestral Sequence"
            variableDefinitions
            areDone
            ;;
        "Genotypes" )
            echo "Genotype Calling"
            variableDefinitions
            areDone
            ;;
        "Thetas" )
            echo "Estimating Thetas"
            variableDefinitions
            areDone
            ;;
        * )
            continue
            ;;
    esac
done
