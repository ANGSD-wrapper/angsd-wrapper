#!/bin/bash

set -e
set -u


echo -e "Zipping advanced arguments onto basic ones...\n" 1>&2
wrapper_args=$1 # Which arguments are set by default?
adv_args=$2 # Which arguments were defined by an advanced user?
# Debugging outputs
#     echo -e "Wrapper's: ${wrapper_args}\n" 1>&2
#     echo -e "User's: ${adv_args}\n" 1>&2

# Initialize lists for flags and values, and the final aggregate string
FINAL_FLAGS=()
FINAL_VALS=()
FINAL_ARGS=""

# Need to check existence of first flag-argument pair before entering while loop (~do-while loop)
k=1
flag=$(echo "${wrapper_args}" | cut -d " " -f $(( k )))
val=$(echo "${wrapper_args}" | cut -d " " -f $(( k+1 )))
while [[ ! -z $flag  && ! -z $val ]]
do
    FINAL_FLAGS+=("${flag}")
    FINAL_VALS+=("${val}")
    let k=k+2
    flag=$(echo "${wrapper_args}" | cut -d " " -f $(( k )))
    val=$(echo "${wrapper_args}" | cut -d " " -f $(( k+1 )))
done

# Need to check existence of first flag-argument pair before entering while loop (~do-while loop)
k=1
flag=$(echo "${adv_args}" | cut -d " " -f $(( k )))
val=$(echo "${adv_args}" | cut -d " " -f $(( k+1 )))
while [[ ! -z "${flag}"  && ! -z "${val}" ]]
do
    # Check for pre-existing flags in the wrapper's args that share the same name
    # If such a flag exists, overwrite the associated value.
    # Otherwise append the pair to the list
    found=false
    for i in "${!FINAL_FLAGS[@]}"
    do
	if [[ "${FINAL_FLAGS[i]}" = "${flag}" ]]
	then
	    echo "Found an overlapping flag '${flag}', overwriting value '${FINAL_VALS[i]}' with '${val}'" 1>&2
	    FINAL_VALS[i]="${val}"
	    found=true
	fi
    done

    if [ $found = false ]
    then
	FINAL_FLAGS+=("${flag}")
	FINAL_VALS+=("${val}")
    fi
    (( k=k+2 ))
    flag=$(echo "${adv_args}" | cut -d " " -f $(( k )))
    val=$(echo "${adv_args}" | cut -d " " -f $(( k+1 )))
done

for i in "${!FINAL_FLAGS[@]}"
do
    # echo -e "FINAL_FLAGS[${i}]: ${FINAL_FLAGS[i]}\t" 1<&2
    FINAL_ARGS+="${FINAL_FLAGS[i]} ${FINAL_VALS[i]} "
done

echo "${FINAL_ARGS}"
