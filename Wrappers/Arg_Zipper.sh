#!/bin/bash

set -e
set -o pipefail

if [ "$#" -lt 2 ]
then
    echo -e "\
    Usage: Arg_Zipper <wrapper args> <advanced args> \n\
    where: <wrapper args> and <advanced args> are ideally both lists of flag-value\n\
           paired strings." >&2
    exit 1
fi

echo -e "\
---------------------------------------------------
WRAPPER: Zipping advanced arguments onto basic ones\n\
---------------------------------------------------" 1>&2
wrapper_args=$1 # Which arguments are set by default?
adv_args=$2 # Which arguments were defined by an advanced user?
# Debugging outputs
# echo -e "Wrapper's: ${wrapper_args}\n" 1>&2
# echo -e "User's: ${adv_args}\n" 1>&2

# Initialize lists for flags and values, and the final aggregate string
FINAL_FLAGS=()
FINAL_VALS=()
FINAL_ARGS=""

# Add the initial wrapper arguments to the lists
k=1
NUM_ARGS=$(echo "${wrapper_args}" | wc -w)
while (( k <= NUM_ARGS ))
do
    flag=$(echo "${wrapper_args}" | cut -d " " -f $(( k )))
    val=$(echo "${wrapper_args}" | cut -d " " -f $(( k+1 )))
    # DEBUGGING: echo "Flag and value pair: $flag '$val' $k" >&2


    # Validate each pair (avoid whitespace errors or missing arguments).
    # If the pair is invalid, discard the only first item and proceed normally.
    if [[ $val = -* || $val = '' ]]
    then
        echo "Found a flag '${flag}' with no associated value, discarding." 1>&2
        let k=k+1
        continue
    fi
    if [[ $flag != -* ]]
    then
        echo "Found a value '${val}' with no associated flag, discarding flag '${flag}'." 1>&2
        let k=k+1
        continue
    fi

    FINAL_FLAGS+=("${flag}")
    FINAL_VALS+=("${val}")
    let k=k+2
done

# Proceed through each of the advanced or user-defined arguments
k=1
NUM_ARGS=$(echo "${adv_args}" | wc -w)
while (( k <= NUM_ARGS ))
do
    flag=$(echo "${adv_args}" | cut -d " " -f $(( k )))
    val=$(echo "${adv_args}" | cut -d " " -f $(( k+1 )))
    # DEBUGGING: echo "Flag and value pair: $flag $val $k" >&2

    # Validate each pair (avoid whitespace errors or missing arguments).
    # If the pair is invalid, discard the only first item and proceed normally.
    if [[ $val = -* || $val = '' ]]
    then
        echo "Found a flag '${flag}' with no associated value, discarding." 1>&2
        let k=k+1
        continue
    fi
    if [[ $flag != -* ]]
    then
        echo "Found a value '${val}' with no associated flag, discarding flag '${flag}'." 1>&2
        let k=k+1
        continue
    fi

    # Check for pre-existing flags in the wrapper's args with the same name
    # If such a flag exists, overwrite the associated value.
    # Otherwise append the pair to the lists
    found=false
    for i in "${!FINAL_FLAGS[@]}"
    do
        if [[ "${FINAL_FLAGS[i]}" = "${flag}" ]]
        then
            echo "Found an overlapping flag '${flag}', overwriting value '${FINAL_VALS[i]}' with '${val}'" 1>&2
            FINAL_VALS[i]="${val}"
            found=true
            break
        fi
    done

    if [[ $found = false ]]
    then
        FINAL_FLAGS+=("${flag}")
        FINAL_VALS+=("${val}")
    fi

    # Increment to grab the next flag-argument pair
    let k=k+2
done

# Merge the flags and associated values into the final string
for i in "${!FINAL_FLAGS[@]}"
do
    FINAL_ARGS+="${FINAL_FLAGS[i]} ${FINAL_VALS[i]} "
done

echo "${FINAL_ARGS}"
