#!/usr/bin/env python3

import sys
if not sys.version_info[0] == 3:
    sys.exit("Please use Python 3 for this program")


import re
import os
import argparse
import itertools
import collections


def create_argument_parser():
    '''Create a parser for the arguments'''
    parser = argparse.ArgumentParser(add_help=True)
    parser.add_argument(
        '--regions',
        type=str,
        default=None,
        dest='regions',
        required=True,
        help="Unsorted regions file"
    )
    parser.add_argument(
        '--fai',
        type=str,
        default=None,
        dest='fai',
        required=True,
        help="Fai index file containing sort order"
    )
    parser.add_argument(
        '--project',
        type=str,
        default="SORTED",
        dest='project',
        required=False,
        help="Name of the project, defaults ot 'SORTED'"
    )
    return(parser)


def get_fai_regions(faifile):
    '''Find the regions within the fai index file'''
    print("Finding regions in " + faifile + "...")
    try:
        fai = open(faifile).read() # Read in the fai file
        fai_regions = re.findall(r'^(.*?)\s', fai, re.MULTILINE) # Take the first column of every line
        print("Found " + str(len(fai_regions)) + " regions in " + faifile)
        return(fai_regions)
    except FileNotFoundError:
        raise


def get_duplicates(regions):
    print("Checking duplicate regions...")
    counter = collections.Counter(regions)
    dups = [r for r, x in counter.items() if x > 1]
    print("Found " + str(len(dups)) + " duplicate regions...")


def check_found_regions(found_regions, regions):
    '''Ensure we've found all of our regions properly'''
    for region in regions:
        try:
            # If we found the region, do nothing
            if region in found_regions:
                pass
        except:
            # Otherwise, let us know that we could not find a region
            print("Failed to find " + region)


def find_regions(fai_regions, regionsfile):
    '''Create a list of all regions shared between the regions file and the fai index'''
    print("Collecting regions from " + regionsfile + "...")
    try:
        regions = open(regionsfile).read().split('\n') # Get the regions from the file
        get_duplicates(regions) # Count the number of duplicates
        region_set = set(regions) # Create a set of the regions
        fai_set = set(fai_regions) # Create a set of the regions from the fai index
        found_regions = {r for r in region_set if r.split(':')[0] in fai_set} # Collect a list of found regions
        found_regions.discard('') # Remove any ''
        check_found_regions(found_regions, regions) # Check to ensure all regions were found
        print("Found " + str(len(found_regions)) + " unique regions")
        return(found_regions)
    except FileNotFoundError:
        raise


def sort_regions(found_regions, fai_regions):
    '''Sort the regions according to the fai index'''
    print("Sorting regions...")
    get_contig = re.compile(r'([0-9])+') # A regex to collect numbers from a contig
    sort_order = lambda contig: int(get_contig.search(contig.split(':')[0]).group()) # Create a lambda expression to handle the sort
    sorted_regions = sorted(found_regions, key=sort_order) # Sort the regions
    return(sorted_regions)


def write_sorted_regions(sorted_regions, project, regionsfile):
    '''Write the sorted regions to an output file'''
    print("Writing the output file...")
    output_directory = os.path.dirname(os.path.abspath(regionsfile)) # Get the directory where our regions file is
    output_name = project + '_SortedRegions.txt' # Create a name for our output file
    output_file = output_directory + '/' + output_name # Append the name to the directory
    out = open(output_file, 'w') # Open the file for writing
    for region in sorted_regions: # Write every region to the file
        out.write(region + '\n')
    out.close() # Close the file
    return(output_file)


def main():
    '''Run the program'''
    parser = create_argument_parser() # Create our argument parser
    if not sys.argv[1:]: # If we don't have any arguments
        parser.print_help() # Run the help message
    args = vars(parser.parse_args()) # Create a dictionary of arguments
    try:
        # Collect regions from the fai index and regions file
        fai_regions = get_fai_regions(args['fai'])
        region_list = find_regions(fai_regions, args['regions'])
    except FileNotFoundError as e: # If either don't exists
        sys.exit("Failed to find " + e.filename) # Exit with error
    sorted_regions = sort_regions(region_list, fai_regions) # Sort our regions
    output_file = write_sorted_regions(sorted_regions, args['project'], args['regions']) # Write our sorted regions to an output file
    print("Sorted regions file can be found at " + output_file)


main() # Run the program
