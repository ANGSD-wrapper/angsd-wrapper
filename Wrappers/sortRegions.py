#!/usr/bin/env python3

#   Import required modules from the standard Python library
import sys
import re
import os
try:
    import argparse
except ImportError:
    sys.exit("Please install argparse or update to a newer version of Python")

Parser = argparse.ArgumentParser(add_help=True)
Parser.add_argument('-r',
    '--regions',
    type=str,
    default=None,
    dest='regions',
    required=True,
    help="Unsorted regions file")
Parser.add_argument('-f',
    '--fai',
    type=str,
    default=None,
    dest='fai',
    required=True,
    help="Fai file for outgroup")
Parser.add_argument('-p',
    '--project',
    type=str,
    default="SORTED",
    dest='project',
    required=False,
    help="Name of the project, defaults to 'SORTED'")


#   Find the regions in the fai file
def findFai(faifile):
    '''Find the regions within the fai file'''
    print("Finding regions in the fai file...")
    fai = open(faifile).read() # Read in the fai file
    faireg = re.findall(r'^(.*?)\s', fai, re.M) # Take the first column of every line
    print("Found " + str(len(faireg)) + " regions in the fai file")
    return(faireg) # Return our list of regions


#   Get the index for every region in the regions file with relation to the fai file
def getIndex(faireg, regionsfile):
    '''Generate a dictionary of all regions with their order in the fai file'''
    print("Sorting the regions file...")
    regDict = dict() # Set up our holding dictionary
    regions = open(regionsfile) # Open the regions file
    for reg in regions: # For every region
        noCol = reg.split(':')[0] # Get rid of the colon
        try: # Try to find the region in the fai file
            index = faireg.index(noCol) # Find the index
            regDict[index] = reg # Store this in the dictionary, with the index as the key
        except ValueError: # If we can't find it
            print("Cannot find " + noCol + " in the fai file!")
            print("Skipping...")
            continue # Skip
    print("Sorted " + str(len(regDict)) + " regions")
    return(regDict) # Return our dictionary


#   Write the outfile in sorted order
def sortRegions(regDict, project, regionsfile):
    '''Sort the regions and write to file'''
    print("Writing the outfile...")
    sortOrder = sorted(regDict.keys()) # Sort the keys
    outfile = os.path.dirname(regionsfile) + '/' + project + '_SortedRegions.txt' # Generate a name
    out = open(outfile, 'w') # Open the outfile for writing
    for index in sortOrder: # For every index that's been sorted
        out.write(regDict[index]) # Write the region
    out.close() # Close the file
    return(outfile) # Return the outname


def main():
    '''Run the program'''
    if not sys.argv[1:]: # If we don't have any arguments
        sys.exit(Parser.print_help()) # Print the help message and exit
    args = vars(Parser.parse_args()) # Generate a dictionary of arguments
    faireg = findFai(args['fai']) # Create our list of regions in the fai file
    regDict = getIndex(faireg, args['regions']) # Create our dictionary of regions with indecies
    outfile = sortRegions(regDict, args['project'], args['regions']) # Write our outfile
    print("Sorted regions file can be found at " + outfile)


main()
