# -*- coding: utf-8 -*-
"""
Created on Mar 2 2022

@author: Mark Calcott
"""

import os

def readPlates(dataDirectory):
    '''
    Read all the plates in a folder
    Return well location, fluorescence and plate name
    '''
    
    # Save folder name containing the files
    filePath = os.getcwd() + "/" + dataDirectory
    
    # Identify all files in the folder
    files = []
    for file in os.listdir(filePath):
        files.append(file)
    
    complete_data = []
    # Read data
    for plate in files:
        plateName = plate.split('.')[0]

        # Find the first well   
        handle = open(filePath + "/" + plate, 'r')
        
        while True:
            words = handle.next().split(',')
            
            if len(words) > 0 and words[0] == "Well":
                break
        
        for i in range(96):
            data = handle.next().split(',')[:2]
            data.append(plateName)
            complete_data.append(data)
        handle.close()
    return complete_data

def create_hit_platemap(dataDirectory, threshold = 10000):
    '''
    In a folder, read all 96 well plate files and identify hits above a threshold
    Save plate maps containing locations to add samples
    '''
    
    C4_hit = 0
    CA_hit = 0
    # Read data
    complete_data = readPlates(dataDirectory)

    # Find all hits above threshold
    hits = []
    for well in complete_data:
        if int(well[1]) > threshold:
            hits.append(well)

            if well[2][:2] == ("C4"):
                C4_hit += 1
            elif well[2][:2] == ("CA"):
                CA_hit += 1

    # saveFile name
    handle = open(dataDirectory + "_fluoro.txt", 'w')
    
    for hit in hits:
        hitName = hit[2] + "_" + hit[0]
        handle.write(hitName + "," + hit[1] + "\n")

    print str(C4_hit), "C4-A10 hits,", str(CA_hit), "C1-A10 hits,", len(hits), "total hits,"
    handle.close()

for name in ("Day1", "Day2", "Day3"):
    create_hit_platemap(name)