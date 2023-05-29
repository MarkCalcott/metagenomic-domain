# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 10:35:46 2022

@author: Mark Calcott
"""
import os
import numpy as np

PYOMZ = {'Gly':1289.58778, 'Ala':1303.60343, 'Ser':1319.59835, 'Pro':1329.61908, 'Val':1331.63473, 'Thr':1333.614, 'Cys':1335.57551, 'Xle':1345.65038, 'Asn':1346.60925, 'Asp':1347.59326, 'Lys':1360.66128, 'Gln':1360.6249, 'Glu':1361.60891, 'Met':1363.60681, 'His':1369.62523, 'Phe':1379.63473, 'Arg':1388.66743, 'Tyr':1395.62965, 'Trp':1418.64563, 'fhOrn':1390.635}

def readPlateMaps(filename):
    '''
    Reads in a platemap file
    Returns a dictionary to convert well locations to samples and the opposite
    '''
    fh = open(filename, 'r')
    
    dict_platemap_hit_to_sample = {}
    dict_platemap_sample_to_hit = {}

    while True:
        # Read the first line and get the plate name
        plate = fh.readline().split(",")[0]
        
        dict_plate = {}
        # Go through each row of the plate
        for i in range(8):
            row = fh.readline().strip().split(",")
            # Go through each well
            for j in range(12):
                # Record the wells name
                well = row[0] + str("%02d" % (j+1,))
                # If the platemap records a sample for the well, add it to dict_platemap_hit_to_sample
                if row[j+1] != "NA":
                    dict_plate[well] = row[j+1]
                    # Add to dict_platemap_sample_to_hit
                    if len(plate) == 4: # All letter variants of day 5 are repeated in other plates
                        dict_platemap_sample_to_hit[row[j+1]] = (plate, well)
        
        dict_platemap_hit_to_sample[plate] = dict_plate
                
        if fh.readline() == "":
            break
        fh.readline()
    return dict_platemap_sample_to_hit, dict_platemap_hit_to_sample

def extractAbsEm():
    '''
    Read all the plates in a folder
    Return well location, fluorescence and plate name
    '''
    # Save folder name containing the files
    
    complete_Abs = {}
    complete_Em = {}
    for folder in ["/AbsEm/Day1", "/AbsEm/Day2", "/AbsEm/Day3", "/AbsEm/Day4", "/AbsEm/Day5", "/AbsEm/Day6"]:
        raw_Abs = {}
        raw_Em = {}
        
        filePath = os.getcwd() + folder
        
        # Identify all files in the folder
        files = []
        for file in os.listdir(filePath):
            files.append(file)

        # Read data
        emControls = {"Positive":[], "Negative":[]}
        absControls = {"Positive":[], "Negative":[]}
        for plate in files:
            plateName = plate.split('_')[0]
    
            # Find the first well   
            handle = open(filePath + "/" + plate, 'r')
            
            while True:
                words = handle.next().split(',')
                
                if len(words) > 0 and words[0] == "Well":
                    break
            
            # Read in absorbance and fluorescence data
            for i in range(96):
                
                data = handle.next().split(',')[:3]

                for i in range(3):
                    handle.next()
                
                # Check fluorescence is greater than 1000
                # Reason: Some wells had no growth and therefore no fluorescence. All deletion strains had at least 1277
                if float(data[1]) > 1000:
                    # Record absorbance
                    key = (plateName, data[0])
                    entry = raw_Abs.get(key, [])
                    entry.append(float(data[2]))
                    raw_Abs[key] = entry           
                    
                    # Record fluorescence emission
                    key = (plateName, data[0])
                    entry = raw_Em.get(key, [])
                    entry.append(float(data[1]))
                    raw_Em[key] = entry    
                
                # The controls were in different locations, so this calcuates the controls data
                if plateName[:4] == "D2P5":
                    if data[0] in ["D01", "D03", "D05"]:
                        # Record the Absorbance positive controls
                        entry = absControls["Positive"]
                        entry.append(float(data[2]))
                        absControls["Positive"] = entry
                        # Record the emission positive controls
                        entry = emControls["Positive"]
                        entry.append(float(data[1]))
                        emControls["Positive"] = entry
                    elif data[0] in ["D02", "D04", "D06"]:
                        # Record the Absorbance positive controls
                        entry = absControls["Negative"]
                        entry.append(float(data[2]))
                        absControls["Negative"] = entry
                        # Record the emission positive controls
                        entry = emControls["Negative"]
                        entry.append(float(data[1]))
                        emControls["Negative"] = entry
        
        # Normalise emission data to controls
        absNegative = sum(absControls["Negative"])/3
        absPositive = sum(absControls["Positive"])/3 - absNegative
        emNegative = sum(emControls["Negative"])/3
        emPositive = sum(emControls["Positive"])/3 - emNegative
        
        for key in raw_Em.keys():
            # Check well contains a sample
            if key[1] in dict_platemap_hit_to_sample[key[0]].keys():
                # Get sample name
                sample = dict_platemap_hit_to_sample[key[0]][key[1]]
                # Add absorbance data
                entry = complete_Abs.get(sample, [])
                
                entry.extend([(x-absNegative)/absPositive*100 for x in raw_Abs[key]])
                complete_Abs[sample] = entry
                # Add fluorescence data
                entry = complete_Em.get(sample, [])
                entry.extend([(x-emNegative)/emPositive*100 for x in raw_Em[key]])
                complete_Em[sample] = entry


    return complete_Abs, complete_Em


def extractMALDI():
    '''
    Reads the CSV files in the MALDI folder and add data to each sample
    '''
    
    # Identify all files in the MALDI folder
    filePath = os.getcwd() + "/MALDI"
    files = []
    for file in os.listdir(filePath):
        files.append(file)
    
    # Read data
    data = {}
    for plate in files:
        plateName = plate.split('.')[0]
        plates = plateName.split('_')[-2:]
        handle = open(filePath + "/" + plate, 'r')
        
        ### Read in the first 5 lines containing calibration spots
        cal_masses = [757.4, 1046.542, 1533.858]
        errors = []
        for i in range(5):
            row = handle.readline().split(",")
            for j in range(len(cal_masses)):
                read = float(row[j*2+1].strip('"').split("(")[0])
                errors.append(read - cal_masses[j])
        error = sum(errors)/len(errors)

        # Read in the value for each well
        for replicate in range(2):
            for column in "ABCDEFGH":
                for row in range(12):
                    for plate in plates:
                        reading = handle.next().split(',')[1].strip('"')
                        if reading != "N":
                            value = float(reading.split("(")[0]) - error
                            well = column + str("%02d" % (row + 1,))

                            sample_name = dict_platemap_hit_to_sample[plate].get(well, None)
                            if sample_name != None:
                                entry = data.get(sample_name, [])
                                entry.append(value)
                                data[sample_name] = entry

    return data

def labelMALDI(mass):
    hit = "NA"
    for pyo in PYOMZ.keys():
        if abs(mass - PYOMZ[pyo]) < 0.5:
            hit = pyo
    return hit

# Create the platemaps to get from sample to hit and vice versa
dict_platemap_sample_to_hit, dict_platemap_hit_to_sample = readPlateMaps("PlateMaps.csv")

# Extract the MALDI data
MALDI_data = extractMALDI()

# Extract the absorbance and emission data
Abs_data, Em_data = extractAbsEm()

# Write the data into a files
fh = open("Screen_data.csv", "w")
fh_MALDI = open("Screen_MALDI_data.csv", "w")
fh_ABS = open("Screen_ABS_data.csv", "w")

fh.write('Unique_ID,Substitution,Library,Sublibrary,ScreeningPlate,ScreeningWell,MALDI,Msamples,Substrate,Abs,Aerror,Asamples,Fluoro,Ferror,Fsamples\n')
fh_MALDI.write('Sample,ScreeningPlate,ScreeningWell,Substrate,ExpectedMALDI,AverageMALDI,Samples,Data\n')
fh_ABS.write('Sample,ScreeningPlate,ScreeningWell,Average,StandardDeviation,Samples,Data\n')
               
for key in dict_platemap_sample_to_hit.keys():
    # Add the sample, motif, library, metagenomic library plate, screening plate and well data
    output = []
    if key == 'WT' or key == 'MUT':
        output.extend([key, 'NA', 'NA', 'NA'])
        
    else:
        Substitution = key[:6]
        Library = key.split('-')[2][:2]
        Sublibrary = key.split('-')[2][-1]
        Plate = key.split('-')[3]
        unique_ID = '_'.join(['H', Substitution, Library, Sublibrary, Plate])
        output.extend([unique_ID,Substitution,Library,Sublibrary])
        
    output.extend(dict_platemap_sample_to_hit[key])

    # Add average MALDI data
    MALDI = MALDI_data.get(key, None)
    if MALDI == None:
        output.extend(["NA", 0, "NA"])
    else:
        output.extend([np.mean(MALDI), len(MALDI), labelMALDI(np.mean(MALDI))])
        
    # Add average absorbance
    output.extend([np.mean(Abs_data[key]), np.std(Abs_data[key]), len(Abs_data[key])])
    
    # Add average emission
    output.extend([np.mean(Em_data[key]), np.std(Em_data[key]), len(Em_data[key])])
    fh.write(','.join(str(x) for x in output))
    fh.write('\n')
    
    # Write the MALDI data
    MALDI = MALDI_data.get(key, None)
    MALDI_output = []
    MALDI_output.append(key)
    MALDI_output.extend(dict_platemap_sample_to_hit[key])
    
    if MALDI == None:
        MALDI_output.extend((["NA", "NA", 0, "NA", "NA", "NA", "NA", "NA"]))
    else:
        substrate = labelMALDI(np.mean(MALDI))
        MALDI_output.append(substrate)

        if substrate == "NA":
            MALDI_output.append("NA")
        else:
            MALDI_output.append(PYOMZ[substrate])
        
        MALDI_output.extend([np.mean(MALDI), len(MALDI)])
        MALDI_output.extend(MALDI)
    
    fh_MALDI.write(','.join(str(x) for x in MALDI_output))
    fh_MALDI.write('\n')

    # Write absorbance data
    if dict_platemap_sample_to_hit[key][0].startswith("D4"):
        ABS_output = []
        ABS_output.append(key)
        ABS_output.extend(dict_platemap_sample_to_hit[key])
        
        # Add average absorbance
        ABS_output.extend([np.mean(Abs_data[key]), np.std(Abs_data[key]), len(Abs_data[key])])
        
        # Add each value
        ABS_output.extend(Abs_data[key])
        fh_ABS.write(','.join(str(x) for x in ABS_output))
        fh_ABS.write('\n')

fh.close()
fh_MALDI.close()
fh_ABS.close()


