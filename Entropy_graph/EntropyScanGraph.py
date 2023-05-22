"""
This script calculates the entropy at each position in an alignment. It then 
creates a graph of the average information over 18 bp along the length of the
alignment
"""


from Bio import SeqIO
import math
import matplotlib.pyplot as plt

def generateProfile(alignmentName):
    """
    >Takes a plain text alignment from Muscle (or other)
    >Returns the frequencies of bases at each position
    >Output as a list of dictionaries
    """
    
    # Open file and add each residue from all sequences to an array
    length = None
    SeqNum = 0 #records the number of sequences
    for record in SeqIO.parse(alignmentName, "fasta"):
        SeqNum += 1
        
        # Create a list that is the length of the alignment and has a 
        # dictionary for each base at each entry
        if length == None:
            length = len(record.seq)
            # Make empty array
            FreqArray = []
            bases = {'A':0, 'C':0, 'G':0, 'T':0, '-':0}
            for i in range(length):
                FreqArray.append(bases.copy())
                
        # Count each residue along the sequence
        for index in range(length):
            base = record.seq[index]
            FreqArray[index][base] += 1
    
    # Calculate the frequency of each residue at each position
    for position in FreqArray:
        for base in position:
            position[base] /= float(SeqNum)
    
    return FreqArray
    
def calculateInformation(profile):
    '''
    Takes a frequency array as input
    Calculates entropy and returns the information present in bits at each position
    '''
    information = [] # list for adding the entropy at each position in a profile
    for residue in profile:
        # Calculate entropy at each position using the equation from Schneider et al (1990)
        entropy = 0
        for value in residue.values():
            if value != 0:
                entropy += value * math.log(value,2)
        entropy *= -1
        # Convert to information present. See Schneider et al (1990)
        information.append(2 - entropy)
    return information

def makeGraph(alignmentFilenames):
    '''
    Draws a graph showing information present averaged over 18 bp
    '''
    plt.clf()
    
    plt.figure(figsize=(10, 3.5))
    
    #Add motif lines - the location of the centre of each motif
    plt.axvline(x=15, color='grey', alpha = 0.3, linestyle='dashed') # C1
    plt.axvline(x=144, color='grey', alpha = 0.3, linestyle='dashed') # C2
    plt.axvline(x=297, color='grey', alpha = 0.3, linestyle='dashed') # C3
    plt.axvline(x=390, color='grey', alpha = 0.3, linestyle='dashed') # C4
    plt.axvline(x=723, color='grey', alpha = 0.3, linestyle='dashed') # C5
    plt.axvline(x=825, color='grey', alpha = 0.3, linestyle='dashed') # C6
    plt.axvline(x=870, color='grey', alpha = 0.3, linestyle='dashed') # C7
    plt.axvline(x=948, color='red', alpha = 0.3) # 1
    plt.axvline(x=993, color='red', alpha = 0.3) # 2
    plt.axvline(x=1068, color='grey', alpha = 0.3) # X
    plt.axvline(x=1200, color='grey', alpha = 0.3, linestyle='dashed') # A1
    plt.axvline(x=1338, color='grey', alpha = 0.3, linestyle='dashed') # A2
    plt.axvline(x=1467, color='grey', alpha = 0.3, linestyle='dashed') # A3
    plt.axvline(x=1587, color='grey', alpha = 0.3, linestyle='dashed') # A4
    plt.axvline(x=1764, color='grey', alpha = 0.3, linestyle='dashed') # A5
    plt.axvline(x=1890, color='grey', alpha = 0.3, linestyle='dashed') # A6
    plt.axvline(x=1971, color='grey', alpha = 0.3, linestyle='dashed') # A7
    plt.axvline(x=2046, color='grey', alpha = 0.3, linestyle='dashed') # A8
    plt.axvline(x=2178, color='grey', alpha = 0.3, linestyle='dashed') # A9
    plt.axvline(x=2229, color='grey', alpha = 0.3, linestyle='dashed') # A10

    # Graph line drawing for each alignment
    for alignmentName in alignmentFilenames:
        
        profile = generateProfile(alignmentName)
        information = calculateInformation(profile)
        data = []
        # Averaging information present
        for i in range(len(information)-18):
            data.append(sum(information[i:i+18])/18)
        
        # Add data 
        wells = range(0, len(data))
        if alignmentName == "Pseudomonas_Streptomyces_gb.fas":
            alpha_value = 1
        else:
            alpha_value = 0.4
        
        plt.plot(wells, data, label = alignmentName, alpha = alpha_value)
    
    axes = plt.gca()
    axes.set_xlim([0, len(data)])
    axes.set_ylim([0, 2])

    plt.legend(loc="upper right")
    
    #Add labels
    params = {'mathtext.default': 'regular' }          
    plt.rcParams.update(params)
    axes.set_xlabel('Alignment position')
    axes.set_ylabel('bits')
    plt.yticks([0, 1, 2])
    plt.savefig('Entropy_scan.svg')


alignmentName = ["Pseudomonas_gb.fas", "Streptomyces_gb.fas", "Pseudomonas_Streptomyces_gb.fas"]
makeGraph(alignmentName)
