import re


def yieldSequences(sequenceFile):
    '''
    Returns each sequence in a file
    '''
    with open(sequenceFile,'rU') as F:
        data = 'placeHolder'
        sequence = []

        #check sequence file starts with a fasta sequence
        firstLine = F.readline().strip()
        if not firstLine[0] == '>':
            raise ValueError('sequence file does not begin with ">"')

        #yield sequences one at a time
        while data:
            data = F.readline().strip()
            
            if data and not data[0] == '>':
                sequence.append(data)
            else:
                degappedSequence = ''.join(sequence)
                degappedSequence = re.sub('-', '', degappedSequence)
                yield degappedSequence
                sequence = []

def hammingDistance(degenerateKmer, kmer):
    '''
    Takes in a list containing all residues at each position in a degenerate kmer and a kmer
    Calculates the hamming distance between these
    '''
    distance = 0
    for residue in range(len(kmer)):
        if not kmer[residue] in degenerateKmer[residue]:
            distance += 1
    return distance

def IUPACtoDegenerate(IUPACprimer):
    '''
    Takes a IUPAC primer and converts it to a degenerate data structure
    '''
    DEGEN = {#Probably double check this
    'A':set(['A']),
    'C':set(['C']),
    'G':set(['G']),
    'T':set(['T']),
    'R':set(['A','G']),
    'Y':set(['C','T']),
    'S':set(['C','G']),
    'W':set(['A','T']),
    'K':set(['G','T']),
    'M':set(['A','C']),
    'D':set(['A','G','T']),
    'H':set(['A','C','T']),
    'V':set(['A','C','G']),
    'B':set(['C','G','T']),
    'N':set(['A','C','G','T'])
    }

    degeneratePrimer = []
    for char in IUPACprimer:
        degeneratePrimer.append(DEGEN[char].copy())
    return degeneratePrimer

def revCompDegenerate(IUPACprimer):
    '''
    Takes a set for a degenerate primer and reverse complements it
    '''
    complementIUPAC = {
    'A':'T',
    'C':'G',
    'G':'C',
    'T':'A',
    'R':'Y',
    'Y':'R',
    'S':'S',
    'W':'W',
    'K':'M',
    'M':'K',
    'D':'H',
    'H':'D',
    'V':'B',
    'B':'V',
    'N':'N'
    }
    degeneratePrimer = []
    for char in range(len(IUPACprimer)-1, -1, -1):
        residue = IUPACprimer[char]
        degeneratePrimer.append(complementIUPAC[residue])

    return ''.join(degeneratePrimer)

def sequenceBound(sequence, degPrimer, mismatches = 1, startPos = 0):
    '''
    Returns whether a primer is found in a sequence and the first location it is found
    '''
    siteFound = False
    bindSite = 0
    for bindSite in range(startPos, len(sequence)-len(degPrimer)+1):
        bindSiteSequence = sequence[bindSite:bindSite+len(degPrimer)]
        if hammingDistance(degPrimer, bindSiteSequence) <= mismatches:
            siteFound = True
            break
    return siteFound, bindSite+len(degPrimer)


def PCRscreen(sequence_file, degFwdPrimer, degRevPrimer, mismatches):
    '''
    Given a file containing FASTA sequences, finds all sequences amplified with a given number of mismatches
    Assumes forward primer is designed to bind the coding strand
    '''
    sequencesAmplified = 0
    
    count = 0
    amplifiedLengths = []
    for sequence in yieldSequences(sequence_file):
        count += 1
        # Find the first forward primer binding site
        fwdSiteFound, location1 = sequenceBound(sequence, degFwdPrimer, mismatches)
        # Find an upstream reverse primer binding site
        if fwdSiteFound:
            revSiteFound, location2 = sequenceBound(sequence, degRevPrimer, mismatches, location1)
            if revSiteFound:
                sequencesAmplified += 1
                amplifiedLengths.append(location2-location1)

    return round(sequencesAmplified / float(count) * 100, 1)

def calculate_degeneracy(degenerate_primer):
    '''
    Takes a degenerate primer as input and returns the level of degeneracy
    '''
    degeneracy = 1
    for i in degenerate_primer:
        degeneracy *= len(i)
    return degeneracy

def wrapper():
    '''
    Reads in primer combinations one at a time
    Calculates the primer degeneracy
    Determenines how many sequences have binding sites in each of the Pseudomonas, Streptomyces and Bacillus files
    Outputs the primer sequences, degeneracy and number of sequences amplified
    '''
    
    fh = open('Input_degenerate_primers.csv', 'r')
    output_file = open('Degenerate_primer_info_table.csv', 'w')
    output_file.write('Primer set, forward, reverse, forward degeneracy, reverse degeneracy, Pseudomonas amplification percentage, Streptomyces amplification percentage, Bacillus amplification\n')
    for x in fh.readlines():
        primer_info = x.strip().split(',')

        # Convert the primers to IUPAC degenerate format for searching
        degFwdPrimer = IUPACtoDegenerate(primer_info[1])
        degRevPrimer = revCompDegenerate(primer_info[2])
        degRevPrimer = IUPACtoDegenerate(degRevPrimer)
        
        # Add the level of degeneracy to primer_info
        primer_info.append(calculate_degeneracy(degFwdPrimer))
        primer_info.append(calculate_degeneracy(degRevPrimer))
        
        # Search each sequence file for the number of binding sites and add to primer_info
        for sequence_file in ["Pseudomonas_sequences.fas", "Streptomyces_sequences.fas", "Bacillus_sequences.fas" ]:
            primer_info.append(PCRscreen(sequence_file, degFwdPrimer, degRevPrimer, mismatches = 1))
        
        # Write primer info to the output file
        primer_info.append('\n')
        output_file.write(','.join(str(x) for x in primer_info))
    fh.close()
    output_file.close()
wrapper()
