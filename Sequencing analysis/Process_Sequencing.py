# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 12:40:21 2023

@author: calcotma
"""

import os
import numpy as np
from Bio import SeqIO, Alphabet, pairwise2
from Bio.SeqUtils import GC
from Bio.SubsMat import MatrixInfo as matlist
from Bio.Seq import Seq

# Codon usage was taken from the Optimizer tool available at http://genomes.urv.es/OPTIMIZER/
# Data was within the page source after codon optimisation
Codons_HEG_PAO1 = {'GCA' :0.078, 'GCC' :1, 'GCG' :0.322, 'GCT' :0.263, 'TGC' :1, 'TGT' :0.04, 'GAC' :1, 'GAT' :0.22, 'GAA' :1, 'GAG' :0.902, 'TTC' :1, 'TTT' :0.017, 'GGA' :0.006, 'GGC' :1, 'GGG' :0.016, 'GGT' :0.379, 'CAC' :1, 'CAT' :0.203, 'ATA' :0.001, 'ATC' :1, 'ATT' :0.055, 'AAA' :0.167, 'AAG' :1, 'TTA' :0.001, 'TTG' :0.014, 'CTA' :0.002, 'CTC' :0.157, 'CTG' :1, 'CTT' :0.012, 'AAC' :1, 'AAT' :0.072, 'CCA' :0.013, 'CCC' :0.159, 'CCG' :1, 'CCT' :0.068, 'CAA' :0.161, 'CAG' :1, 'AGA' :0.001, 'AGG' :0.001, 'CGA' :0.004, 'CGC' :1, 'CGG' :0.029, 'CGT' :0.631, 'AGC' :0.679, 'AGT' :0.045, 'TCA' :0.002, 'TCC' :1, 'TCG' :0.514, 'TCT' :0.037, 'ACA' :0.009, 'ACC' :1, 'ACG' :0.028, 'ACT' :0.076, 'GTA' :0.218, 'GTC' :1, 'GTG' :0.761, 'GTT' :0.291, 'TAC' :1, 'TAT' :0.107}
Codons_RPG_PAO1 = {'GCA' :0.195, 'GCC' :1, 'GCG' :0.326, 'GCT' :0.669, 'TGC' :1, 'TGT' :0.051, 'GAC' :1, 'GAT' :0.397, 'GAA' :1, 'GAG' :0.7, 'TTC' :1, 'TTT' :0.083, 'GGA' :0.014, 'GGC' :1, 'GGG' :0.03, 'GGT' :0.637, 'CAC' :1, 'CAT' :0.258, 'ATA' :0.002, 'ATC' :1, 'ATT' :0.14, 'AAA' :0.374, 'AAG' :1, 'TTA' :0.002, 'TTG' :0.02, 'CTA' :0.007, 'CTC' :0.173, 'CTG' :1, 'CTT' :0.029, 'AAC' :1, 'AAT' :0.142, 'CCA' :0.037, 'CCC' :0.201, 'CCG' :1, 'CCT' :0.177, 'CAA' :0.269, 'CAG' :1, 'AGA' :0.006, 'AGG' :0.003, 'CGA' :0.006, 'CGC' :0.819, 'CGG' :0.026, 'CGT' :1, 'AGC' :0.801, 'AGT' :0.062, 'TCA' :0.021, 'TCC' :1, 'TCG' :0.356, 'TCT' :0.151, 'ACA' :0.014, 'ACC' :1, 'ACG' :0.038, 'ACT' :0.176, 'GTA' :0.341, 'GTC' :1, 'GTG' :0.645, 'GTT' :0.495, 'TAC' :1, 'TAT' :0.117}

# PvdD sequences for alignment
PvdD_C1_A10	= "TGGCAACTGGAGCCGGAAAGCGCGGCCTACCATATTCCGAGTGCCTTGCGCCTACGCGGGCGGCTGGACGTGGATGCCTTGCAACGCAGCTTCGACAGCCTGGTCGCGCGGCATGAAACCTTGCGTACCCGCTTCCGGCTGGAGGGAGGGCGTTCGTACCAGCAGGTACAACCTGCGGTTAGCGTTTCCATCGAGCGGGAACAGTTCGGTGAAGAAGGCCTGATCGAACGGATACAGGCCATCGTTGTGCAGCCATTCGACCTGGAACGGGGGCCGCTGCTGCGGGTGAACCTGTTGCAACTGGCCGAGGACGACCATGTACTGGTGCTGGTCCAGCACCACATCGTGTCCGATGGTTGGTCGATGCAGGTGATGGTCGAGGAACTGGTCCAGCTCTATGCCGCCTATAGCCAAGGGCTCGACGTGGTGTTGCCAGCCCTGCCGATCCAGTACGCGGACTACGCCCTGTGGCAGCGCAGCTGGATGGAGGCGGGGGAAAAGGAGCGCCAGTTGGCGTACTGGACCGGCCTGCTGGGCGGCGAGCAGCCGGTGCTGGAGTTGCCGTTCGATCGGCCGCGTCCGGCCCGGCAGAGCCATCGTGGCGCGCAGTTGGGTTTCGAGCTATCGCGGGAACTGGTCGAGGCCGTGAGAGCCTTGGCCCAGCGTGAAGGCGCCAGTAGTTTCATGCTGTTGCTGGCCTCGTTCCAGGCGCTGTTGTATCGCTACAGCGGGCAGGCGGATATCCGTGTCGGTGTGCCGATCGCCAATCGCAACCGCGTGGAGACCGAGCGGCTGATCGGCTTCTTCGTCAACACCCAGGTGCTCAAGGCCGACCTGGACGGTCGGATGGGCTTCGACGAGCTGCTGGCCCAGGCCCGCCAACGCGCGCTGGAGGCCCAGGCGCACCAGGACCTGCCGTTCGAGCAACTGGTGGAAGCCTTGCAGCCGGAGCGCAATGCCAGCCACAACCCACTGTTCCAGGTGCTGTTCAACCATCAGAGCGAGATACGCTCGGTGACGCCCGAGGTTCAGTTGGAGGACCTGCGTCTGGAAGGCCTGGCCTGGGACGGCCAGACTGCGCAGTTCGACCTGACGCTGGATATTCAGGAAGACGAAAACGGCATCTGGGCCTCCTTCGACTATGCCACCGATCTGTTCGACGCCTCCACCGTGGAACGCCTGGCCGGCCATTGGCGCAACCTGTTGCGCGGCATCGTCGCCAACCCACGACAGCGGCTCGGCGAGTTGCCGCTGCTGGATGCGCCGGAGCGCCGGCAGACCCTCTCCGAATGGAACCCGGCCCAGCGCGAGTGCGCGGTGCAGGGCACCTTGCAGCAGCGTTTCGAGGAGCAGGCGCGGCAACGGCCACAGGCGGTTGCGCTGATCCTCGACGAACAACGGTTGAGCTACGGCGAACTGAATGCGCGGGCCAATCGCCTGGCGCACTGCCTGATCGCTCGCGGCGTTGGCGCGGACGTGCCGGTCGGGCTGGCGCTGGAGCGTTCGCTGGACATGCTGGTCGGCTTGCTGGCGATCCTCAAGGCCGGCGGCGCCTACCTGCCGTTGGACCCGGCGGCGCCAGAGGAGCGCCTGGCGCATATCCTCGACGACAGTGGGGTACGGCTGCTGCTGACCCAGGGGCATCTGCTCGAGCGCCTGCCGCGGCAGGCGGGGGTGGAGGTGCTGGCCATCGACGGACTGGTGCTGGACGGCTACGCCGAGAGCGATCCGCTCCCGACGCTATCGGCGGACAACCTGGCCTACGTGATCTATACCTCGGGCTCGACCGGCAAGCCCAAGGGCACGTTGCTCACCCACCGCAACGCGCTGCGCCTGTTCAGCGCCACCGAGGCCTGGTTCGGCTTCGACGAGCGGGACGTGTGGACGTTGTTCCATTCCTACGCCTTCGATTTCTCGGTCTGGGAAATCTTCGGCGCGCTGCTCTATGGCGGGCGCCTGGTGATCGTGCCGCAATGGGTGAGCCGTTCGCCGGAAGACTTCTACCGTCTGCTGTGCCGCGAAGGCGTGACGGTGCTCAACCAGACGCCGTCGGCGTTCAAGCAACTGATGGCGGTGGCCTGTTCCGCCGACATGGCGACGCAGCAGCCGGCGCTGCGCTACGTGATCTTCGGTGGCGAGGCGCTGGATCTGCAGAGCCTGCGGCCGTGGTTCCAGCGCTTTGGCGATCGCCAGCCGCAACTGGTGAACATGTACGGCATCACCGAGACCACGGTACACGTAACCTACCGTCCGGTGAGCGAAGCCGACCTGAAGGGTGGCCTGGTCAGTCCGATCGGCGGGACCATCCCGGACCTGTCCTGGTACATCCTCGACCGTGACCTGAACCCGGTGCCGCGCGGCGCGGTGGGCGAGCTGTACATCGGTCGCGCCGGTCTGGCGCGCGGCTACCTGAGGCGGCCCGGGTTGAGTGCCACCCGCTTCGTGCCGAACCCGTTCCCCGGCGGTGCCGGCGAGCGGCTGTACCGTACCGGCGACCTGGCACGGTTCCAGGCGGATGGCAATATCGAGTACATCGGGCGTATCGACCACCAGGTGAAGGTTCGCGGCTTCCGTATCGAACTGGGTGAGATCGAAGCGGCGCTCGCCGGTCTGGCCGGGGTACGCGATGCCGTGGTGCTGGCCCATGACGGGGTCGGCGGCACGCAACTGGTGGGATACGTGGTGGCGGACTCGGCGGAGGATGCCGAGCGTCTGCGGGAGTCGCTGCGGGAGTCGCTGAAGCGGCACCTGCCGGACTACATGGTGCCGGCGCACCTGATGCTGCTGGAGCGGATGCCGCTGACGGTC"
PvdD_C4_A10	= "AGCTGGATGGAGGCGGGGGAAAAGGAGCGCCAGTTGGCGTACTGGACCGGCCTGCTGGGCGGCGAGCAGCCGGTGCTGGAGTTGCCGTTCGATCGGCCGCGTCCGGCCCGGCAGAGCCATCGTGGCGCGCAGTTGGGTTTCGAGCTATCGCGGGAACTGGTCGAGGCCGTGAGAGCCTTGGCCCAGCGTGAAGGCGCCAGTAGTTTCATGCTGTTGCTGGCCTCGTTCCAGGCGCTGTTGTATCGCTACAGCGGGCAGGCGGATATCCGTGTCGGTGTGCCGATCGCCAATCGCAACCGCGTGGAGACCGAGCGGCTGATCGGCTTCTTCGTCAACACCCAGGTGCTCAAGGCCGACCTGGACGGTCGGATGGGCTTCGACGAGCTGCTGGCCCAGGCCCGCCAACGCGCGCTGGAGGCCCAGGCGCACCAGGACCTGCCGTTCGAGCAACTGGTGGAAGCCTTGCAGCCGGAGCGCAATGCCAGCCACAACCCACTGTTCCAGGTGCTGTTCAACCATCAGAGCGAGATACGCTCGGTGACGCCCGAGGTTCAGTTGGAGGACCTGCGTCTGGAAGGCCTGGCCTGGGACGGCCAGACTGCGCAGTTCGACCTGACGCTGGATATTCAGGAAGACGAAAACGGCATCTGGGCCTCCTTCGACTATGCCACCGATCTGTTCGACGCCTCCACCGTGGAACGCCTGGCCGGCCATTGGCGCAACCTGTTGCGCGGCATCGTCGCCAACCCACGACAGCGGCTCGGCGAGTTGCCGCTGCTGGATGCGCCGGAGCGCCGGCAGACCCTCTCCGAATGGAACCCGGCCCAGCGCGAGTGCGCGGTGCAGGGCACCTTGCAGCAGCGTTTCGAGGAGCAGGCGCGGCAACGGCCACAGGCGGTTGCGCTGATCCTCGACGAACAACGGTTGAGCTACGGCGAACTGAATGCGCGGGCCAATCGCCTGGCGCACTGCCTGATCGCTCGCGGCGTTGGCGCGGACGTGCCGGTCGGGCTGGCGCTGGAGCGTTCGCTGGACATGCTGGTCGGCTTGCTGGCGATCCTCAAGGCCGGCGGCGCCTACCTGCCGTTGGACCCGGCGGCGCCAGAGGAGCGCCTGGCGCATATCCTCGACGACAGTGGGGTACGGCTGCTGCTGACCCAGGGGCATCTGCTCGAGCGCCTGCCGCGGCAGGCGGGGGTGGAGGTGCTGGCCATCGACGGACTGGTGCTGGACGGCTACGCCGAGAGCGATCCGCTCCCGACGCTATCGGCGGACAACCTGGCCTACGTGATCTATACCTCGGGCTCGACCGGCAAGCCCAAGGGCACGTTGCTCACCCACCGCAACGCGCTGCGCCTGTTCAGCGCCACCGAGGCCTGGTTCGGCTTCGACGAGCGGGACGTGTGGACGTTGTTCCATTCCTACGCCTTCGATTTCTCGGTCTGGGAAATCTTCGGCGCGCTGCTCTATGGCGGGCGCCTGGTGATCGTGCCGCAATGGGTGAGCCGTTCGCCGGAAGACTTCTACCGTCTGCTGTGCCGCGAAGGCGTGACGGTGCTCAACCAGACGCCGTCGGCGTTCAAGCAACTGATGGCGGTGGCCTGTTCCGCCGACATGGCGACGCAGCAGCCGGCGCTGCGCTACGTGATCTTCGGTGGCGAGGCGCTGGATCTGCAGAGCCTGCGGCCGTGGTTCCAGCGCTTTGGCGATCGCCAGCCGCAACTGGTGAACATGTACGGCATCACCGAGACCACGGTACACGTAACCTACCGTCCGGTGAGCGAAGCCGACCTGAAGGGTGGCCTGGTCAGTCCGATCGGCGGGACCATCCCGGACCTGTCCTGGTACATCCTCGACCGTGACCTGAACCCGGTGCCGCGCGGCGCGGTGGGCGAGCTGTACATCGGTCGCGCCGGTCTGGCGCGCGGCTACCTGAGGCGGCCCGGGTTGAGTGCCACCCGCTTCGTGCCGAACCCGTTCCCCGGCGGTGCCGGCGAGCGGCTGTACCGTACCGGCGACCTGGCACGGTTCCAGGCGGATGGCAATATCGAGTACATCGGGCGTATCGACCACCAGGTGAAGGTTCGCGGCTTCCGTATCGAACTGGGTGAGATCGAAGCGGCGCTCGCCGGTCTGGCCGGGGTACGCGATGCCGTGGTGCTGGCCCATGACGGGGTCGGCGGCACGCAACTGGTGGGATACGTGGTGGCGGACTCGGCGGAGGATGCCGAGCGTCTGCGGGAGTCGCTGCGGGAGTCGCTGAAGCGGCACCTGCCGGACTACATGGTGCCGGCGCACCTGATGCTGCTGGAGCGGATGCCGCTGACGGTC"


def trim_expected_error(seq_record): 
    """Trims the sequence to have the longest sequence with a given probability of less than errors""" 
    
    # Calculate a list of probabilities of error from Q scores    
    qualList = [(10 ** (qual / -10.0)) for qual in 
            seq_record.letter_annotations['phred_quality']]
    
    # Calculate the expected error of an increasingly larger sequence
    end = 1
    probError = qualList[end]
    while end < len(qualList)-1:
        # Calculate error if sequence is extended by 1 bp
        probError += qualList[end+1] 
        # Exit if the error is equal to or greater than a specific threshold
        if probError >=  0.5:
            break           
        end+=1

    return seq_record[:end]

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
    DEGEN = {
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


def sequenceBound(sequence, degPrimer, mismatches = 0, startPos = 0):
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
    return siteFound, bindSite + len(degPrimer)

def calculate_CAI(sequence, Codon_dict):
    '''
    Returns the CAI for a sequence
       
    Uses the formula from:
        Sharp, P. & Li, W. (1987). The codon adaptation index - a measure of 
        directional synonymous codon usage bias, and its potential 
        applications. Nucleic Acids Research, 15(3): 1281-1295
    '''
    data = []
    for i in range(len(sequence)/3):
        codon = sequence[i*3:i*3+3]
        if codon == 'ATG' or codon == 'TGG':
            continue
        data.append(Codon_dict[codon])
        
    data_array = np.array(data)
    return data_array.prod()**(1.0/len(data_array))

def processSequences(targetDir, output_folder):
    '''
    Reads all the sequencing AB1 files from a target directory
    Performs QC and saves the fasta files that pass in a clustering folder
    '''
    dict_sequence_data = {} # A dictionary to record GC content, CAI, etc
    
    fileList = os.listdir(targetDir) # Get the list of files from the target directory
    
    # Open files for output
    forward_fasta = open(output_folder + 'processed_forward_sequences.fasta', 'w')
    reverse_fasta = open(output_folder + 'processed_reverse_sequences.fasta', 'w')
    sequence_data = open('Sequence_data.csv', 'w')
    
    # Variables to print output
    no_primer = 0
    short_insert = 0
    stop_codon = 0
    short_sequence = 0
    
    for fileName in fileList:
        if fileName[-4:] == '.ab1':
            # Open ab1 file
            file_handle = open(os.path.join(targetDir, fileName),'rb')
            seq_record = SeqIO.read(file_handle, 'abi')
            seq_record.seq.alphabet = Alphabet.IUPAC.ambiguous_dna           
            
            # Select the sequences for the vector, primer and short insert check
            if str(fileName.split('_')[-1]) == 'rev.ab1':
                primer_end = IUPACtoDegenerate("BVCGGTCSABCTTGCCGTT") # CA and C4 primer - 22R
                if fileName[1] == '4':
                    short_insert_check = IUPACtoDegenerate("CAACACCACGTCGAGCCC") # Vector on the other site of the insert
                elif fileName[1] == 'A':
                    short_insert_check = IUPACtoDegenerate("CTGCCGATCGGCAAGCGG") # Vector on the other site of the insert
            else:
                short_insert_check = IUPACtoDegenerate("TCGCAGCAGGCCTATCGA") # Vector on the other site of the insert
                if fileName[1] == '4':
                    primer_end = IUPACtoDegenerate("AYTWYRCSSTSTGGCAGCGN") # C4 primer - 10F
                    
                elif fileName[1] == 'A':
                    primer_end = IUPACtoDegenerate("AGSRVCGSMWGTGGTTCCTN") # CA primer - 22F
            
            # Find relevant sequences
            primer_sequence_found = sequenceBound(seq_record.seq, primer_end, mismatches = 1)
            vector_found = sequenceBound(seq_record.seq, short_insert_check, mismatches = 1)
            
            # Discard sequences missing a primer binding site
            if primer_sequence_found[0] == False:
                no_primer += 1
                
            # Discard sequences that contain short inserts
            elif vector_found[0] == True:
                short_insert += 1
            
            else:
                # Trim sequence               
                seq_record = trim_expected_error(seq_record[primer_sequence_found[1]:])
                
                # Trim until only one reading frame
                seq_record = seq_record[:len(seq_record)/3*3]
                
                read_direction = fileName.split('_')[-1][:3] # Fwd vs rev read
                if read_direction == 'rev':
                    seq_record = seq_record.reverse_complement()
                
                # Check for stop codons
                if "*" in seq_record.translate().seq:
                    stop_codon += 1
                
                # Check sequence is at least 100 bp
                elif len(seq_record) < 100:
                    short_sequence += 1
                
                else:
                    ### Calculate sequence parameters and save
                    # Is it a hit?
                    if fileName[0] == "H":
                        plate = fileName.split('_')[1]
                    else:
                        plate = str(0)
                    
                    # The type of substitution boundaries
                    if fileName[1] == "A":
                        substitution = "C1-A10"
                        PvdD_Seq = PvdD_C1_A10
                    else:
                        substitution = "C4-A10"
                        PvdD_Seq = PvdD_C4_A10
                    
                    well = fileName.split('_')[-2] # The well it is located in
                    library = fileName.split('_')[0][2:4] # The library
                    sublibrary = fileName.split('_')[0][4] # The sublibrary
                    
                    unique_ID = '_'.join([fileName[0], substitution, library, sublibrary, plate, well])

                    # Save files that pass in a fasta file for the forward or reverse sequences
                    sequence_name = '>' + fileName.split('_')[0] + ';' + plate + '_' + well + '_' + read_direction
                    
                    if read_direction == 'rev':
                        reverse_fasta.write(sequence_name + '\n')
                        reverse_fasta.write(str(seq_record.seq) + '\n')
                        PvdD_Seq = Seq(PvdD_Seq[len(PvdD_Seq) - len(seq_record.seq):], Alphabet.generic_dna)

                        
                    elif read_direction == 'fwd':
                        forward_fasta.write(sequence_name + '\n')
                        forward_fasta.write(str(seq_record.seq) + '\n')
                        PvdD_Seq = Seq(PvdD_Seq[:len(seq_record.seq)], Alphabet.generic_dna)

                    # Calculate CAI for highly expressed and all ribosomal genes
                    HEG = calculate_CAI(seq_record.seq, Codons_HEG_PAO1)
                    RPG = calculate_CAI(seq_record.seq, Codons_RPG_PAO1)
                    
                    # Calculate GC content
                    GC_content = GC(seq_record.seq)

                    # Run a Needleman-Wunsch algorithm alignment between DNA sequence of PvdD and the sequence
                    # Select the highest percentage identify
                    alignment = pairwise2.align.globalms(PvdD_Seq, seq_record.seq, 5, -4, -10, -.5, 
                                                         one_alignment_only = True)[0]
                    
                    PID = 0.0
                    for i in range(len(alignment[0])):
                        if alignment[0][i] == alignment[1][i]:
                            PID += 1 
                    
                    PID = PID/alignment[-1] * 100
                    
                    # Do an alignment with translated protein                   
                    matrix = matlist.blosum62
                    alignment = pairwise2.align.globalds(PvdD_Seq.translate(), seq_record.translate().seq, matrix, -10, -.5, one_alignment_only = True)[0]
                    aa_PID = 0.0
                    for i in range(len(alignment[0])):
                        if alignment[0][i] == alignment[1][i]:
                            aa_PID += 1 
                    
                    aa_PID = aa_PID/alignment[-1] * 100

                    # Record HEG, RPG, GC, and PID
                    entry = dict_sequence_data.get(unique_ID, {})
                    entry[read_direction] = [str(x) for x in [HEG, RPG, GC_content, PID, aa_PID]]
                    dict_sequence_data[unique_ID] = entry


    # Add HEG, RPG, GC, and PID to an output file
    sequence_data.write(','.join(['Unique_ID', 'Substitution', 'Library', 'Sublibrary', 'Plate', 'Well', 'fwd_HEG', 'fwd_RPG', 'fwd_GC', 'fwd_PID', 'fwd_aaPID', 'rev_HEG', 'rev_RPG', 'rev_GC', 'rev_PID', 'rev_aaPID']))
    sequence_data.write('\n')
    
    for key in dict_sequence_data:
        output = []
        info = key.split('_')
        output.extend([key, info[1], info[2], info[3], info[4], info[5]])

        fwd_data = dict_sequence_data[key].get('fwd', None)
        if fwd_data == None:
            output.extend(['NA', 'NA', 'NA', 'NA', 'NA'])
        else:
            output.extend(fwd_data)
            
        rev_data = dict_sequence_data[key].get('rev', None)
        if rev_data == None:
            output.extend(['NA', 'NA', 'NA', 'NA', 'NA'])
        else:
            output.extend(rev_data)
        
        sequence_data.write(','.join(output))
        sequence_data.write('\n')


    print "There are " + str(len(fileList)) + " files in the target directory"
    print str(no_primer) + " sequences were missing a primer binding site"
    print str(short_insert) + " sequences contained a short insert"
    print str(stop_codon) + " sequences contained a stop codon"
    print str(short_sequence) + " sequences were below 100 bp"

    forward_fasta.close()
    reverse_fasta.close()
    sequence_data.close()

raw_sequence_folder = os.getcwd() + "/Raw_sequences/" 
output_folder = os.getcwd() + '/Clustering/' #desired output folder for sequence fasta files

processSequences(raw_sequence_folder, output_folder)  