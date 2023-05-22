import schema
import matplotlib.pyplot as plt


def makeIndividualBarGraph(complete_data, name):
    '''
    Draws a graph containing the ordered data and OD600 from each well
    Saves based on sample number
    '''
    plt.clf()
    plt.figure(figsize=(10, 3))

    #Add motif lines - the location of the centre of each motif
    plt.axvline(x=5, color='grey', alpha = 0.3, linestyle='dashed') #C1
    plt.axvline(x=51, color='grey', alpha = 0.3, linestyle='dashed') #C2
    plt.axvline(x=127, color='grey', alpha = 0.3, linestyle='dashed') #C3
    plt.axvline(x=164, color='grey', alpha = 0.3, linestyle='dashed') #C4
    plt.axvline(x=281, color='grey', alpha = 0.3, linestyle='dashed') #C5
    plt.axvline(x=315, color='grey', alpha = 0.3, linestyle='dashed') #C6
    plt.axvline(x=332, color='grey', alpha = 0.3, linestyle='dashed') #C7
    plt.axvline(x=481, color='grey', alpha = 0.3, linestyle='dashed') #A1
    plt.axvline(x=528, color='grey', alpha = 0.3, linestyle='dashed') #A2
    plt.axvline(x=603, color='grey', alpha = 0.3, linestyle='dashed') #A3
    plt.axvline(x=648, color='grey', alpha = 0.3, linestyle='dashed') #A4
    plt.axvline(x=750, color='grey', alpha = 0.3, linestyle='dashed') #A5
    plt.axvline(x=807, color='grey', alpha = 0.3, linestyle='dashed') #A6
    plt.axvline(x=840, color='grey', alpha = 0.3, linestyle='dashed') #A7
    plt.axvline(x=866, color='grey', alpha = 0.3, linestyle='dashed') #A8
    plt.axvline(x=930, color='grey', alpha = 0.3, linestyle='dashed') #A9
    plt.axvline(x=950, color='grey', alpha = 0.3, linestyle='dashed') #A10
    plt.axvline(x=372, color='red', alpha = 0.3) #site 1
    plt.axvline(x=391, color='red', alpha = 0.3) #site 2
    
    plt.axvline(x=419, color='grey', alpha = 0.3) #X
    

    #Add data and error bars
    n = 0
    for i in range(len(complete_data)):
        n += 1
        label_dict = {1 : 'Gly', 2 : 'fhOrn', 3 : 'Ser', 4 : 'Ala', 5 : 'Arg1', 6 : 'Arg2', 7 : 'Glu', 8 : 'Phe'}
        data = complete_data[i]
        wells = range(0, len(data))
        plt.plot(wells, data, label = label_dict[n])

        axes = plt.gca()
        axes.set_xlim([0, len(data)])


    plt.legend(prop={"size":7}, loc="upper left")
    
    #Add labels
    params = {'mathtext.default': 'regular' }          
    plt.rcParams.update(params)
    axes.set_xlabel('Alignment position')
    axes.set_ylabel('Average number of clashes')
    axes.tick_params(axis='x')
    axes.tick_params(axis='y')
    
    plt.savefig('Profile_' + str(name) + "_Txo2_labelled.svg")


def wrapperIndividualGraph(parentFile, contactFile):
    '''
    Considers recombination separately for each parent/child pair. Calculates 
    the number of clashes and draws it on a graph.
    '''
    pdbName = contactFile.split('_')[0][-4:]
    parent_list = schema.readMultipleSequenceAlignmentFile(file(parentFile, 'r'))
    parents = [p for (k,p) in parent_list]

    pdb_contacts = schema.readContactFile(file(contactFile, 'r'))

    clash_data = [[] for sample in range(1, len(parents))]
    for i in range(1, len(parents)):
        print i

        newList = [parents[0], parents[i]]

        for residue in range(0, len(parents[0])):
            crossovers = [residue]  
            contacts = schema.getSCHEMAContactsWithCrossovers(pdb_contacts, newList, crossovers)
            fragments = schema.getFragments(crossovers, parents[0])
            
            clash_data[i-1].append(schema.getChimeraDisruption('21', contacts, fragments, newList))

    makeIndividualBarGraph(clash_data, pdbName)

wrapperIndividualGraph('Pa11_MSA_C1-A10.aln', 'Pa11_contacts.txt')

