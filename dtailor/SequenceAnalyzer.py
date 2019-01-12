'''
Created on Nov 1, 2012

@author: jcg
'''

import os
import sys
from dtailor.Solution import Solution
from csv import DictReader
from Bio import SeqIO
from tqdm import tqdm


class SequenceAnalyzer(object):
    
    '''
    Initializes class that analyzes sequence features 
    '''
    
    def __init__(self, input_file, input_type, root_dir, sep = ","):

        self.root_dir = root_dir
        
        if input_type == "CSV":
            self.list_of_input_sequences = self.readCSV(input_file,sep)
        elif input_type == "FASTA": 
            self.list_of_input_sequences = self.readFASTA(input_file)
        elif input_file == "GENBANK":
            self.list_of_input_sequences = self.readGENBANK(input_file)
        else:
            sys.stderr.write("The input type entered is not supported, please use one of the following: [CSV,FASTA,GENBANK]") 
    
    def readCSV(self,input_file,sep=","):
        pass
        list_seq = []
        
        reader = DictReader(open(input_file), delimiter=sep, quotechar='"')
        
        for l in reader:            
            list_seq.append(l)
            
        return list_seq

    def readFASTA(self,input_file):
        # pass
        # list_seq = []
        #
        # reader = open(input_file)
        #
        # name = ""
        # seq  = ""
        #
        # for l in reader:
        #     l=l.rstrip()
        #
        #     if l[0] == ">":
        #         name = l.split(' ')[0][1:]
        #     else:
        #         seq = l
        #         list_seq.append({ 'name' : name , 'sequence' : seq})
        #
        # return list_seq
        list_seq = []
        with open(input_file) as f:
            for s in SeqIO.parse(f, 'fasta'):
                list_seq.append({'name': s.id, 'sequence': str(s.seq)})
        return list_seq

    def readGENBANK(self):
        pass        
            
    def configureSolution(self, solution):
        pass
    
    def outputStart(self):
        pass
    
    def output(self,solution):
        pass

    def run(self):
        
        # self.outputStart()
        final_output = {}
        
        for sequence in tqdm(self.list_of_input_sequences):
            sol_id = sequence['name']
            seq = sequence['sequence']
            
            solution = Solution(sol_id = sol_id, sequence = seq, project_dir=os.path.join(self.root_dir, sol_id))
            self.configureSolution(solution)

            final_output[sol_id] = self.output(solution)

        return final_output