'''
Created on Nov 16, 2011

@author: jcg
'''

from dtailor.Features.Feature import Feature
from dtailor import Functions, Solution
from uuid import uuid4

class Structure(Feature):
    """
    Structure Feature
        solution - solution where structure should be computed
        label - some label to append to the name of structure file
        structure_range - start and end position to calculate structure - a tuple in the form (start, end)  
        mutable_region - a list with all bases that can be mutated
        mutable_region - a list with all bases that can be mutated
        cds_region - a pair with begin and end of CDSs - example: (0,100)
        keep_aa - boolean option indicating if in the design mode amino acids should be kept
    """    
    def __init__(self, structureObject = None, solution = None, label="", args = { 'structure_range' : (0,59),
                                                                                   'mutable_region' : None, 
                                                                                   'cds_region' : None, 
                                                                                   'keep_aa' : True }):
        
        if structureObject == None: #create new instance
            #General properties of feature
            Feature.__init__(self, solution=solution, label=label)
            #Specifics of this Feature
            self.structurefile      = solution.solid + label
            self.structure_range    = args['structure_range']        
            self.sequence           = solution.sequence[self.structure_range[0]:self.structure_range[1]]
            self.mutable_region     = args['mutable_region']   if 'mutable_region' in args else solution.mutable_region
            self.cds_region         = args['cds_region']       if 'cds_region' in args else solution.cds_region
            self.keep_aa            = args['keep_aa']          if 'keep_aa' in args else solution.keep_aa
            self.set_scores()
            self.set_level()                    
        else: #copy instance
            Feature.__init__(self, structureObject)
            self.structurefile      = structureObject.structurefile
            self.structure_range    = structureObject.structure_range         
            self.sequence           = structureObject.sequence
            self.mutable_region     = structureObject.mutable_region
            self.cds_region         = structureObject.cds_region
            self.keep_aa            = structureObject.keep_aa
            self.scores             = structureObject.scores

                            
    def set_scores(self, scoring_function=Functions.analyze_structure): 
        scoring_function(seq=self.sequence, filename=self.structurefile, project_dir=self.solution.project_dir)
                                                                     
    def mutate(self, operator=Functions.SimpleStructureOperator):        
        if not self.targetInstructions:
            return None        
        ss_bases = None if not self.label+'StructureSingleStrandedBasesList' in self.scores else self.scores[self.label+'StructureSingleStrandedBasesList']
        ds_bases = None if not self.label+'StructureDoubleStrandedBasesList' in self.scores else self.scores[self.label+'StructureDoubleStrandedBasesList']
                    
        new_seq = operator(self.solution.sequence, self.structurefile,  self.structure_range, self.mutable_region, self.cds_region, self.targetInstructions['direction'], ss_bases=ss_bases, ds_bases=ds_bases)
        if not new_seq:
            return None                 
        return Solution.Solution(sol_id=str(uuid4().int), sequence=new_seq, cds_region = self.cds_region, mutable_region = self.mutable_region, parent=self.solution, design=self.solution.designMethod)
    
class StructureMFE(Structure):
    """
    Manipulate the structure MFE
    """
    def __init__(self, structureObject, label = ""):
        Structure.__init__(self, structureObject=structureObject)
        self.label = self.label + label
        self.set_scores()
        self.set_level()      
        
    def set_scores(self, scoring_function=Functions.analyze_structure_mfe):    
        self.scores.update(Functions.appendLabelToDict(scoring_function(filename=self.structurefile, project_dir=self.solution.project_dir), self.label))


class StructureSingleStranded(Structure):
    """
    Manipulate the structure single stranded bases
    """
    def __init__(self, structureObject, label = ""):
        Structure.__init__(self,structureObject)
        self.label = self.label + label
        self.set_scores()
        self.set_level()

    def set_scores(self, scoring_function=Functions.analyze_structure_ss):    
        self.scores.update(Functions.appendLabelToDict(scoring_function(self.structurefile), self.label))
        
class StructureDoubleStranded(Structure):
    """
    Manipulate the structure double stranded bases
    """
    def __init__(self, structureObject, label = ""):
        Structure.__init__(self,structureObject)
        self.label = self.label + label
        self.set_scores()
        self.set_level()

    def set_scores(self, scoring_function=Functions.analyze_structure_ds):    
        self.scores.update(Functions.appendLabelToDict(scoring_function(self.structurefile), self.label))
        
    def defineTarget(self,desiredSolution):
        '''
        Function that determines if a target wasn't hit, and if not updates targetDirections 
        '''
        if desiredSolution == None:
            return True
        
        if Structure.defineTarget(self, desiredSolution):
            if self.targetInstructions['direction'] == '+':
                self.targetInstructions['direction'] = '-' #decrease number of single bases (increase double stranded bases)
            else:
                self.targetInstructions['direction'] = '+' #increase number of single bases (decrease double stranded bases)                
