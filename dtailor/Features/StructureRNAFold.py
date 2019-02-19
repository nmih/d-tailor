'''
Created on Nov 16, 2011

@author: jcg
'''

from dtailor.Features.Feature import Feature
from dtailor import Functions
from uuid import uuid4
import logging
logger = logging.getLogger(__name__)


class StructureRNAFold(Feature):
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
            logger.debug('RNA range set to: {}:{}'.format(self.structure_range[0], self.structure_range[1]))
            self.sequence           = solution.sequence[self.structure_range[0]:self.structure_range[1]]
            self.mutable_region     = args['mutable_region']   if 'mutable_region' in args else solution.mutable_region
            self.cds_region         = args['cds_region']       if 'cds_region' in args else solution.cds_region
            self.keep_aa            = args['keep_aa']          if 'keep_aa' in args else solution.keep_aa
            self.set_scores()
            self.set_level()                    
        else: #copy instance
            Feature.__init__(self, featureObject=structureObject)
            self.structurefile      = structureObject.structurefile
            self.structure_range    = structureObject.structure_range         
            self.sequence           = structureObject.sequence
            self.mutable_region     = structureObject.mutable_region
            self.cds_region         = structureObject.cds_region
            self.keep_aa            = structureObject.keep_aa
            self.scores             = structureObject.scores
                            
    def set_scores(self, scoring_function=Functions.analyze_structure_rnafold):
        logger.debug('Scoring RNA structure for sequence')#: {}'.format(self.sequence))
        scoring_function(seq=self.sequence, filename=self.structurefile, project_dir=self.solution.project_dir)
                                                                     
    def mutate(self):        
        return Feature.randomMutation(self, mutable_region=self.mutable_region)
    
class StructureRNAFoldMFE(StructureRNAFold):
    """
    Manipulate the structure MFE
    """
    def __init__(self, structureObject, label = ""):
        StructureRNAFold.__init__(self, structureObject=structureObject)
        self.label = self.label + label
        self.set_scores()
        self.set_level()
        
    def set_scores(self, scoring_function=Functions.analyze_structure_mfe_rnafold):    
        self.scores.update(Functions.appendLabelToDict(scoring_function(filename=self.structurefile, project_dir=self.solution.project_dir), self.label))
