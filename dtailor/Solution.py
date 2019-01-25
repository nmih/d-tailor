'''
Created on Nov 1, 2011

@author: jcg
'''

import sys
from random import choice
from dtailor.Functions import randomMutationOperator
from uuid import uuid4

class Solution:
    '''
    A Solution encapsulates a sequence and their inherent attributes:
        sol_id - ID for Solution
        seqence - sequence for Solution
        cds_region - a tuple indicating the location of (Begin,End) of CDS sequence (this will be necessary in the design mode if one want to contrain mutations).
        mutable_region - a list with all positions that can be mutated
        parent - Solution from which the current Solution was derived
        
    '''
    # def __init__(self, sol_id = 0, sequence="", cds_region = None, keep_aa = False, mutable_region = None, parent=None, design = None):
    #
    #     if sequence == None:
    #         sys.stderr.write("Tried to create a solution with sequence NULL\n")
    #         self.sequence = None
    #         return None

    def __init__(self, sol_id, sequence, project_dir, cds_region=None, keep_aa=False, mutable_region=None, parent=None,
                 design=None):

        #check if solution is in DB        
        self.mutable_region     = mutable_region                        
        self.cds_region         = cds_region
        self.keep_aa            = keep_aa
        self.solid              = sol_id
        self.parent             = parent
        self.sequence           = sequence.lower()
        self.scores             = {}
        self.levels             = {}
        self.features           = {}
        self.designMethod       = design
        self.valid = True
        self.project_dir = project_dir

        self.new_features_list = []
        if design:
            self.new_features_list = [feature + param['feattype'] for feature, param in self.designMethod.features.items()]
        self.cai_table = None

    # def add_feature(self, feature):
    #     # featureLabel = feature.label + feature.__class__.__name__
    #     if not self.features.has_key(featureLabel):
    #         self.features[featureLabel] = feature
    #         # update scores
    #         self.scores.update(feature.scores)
    #         # update levels
    #         if feature.level != None:
    #             self.levels[featureLabel + "Level"] = feature.level
    #         for subfeature in feature.subfeatures.values():
    #             self.add_feature(subfeature)
    #     else:
    #         sys.stderr.write("Feature label already exists!")
    #
    #     return

    def add_feature(self, feature):
        featureLabel = feature.label + feature.__class__.__name__
        if not featureLabel in self.features:
            self.features[featureLabel] = feature
            # update scores
            self.scores.update(feature.scores)
            # update levels
            if feature.level != None:
                self.levels[featureLabel + "Level"] = feature.level
            for subfeature in feature.subfeatures.values():
                self.add_feature(subfeature)
        else:
            sys.stderr.write("Feature label already exists!")

        if feature.__class__.__name__ == 'CAI':
            self.cai_table = feature.cai_table

        return

    def checkSolution(self, desiredSolution):
        
        if desiredSolution == None:
            return False
        
        same = True
        for feature in self.new_features_list:
            key = feature+"Level"            
            same = same & (desiredSolution[key]==0 or desiredSolution[key]==self.levels[key]) 
            
        return same            
        
    def mutate(self,desiredSolution=None,random=False):
                
        if desiredSolution==None or random or self.designMethod.listDesigns == [] or self.features == {}:
            return self.randomMutation()
        else:        
            # get features with targets
            mutable = []
            for feature in self.features.values():
                if feature.defineTarget(desiredSolution):
                    mutable.append(feature)
        
            if mutable == []:
                return None                    
            
            rm =  choice(mutable)
            #tomutatefeatures = [k.label+k.__class__.__name__ for k in mutable]                        
            #print(tomutatefeatures)
            #print([self.scores[k] for k in tomutatefeatures])
            #print([self.levels[k+"Level"] for k in tomutatefeatures])
            #print([desiredSolution[k+"Level"] for k in tomutatefeatures])
            #print("mutating... " + rm.label+rm.__class__.__name__)
            
            return rm.mutate()
            #return choice(mutable).randomMutation()
            #return self.randomMutation()
        
    def randomMutation(self,pos=None,n_mut=[1,2]):
        new_seq = randomMutationOperator(self.sequence, self.keep_aa, self.mutable_region, self.cds_region, pos,n_mut)     
        return Solution(sol_id=str(uuid4().int), sequence=new_seq, cds_region = self.cds_region, keep_aa = self.keep_aa, mutable_region = self.mutable_region, parent=self, design=self.designMethod)