'''
Created on Dec 22, 2012

@author: jcg
'''

import sys
from dtailor.SequenceDesigner import SequenceDesigner
from dtailor.Features.CAI import CAI
from dtailor.Features.StructureRNAFold import StructureRNAFold, StructureRNAFoldMFE
from dtailor.Features.RNADuplexRNAFold import RNADuplexRNAFold, RNADuplexRNAFoldRibosome, RNADuplexRNAFoldMFE
from dtailor.Functions import validateCDS


class AllAndRampDesigner(SequenceDesigner):
    
    def __init__(self, name, seed, design, dbfile, cai_table, ramp_from_to, root_dir, createDB=True):
        SequenceDesigner.__init__(self, name, seed, design, dbfile, root_dir, createDB)
        self.cai_table = cai_table
        self.ramp_from_to = ramp_from_to
        
    def configureSolution(self, solution):
        '''
        Solution configuration
        '''
                
        if not solution.sequence:
            return 0

        # Populate solution with desired features
        # Everything can be mutated, everything is coding
        solution.mutable_region = range(0, len(solution.sequence))
        solution.cds_region = (0, len(solution.sequence))

        # This should always be true unless you want to make mutant proteins
        solution.keep_aa = True

        # CAI - entire sequence
        cai_all_obj = CAI(solution=solution,
                          label="cai_all",
                          cai_table=self.cai_table,
                          args={'cai_range'     : (0, len(solution.sequence)),
                                'mutable_region': range(0, len(solution.sequence))})
        solution.add_feature(cai_all_obj)

        # CAI - ramp only
        cai_ramp_obj = CAI(solution=solution,
                           label="cai_ramp",
                           cai_table=self.cai_table,
                           args={'cai_range'     : (0, self.ramp_from_to[1]),
                                 'mutable_region': range(0, self.ramp_from_to[1])})
        solution.add_feature(cai_ramp_obj)

        # CAI - not ramp only
        cai_rest_obj = CAI(solution=solution,
                           label="cai_rest",
                           cai_table=self.cai_table,
                           args={'cai_range'     : (self.ramp_from_to[1], len(solution.sequence)),
                                 'mutable_region': range(self.ramp_from_to[1], len(solution.sequence))})
        solution.add_feature(cai_rest_obj)

        # MFE - ramp only
        st1_obj = StructureRNAFold(solution=solution,
                                   label="mfe",
                                   args={'structure_range': (self.ramp_from_to[0], self.ramp_from_to[1]),
                                         'mutable_region' : range(self.ramp_from_to[0], self.ramp_from_to[1])})
        st_mfe = StructureRNAFoldMFE(structureObject=st1_obj)
        st1_obj.add_subfeature(st_mfe)
        solution.add_feature(st1_obj)
    
    def validateSolution(self, solution):
        '''
        Solution validation tests
        '''
        if solution.sequence == None or ('?' in solution.levels.values()):
            sys.stderr.write("SolutionValidator: Level unknown - "+str(solution.levels)+"\n")                        
            solution.valid = False
            return 0

        # Basic validation
        designed_region = solution.sequence
        valid = validateCDS(designed_region)

        # # No internal Promoters - e coli specific
        # (score, _, _) = Functions.look_for_promoters(designed_region)
        # if score >= 15.3990166: #0.95 percentile for Promoter PWM scores
        #     valid = False
        #     sys.stderr.write("SolutionValidator: High Promoter score: "+str(score)+"\n")
        
        # # No internal Terminator - only for bacterial genomes?
        # score = Functions.look_for_terminators(seq=designed_region, outdir=self.root_dir)
        # if score >= 90: #90% confidence from transtermHP
        #     valid = False
        #     sys.stderr.write("SolutionValidator: High Terminator score\n")
            
        # # No restriction enzymes
        # if 'ggtctc' in designed_region or 'gagacc' in designed_region:
        #    sys.stderr.write("SolutionValidator: Restriction enzyme found\n")
        #    valid = False        
        
        solution.valid = valid
        
        return valid
