'''
Created on Dec 22, 2012

@author: jcg
'''

import sys
import os.path as op
from dtailor.SequenceDesigner import SequenceDesigner
from dtailor.Features.CAI import CAI
from dtailor.Features.StructureRNAFold import StructureRNAFold, StructureRNAFoldMFE
from dtailor.Features.RNADuplexRNAFold import RNADuplexRNAFold, RNADuplexRNAFoldRibosome, RNADuplexRNAFoldMFE
from dtailor.Functions import validateCDS


class GenericDesigner(SequenceDesigner):
    
    def __init__(self, name, seed, design, cai_table, root_dir,
                 mutable_region=None, cds_region=None, keep_aa=True,
                 createDB=True):
        """Initialize a project to design a sequence according to design parameters.

        Args:
            name (str): Identifier of the sequence
            seed (str): DNA sequence
            design (dict): Parameters for design
            cai_table (dict): Precalculated CAI table for your organism
            root_dir (str): Output directory for SQL and temp files
            mutable_region (tuple): What can be mutated (0-index)
            cds_region (tuple): What is coding TODO: does this really impact anything?
            createDB (bool): TODO: Usually True? What happens if it's not?
        """
        if not op.exists(root_dir):
            raise IOError('{} does not exist, please make it first'.format(root_dir))

        dbfile = op.join(root_dir, name)
        SequenceDesigner.__init__(self, name, seed, design, dbfile, root_dir, createDB)

        self.cai_table = cai_table
        self.keep_aa = keep_aa

        if not mutable_region:
            self.mutable_region = range(0, len(seed))
        else:
            self.mutable_region = range(mutable_region[0], mutable_region[1])

        if not cds_region:
            self.cds_region = (0, len(seed))
        else:
            self.cds_region = cds_region

    def configureSolution(self, solution):
        '''
        Solution configuration
        '''
                
        if not solution.sequence:
            return 0

        # Populate solution with desired features
        solution.mutable_region = self.mutable_region
        solution.cds_region = self.cds_region

        # This should always be true unless you want to make mutant proteins
        solution.keep_aa = self.keep_aa

        # TODO: loop over self.features (input as "design" and saved when SequenceDesigner is init and add features using the info in that
        # Design param will have more stuff in it like so
        # design_param = {
        #     "mfeStructureRNAFoldMFE": {
        #         'feattype'      : 'MFE',
        #         'type'          : 'REAL',
        #         'mutable_region': (0, 30),
        #         'thresholds'    : mfe_levels},
        #     "cai_rampCAI"           : {
        #         'feattype'      : 'CAI',
        #         'type'          : 'REAL',
        #         'mutable_region': (0, 30),
        #         'thresholds'    : cai_levels},
        #     "cai_restCAI"           : {
        #         'feattype'      : 'CAI',
        #         'type'          : 'REAL',
        #         'mutable_region': (30, len(sequence)),
        #         'thresholds'    : cai_levels}
        # }

        for feat, params in self.features:
            if params['feattype'] == 'CAI':
                a_feature = CAI(solution=solution,
                                label=feat,
                                cai_table=self.cai_table,
                                args={'cai_range'     : (params['mutable_region'][0],
                                                         params['mutable_region'][1]),
                                      'mutable_region': range(params['mutable_region'][0],
                                                              params['mutable_region'][1])})
            elif params['feattype'] = 'MFE':
                a_feature = StructureRNAFold(solution=solution,
                                             label=feat,
                                             args={'structure_range': (params['mutable_region'][0],
                                                                       params['mutable_region'][1]),
                                                   'mutable_region' : range(params['mutable_region'][0],
                                                                            params['mutable_region'][1])})
                st_mfe = StructureRNAFoldMFE(structureObject=a_feature)
                a_feature.add_subfeature(st_mfe)
            else:
                print('Feature {} not supported yet'.format(feat))

            solution.add_feature(a_feature)


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
