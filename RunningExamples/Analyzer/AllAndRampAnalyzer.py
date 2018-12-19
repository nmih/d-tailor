from SequenceAnalyzer import SequenceAnalyzer

from Features.CAI import CAI
from Features.StructureRNAFold import StructureRNAFold, StructureRNAFoldMFE
from Features.RNADuplexRNAFold import RNADuplexRNAFold, RNADuplexRNAFoldRibosome, RNADuplexRNAFoldMFE

from Functions import validateCDS


class AllAndRampAnalyzer(SequenceAnalyzer):
    '''
    Class to analyze CAI and RNA structure
    '''

    def __init__(self, input_file, input_type, cai_table, ramp_from_to, root_dir):
        SequenceAnalyzer.__init__(self,
                                  input_file,
                                  input_type,
                                  root_dir=root_dir)
        self.cai_table = cai_table
        self.ramp_from_to = ramp_from_to

    def configureSolution(self, solution):

        solution.valid = validateCDS(solution.sequence)

        if solution.valid:

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

    def outputStart(self):
        # print "gene_name,sd_hyb_energy,mfe_structure,cai"
        pass

    def output(self, solution):
        if solution.valid:
            # print solution.solid, ',',
            # # print solution.scores['sd16sRNADuplexRNAFoldMFE'], ',',
            # print solution.scores['utrStructureRNAFoldMFE'], ',',
            # print solution.scores['cdsCAI']

            return {'mfeStructureRNAFoldMFE': solution.scores['mfeStructureRNAFoldMFE'],
                    'cai_allCAI': solution.scores['cai_allCAI'],
                    'cai_rampCAI': solution.scores['cai_rampCAI'],
                    'cai_restCAI': solution.scores['cai_restCAI']}
