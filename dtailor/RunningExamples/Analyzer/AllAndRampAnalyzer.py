from dtailor.SequenceAnalyzer import SequenceAnalyzer

from dtailor.Features.CAI import CAI
from dtailor.Features.StructureRNAFold import StructureRNAFold, StructureRNAFoldMFE
from dtailor.Functions import validateCDS
import logging
logger = logging.getLogger(__name__)


class AllAndRampAnalyzer(SequenceAnalyzer):
    """Class to analyze the entire sequence, the ramp (defined in ramp_from_to), and the rest of the sequence. CAI is
    calculated for all three types, while MFE (minimum folding energy of the RNA secondary structure) is calculated
    for the ramp only.
    """

    def __init__(self, input_file, input_type, cai_table, root_dir, ramp=True, ramp_from_to=None,
                 check_frame=True, check_start=True, check_end_stop=True, check_within_stop=True):
        """Initialize the analyzer.

        Args:
            input_file (str): Path to genome file
            input_type (str): File type of the genome file
            cai_table (dict): Precalculated CAI table for this genome file
            ramp (bool): If a ramp should be analyzed
            ramp_from_to (tuple): Start and end of the ramp, 0-indexed - SO this means you should do (0,39) if you want
                positions 1-39 to be analyzed
            root_dir (str): Path to directory where output files will be stored (mostly temporary files for RNA calcs)
        """

        SequenceAnalyzer.__init__(self,
                                  input_file,
                                  input_type,
                                  root_dir=root_dir)
        self.cai_table = cai_table

        self.ramp = ramp
        self.ramp_from_to = ramp_from_to

        self.check_frame = check_frame
        self.check_start = check_start
        self.check_end_stop = check_end_stop
        self.check_within_stop = check_within_stop

        if not self.ramp:
            logger.info('Ramp set to False, turning off start codon checker')
            self.check_start = False
        else:
            if (ramp_from_to[1]-ramp_from_to[0]) % 3 != 0:
                raise ValueError('Invalid ramp range')

    def configureSolution(self, solution):

        solution.valid = validateCDS(solution.sequence,
                                     check_frame=self.check_frame, check_start=self.check_start,
                                     check_end_stop=self.check_end_stop, check_within_stop=self.check_within_stop)

        if solution.valid:

            if self.ramp:
                # CAI - entire sequence
                cai_all_obj = CAI(solution=solution,
                                  label="all",
                                  cai_table=self.cai_table,
                                  args={'cai_range'     : (0, len(solution.sequence)),
                                        'mutable_region': range(0, len(solution.sequence))})
                solution.add_feature(cai_all_obj)

                # CAI - ramp only
                cai_ramp_obj = CAI(solution=solution,
                                   label="ramp",
                                   cai_table=self.cai_table,
                                   args={'cai_range'     : (0, self.ramp_from_to[1]),
                                         'mutable_region': range(0, self.ramp_from_to[1])})
                solution.add_feature(cai_ramp_obj)

                # CAI - not ramp only
                cai_rest_obj = CAI(solution=solution,
                                   label="rest",
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
            else:
                # CAI - entire sequence
                cai_all_obj = CAI(solution=solution,
                                  label="all",
                                  cai_table=self.cai_table,
                                  args={'cai_range'     : (0, len(solution.sequence)),
                                        'mutable_region': range(0, len(solution.sequence))})
                solution.add_feature(cai_all_obj)

    def outputStart(self):
        # print "gene_name,sd_hyb_energy,mfe_structure,cai"
        pass

    def output(self, solution):
        if solution.valid:
            # print solution.solid, ',',
            # # print solution.scores['sd16sRNADuplexRNAFoldMFE'], ',',
            # print solution.scores['utrStructureRNAFoldMFE'], ',',
            # print solution.scores['cdsCAI']

            return {k: v for k, v in solution.scores.items()}
            # {'mfeStructureRNAFoldMFE'    : solution.scores['mfeStructureRNAFoldMFE'],
            #         'allCAI' : solution.scores['allCAI'],
            #         'rampCAI': solution.scores['rampCAI'],
            #         'restCAI': solution.scores['restCAI']}
