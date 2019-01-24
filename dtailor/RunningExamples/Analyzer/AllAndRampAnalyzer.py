from dtailor.SequenceAnalyzer import SequenceAnalyzer

from dtailor.Features.CAI import CAI
from dtailor.Features.StructureRNAFold import StructureRNAFold, StructureRNAFoldMFE
from dtailor.Functions import validateCDS


class AllAndRampAnalyzer(SequenceAnalyzer):
    """Class to analyze the entire sequence, the ramp (defined in ramp_from_to), and the rest of the sequence. CAI is
    calculated for all three types, while MFE (minimum folding energy of the RNA secondary structure) is calculated
    for the ramp only.
    """

    def __init__(self, input_file, input_type, cai_table, ramp_from_to, root_dir):
        """Initialize the analyzer.

        Args:
            input_file (str): Path to genome file
            input_type (str): File type of the genome file
            cai_table (dict): Precalculated CAI table for this genome file
            ramp_from_to (tuple): Start and end of the ramp, 0-indexed
            root_dir (str): Path to directory where output files will be stored (mostly temporary files for RNA calcs)
        """

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

    def outputStart(self):
        # print "gene_name,sd_hyb_energy,mfe_structure,cai"
        pass

    def output(self, solution):
        if solution.valid:
            # print solution.solid, ',',
            # # print solution.scores['sd16sRNADuplexRNAFoldMFE'], ',',
            # print solution.scores['utrStructureRNAFoldMFE'], ',',
            # print solution.scores['cdsCAI']

            return {'mfe'    : solution.scores['mfeStructureRNAFoldMFE'],
                    'allCAI' : solution.scores['allCAI'],
                    'rampCAI': solution.scores['rampCAI'],
                    'restCAI': solution.scores['restCAI']}
