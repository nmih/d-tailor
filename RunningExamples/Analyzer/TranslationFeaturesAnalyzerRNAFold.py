from SequenceAnalyzer import SequenceAnalyzer

from Features.CAI import CAI
from Features.StructureRNAFold import StructureRNAFold, StructureRNAFoldMFE
from Features.RNADuplexRNAFold import RNADuplexRNAFold, RNADuplexRNAFoldRibosome, RNADuplexRNAFoldMFE

from Functions import validateCDS


class TranslationFeaturesAnalyzer(SequenceAnalyzer):
    '''
    Class to analyze CAI and RNA structure
    '''

    def __init__(self, input_file, input_type, cai_table, mfe_from_to, root_dir):
        SequenceAnalyzer.__init__(self,
                                  input_file,
                                  input_type,
                                  root_dir=root_dir)
        self.cai_table = cai_table
        self.mfe_from_to = mfe_from_to

    def configureSolution(self, solution):

        # TODO: what's validateCDS do?
        solution.valid = validateCDS(solution.sequence)

        if solution.valid:

            # TODO: read about ranges selected...

            # CAI
            cai_obj = CAI(solution=solution,
                          label="cds",
                          cai_table=self.cai_table,
                          args={'cai_range': (0, len(solution.sequence))})
            solution.add_feature(cai_obj)

            # MFE entire sequence
            st1_obj = StructureRNAFold(solution=solution,
                                       label="utr",
                                       args={'structure_range': (self.mfe_from_to[0], self.mfe_from_to[1])})
            st_mfe = StructureRNAFoldMFE(structureObject=st1_obj)
            st1_obj.add_subfeature(st_mfe)
            solution.add_feature(st1_obj)

            # # RBS
            # dup_obj1 = RNADuplexRNAFoldRibosome(solution1=solution,
            #                                     label="sd16s",
            #                                     args={'rnaMolecule1region': (25, 48)})
            # dup_mfe = RNADuplexRNAFoldMFE(dup_obj1)
            # dup_obj1.add_subfeature(dup_mfe)
            # solution.add_feature(dup_obj1)

    def outputStart(self):
        # print "gene_name,sd_hyb_energy,mfe_structure,cai"
        print "gene_name,mfe_structure,cai"

    def output(self, solution):
        if solution.valid:
            # print solution.solid, ',',
            # # print solution.scores['sd16sRNADuplexRNAFoldMFE'], ',',
            # print solution.scores['utrStructureRNAFoldMFE'], ',',
            # print solution.scores['cdsCAI']

            return {'utrStructureRNAFoldMFE': solution.scores['utrStructureRNAFoldMFE'],
                    'cdsCAI': solution.scores['cdsCAI']}