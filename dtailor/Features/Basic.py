from dtailor.Features.Feature import Feature
from dtailor import Functions, Solution
from uuid import uuid4
import Bio.SeqUtils
import random
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import os
import logging
logger = logging.getLogger(__name__)


class LongestRepeatedSubstr(Feature):
    """Basic boolean check for if there are subsequences that """

    def __init__(self, featureObject=None, solution=None, label="", args={'feature_range' : (0, 59),
                                                                          'mutable_region': None,
                                                                          'cds_region': None,
                                                                          'keep_aa'   : True}):
        if featureObject == None:  # create new instance
            # General properties of feature
            Feature.__init__(self, solution=solution, label=label)
            # Specifics of this Feature
            self.feature_range = args['feature_range']
            self.sequence = solution.sequence[self.feature_range[0]:self.feature_range[1]]
            self.mutable_region = args['mutable_region'] if 'mutable_region' in args else solution.mutable_region
            self.cds_region = args['cds_region'] if 'cds_region' in args else solution.cds_region
            self.keep_aa = args['keep_aa'] if 'keep_aa' in args else solution.keep_aa
            self.set_scores()
            self.set_level()
        else:  # copy instance
            Feature.__init__(self, featureObject)
            self.feature_range = featureObject.feature_range
            self.sequence = featureObject.sequence
            self.mutable_region = featureObject.mutable_region
            self.cds_region = featureObject.cds_region
            self.keep_aa = featureObject.keep_aa
            self.scores = featureObject.scores

    def analyze_longest_repeated_substring(self):
        """https://www.geeksforgeeks.org/longest-repeating-and-non-overlapping-substring/"""

        n = len(self.sequence)
        LCSRe = [[0 for x in range(n + 1)]
                    for y in range(n + 1)]

        res = "" # To store result
        res_length = 0 # To store length of result

        # building table in bottom-up manner
        index = 0
        for i in range(1, n + 1):
            for j in range(i + 1, n + 1):

                # (j-i) > LCSRe[i-1][j-1] to remove
                # overlapping
                if self.sequence[i - 1] == self.sequence[j - 1] and LCSRe[i - 1][j - 1] < (j - i):
                    LCSRe[i][j] = LCSRe[i - 1][j - 1] + 1

                    # updating maximum length of the
                    # substring and updating the finishing
                    # index of the suffix
                    if (LCSRe[i][j] > res_length):
                        res_length = LCSRe[i][j]
                        index = max(i, index)

                else:
                    LCSRe[i][j] = 0

        # If we have non-empty result, then insert
        # all characters from first character to
        # last character of string
        if (res_length > 0):
            for i in range(index - res_length + 1,
                                        index + 1):
                res = res + self.sequence[i - 1]

        # Set mutable region to this repeated substring
        mutstart = self.sequence.index(res)
        self.mutable_region = range(mutstart, mutstart + len(res))
        logger.debug('Repeated substring at index {} of len {}'.format(mutstart, len(res)))

        return len(res)

    def set_scores(self):
        # TODO: #later: This part takes the longest to score (not that long really) and can be sped up with a suffix tree or suffix array
        logger.debug('Scoring longest repeated sequences for sequence: {}'.format(self.sequence))
        self.scores[self.label + "LongestRepeatedSubstr"] = self.analyze_longest_repeated_substring()

    def mutate(self):
        return Feature.randomMutation(self, mutable_region=self.mutable_region)


class LongestHomopolymer(Feature):
    """Find the longest length of a homopolymer repeats (ie. ATAAAAAC -> 5 (AAAAA))"""

    def __init__(self, featureObject=None, solution=None, label="", args={'feature_range' : (0, 59),
                                                                          'mutable_region': None,
                                                                          'cds_region': None,
                                                                          'keep_aa'   : True}):
        if featureObject == None:  # create new instance
            # General properties of feature
            Feature.__init__(self, solution=solution, label=label)
            # Specifics of this Feature
            self.feature_range = args['feature_range']
            self.sequence = solution.sequence[self.feature_range[0]:self.feature_range[1]]
            self.mutable_region = args['mutable_region'] if 'mutable_region' in args else solution.mutable_region
            self.cds_region = args['cds_region'] if 'cds_region' in args else solution.cds_region
            self.keep_aa = args['keep_aa'] if 'keep_aa' in args else solution.keep_aa
            self.set_scores()
            self.set_level()
        else:  # copy instance
            Feature.__init__(self, featureObject)
            self.feature_range = featureObject.feature_range
            self.sequence = featureObject.sequence
            self.mutable_region = featureObject.mutable_region
            self.cds_region = featureObject.cds_region
            self.keep_aa = featureObject.keep_aa
            self.scores = featureObject.scores

    def analyze_longest_homopolymer(self):
        """Get the length of the longest homopolymer within a sequence.

        https://stackoverflow.com/questions/31084639/finding-the-length-of-longest-repeating

        Returns:
            int: Length of longest homopolymer

        """
        len_substring = 0
        longest = 0
        for i in range(len(self.sequence)):
            if i > 1:
                if self.sequence[i] != self.sequence[i - 1]:
                    len_substring = 0
            len_substring += 1
            if len_substring > longest:
                longest = len_substring
        return longest

    def set_scores(self):
        logger.debug('Scoring longest homopolymer for sequence')#: {}'.format(self.sequence))
        self.scores[self.label + "LongestHomopolymer"] = self.analyze_longest_homopolymer()

    def mutate(self):
        return Feature.randomMutation(self, mutable_region=self.mutable_region)


class GlobalGC(Feature):
    """Basic boolean check for if there are subsequences that """

    def __init__(self, featureObject=None, solution=None, label="", args={'feature_range' : (0, 59),
                                                                          'mutable_region': None,
                                                                          'cds_region': None,
                                                                          'keep_aa'   : True}):
        if featureObject == None:  # create new instance
            # General properties of feature
            Feature.__init__(self, solution=solution, label=label)
            # Specifics of this Feature
            self.feature_range = args['feature_range']
            self.sequence = solution.sequence[self.feature_range[0]:self.feature_range[1]]
            self.mutable_region = args['mutable_region'] if 'mutable_region' in args else solution.mutable_region
            self.cds_region = args['cds_region'] if 'cds_region' in args else solution.cds_region
            self.keep_aa = args['keep_aa'] if 'keep_aa' in args else solution.keep_aa
            self.set_scores()
            self.set_level()
        else:  # copy instance
            Feature.__init__(self, featureObject)
            self.feature_range = featureObject.feature_range
            self.sequence = featureObject.sequence
            self.mutable_region = featureObject.mutable_region
            self.cds_region = featureObject.cds_region
            self.keep_aa = featureObject.keep_aa
            self.scores = featureObject.scores

    def analyze_global_gc(self):
        """Returns a float from 0 to 100 representing percent GC content"""
        return Bio.SeqUtils.GC(self.sequence)

    def set_scores(self):
        logger.debug('Scoring global GC content for sequence')#: {}'.format(self.sequence))
        self.scores[self.label + "GlobalGC"] = self.analyze_global_gc()

    def mutate(self):
        return Feature.randomMutation(self, mutable_region=self.mutable_region)


class LocalGC(Feature):
    """Basic boolean check for if there are subsequences that """

    def __init__(self, window_size=50,
                 featureObject=None,
                 solution=None,
                 label="",
                 args={'feature_range' : (0, 59),
                       'mutable_region': None,
                       'cds_region'    : None,
                       'keep_aa'       : True}):
        if featureObject == None:  # create new instance
            # General properties of feature
            Feature.__init__(self, solution=solution, label=label)
            # Specifics of this Feature
            self.window_size = window_size
            self.feature_range = args['feature_range']
            self.sequence = solution.sequence[self.feature_range[0]:self.feature_range[1]]
            self.mutable_region = args['mutable_region'] if 'mutable_region' in args else solution.mutable_region
            self.cds_region = args['cds_region'] if 'cds_region' in args else solution.cds_region
            self.keep_aa = args['keep_aa'] if 'keep_aa' in args else solution.keep_aa
            self.set_scores()
            self.set_level()
        else:  # copy instance
            Feature.__init__(self, featureObject)
            self.feature_range = featureObject.feature_range
            self.sequence = featureObject.sequence
            self.mutable_region = featureObject.mutable_region
            self.cds_region = featureObject.cds_region
            self.keep_aa = featureObject.keep_aa
            self.scores = featureObject.scores

        self.window_start_index = 0

    def analyze_local_gc(self):
        """Return a percentage (0 to 100) of the window with the highest local GC%"""
        max_gc = 0
        if len(self.sequence) <= self.window_size:
            return Bio.SeqUtils.GC(self.sequence)

        for x in range(len(self.sequence)):
            if x <= len(self.sequence) - self.window_size:
                gc = Bio.SeqUtils.GC(self.sequence[x:x + self.window_size])

                if gc > max_gc:
                    max_gc = gc
                    self.window_start_index = x

        # Change the mutable region to be in this window of high GC%
        self.mutable_region = range(self.window_start_index, self.window_start_index+self.window_size)
        logger.debug('Local GC content of percentage {} at index {} to {}'.format(max_gc, self.window_start_index, self.window_start_index+self.window_size))
        return max_gc

    def set_scores(self):
        logger.debug('Scoring local GC content for sequence')#: {}'.format(self.sequence))
        self.scores[self.label + "LocalGC"] = self.analyze_local_gc()

    def mutate(self):
        return Feature.randomMutation(self, mutable_region=self.mutable_region)


class SmallRepeatPercentage(Feature):
    """Basic boolean check for if there are subsequences that """

    def __init__(self, featureObject=None,
                 solution=None,
                 label="",
                 args={'feature_range' : (0, 59),
                       'mutable_region': None,
                       'cds_region'    : None,
                       'keep_aa'       : True}):
        if featureObject == None:  # create new instance
            # General properties of feature
            Feature.__init__(self, solution=solution, label=label)
            # Specifics of this Feature
            self.feature_range = args['feature_range']
            self.sequence = solution.sequence[self.feature_range[0]:self.feature_range[1]]
            self.mutable_region = args['mutable_region'] if 'mutable_region' in args else solution.mutable_region
            self.cds_region = args['cds_region'] if 'cds_region' in args else solution.cds_region
            self.keep_aa = args['keep_aa'] if 'keep_aa' in args else solution.keep_aa
            self.set_scores()
            self.set_level()
        else:  # copy instance
            Feature.__init__(self, featureObject)
            self.feature_range = featureObject.feature_range
            self.sequence = featureObject.sequence
            self.mutable_region = featureObject.mutable_region
            self.cds_region = featureObject.cds_region
            self.keep_aa = featureObject.keep_aa
            self.scores = featureObject.scores

    def write_tmp_fasta_file(self, alphabet=IUPAC.IUPACUnambiguousDNA):
        """Write a temporary FASTA file"""
        sr = SeqRecord(Seq(self.sequence, alphabet), id="tmp_id", name="tmp_name", description="tmp_desc", dbxrefs=None,
                       features=None, annotations=None, letter_annotations=None)
        outfile = 'tmp_repfind.fasta'
        SeqIO.write(sr, outfile, "fasta")
        return outfile

    def run_repfind(self, in_fasta, min_repeat_len=8, max_repeat_len=20):
        """Run the REPFIND program and load results into memory"""
        executed = subprocess.run(
                ['repfind', '-s {}'.format(min_repeat_len), '-l {}'.format(max_repeat_len), '-f', in_fasta],
                stdout=subprocess.PIPE)
        if executed.returncode != 0:
            raise Exception('Error code {}'.format(executed.returncode))
        result = executed.stdout.decode('utf-8')
        result = result.splitlines()[14:]

        pattern_to_locations = {}

        for x in range(0, len(result), 4):
            pattern = result[x].replace('Word: ', '')
            locations = [int(y) for y in result[x + 1].replace('Locations: ', '').replace('|', '').split()]
            pattern_to_locations[pattern] = locations

        os.remove(in_fasta)

        return pattern_to_locations

    def get_repeat_percent_composition(self, repeat_dict):
        """Get the percentage of the sequence composed of small repeats.
        NOTE: this will actually overestimate the repeats as it doesn't care about overlaps - however that's ok for
        current purposes"""
        total_num_repeated_bases = 0
        for pat, locs in repeat_dict.items():
            pat_total_len = len(pat) * len(locs)
            total_num_repeated_bases += pat_total_len

        # Randomly select one of these repeats to target for mutagenesis
        random_key = random.choice(list(repeat_dict.keys()))
        if len(repeat_dict[random_key]) != 0:
            random_loc = random.choice(repeat_dict[random_key])
            self.mutable_region = range(random_loc, random_loc + len(random_key))

        return (total_num_repeated_bases / len(self.sequence))*100

    def set_scores(self):
        """Write the FASTA, run REPFIND, and get the percentage"""
        logger.debug('Scoring percentage of small repeats for sequence')  #: {}'.format(self.sequence))
        tmp_fasta = self.write_tmp_fasta_file()
        result = self.run_repfind(tmp_fasta)
        percent = self.get_repeat_percent_composition(result)
        self.scores[self.label + "SmallRepeatPercentage"] = percent

    def mutate(self):
        return Feature.randomMutation(self, mutable_region=self.mutable_region)


