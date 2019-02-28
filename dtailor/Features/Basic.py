from dtailor.Features.Feature import Feature
from dtailor import Functions, Solution
from uuid import uuid4
import Bio.SeqUtils
import logging
logger = logging.getLogger(__name__)


class LongestRepeatedSubseq(Feature):
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

    def analyze_longest_repeated_subsequence(self):
        """https://www.geeksforgeeks.org/longest-repeated-subsequence/"""

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

        return len(res)

    def set_scores(self):
        # TODO: definitely the slowest link, optimize this with a rolling hash table?
        logger.debug('Scoring longest repeated sequences for sequence')#: {}'.format(self.sequence))
        self.scores[self.label + "LongestRepeatedSubseq"] = self.analyze_longest_repeated_subsequence()

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

                    # logger.debug('x {} localgcmax {}'.format(x, max_gc))
        return max_gc

    def set_scores(self):
        logger.debug('Scoring local GC content for sequence')#: {}'.format(self.sequence))
        self.scores[self.label + "LocalGC"] = self.analyze_local_gc()

    def mutate(self):
        # TODO: smarter mutation would be in the window of non desired GC content
        # self.analyze_local_gc()
        # return Feature.randomMutation(self, mutable_region=(self.window_start_index-30, self.window_start_index+self.window_size+30))
        return Feature.randomMutation(self, mutable_region=self.mutable_region)
