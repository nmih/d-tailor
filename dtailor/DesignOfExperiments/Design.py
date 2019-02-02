"""
Created on Nov 11, 2012

@author: jcg
"""

from itertools import product
from collections import OrderedDict
import logging
logger = logging.getLogger(__name__)


class Design(object):
    """
    Abstract Class defining experimental design to be used
    featuresObj is a dictionary containing relevant info about design such as variables to use and thresholds for feature values
    """

    def __init__(self, featuresObj):
        # validate if ids only have one character
        for feat in featuresObj.keys():
            for l_id in featuresObj[feat]['thresholds']:
                if len(str(l_id)) != 1:
                    raise Exception("Please use only one digit/character as the level identifier!\n")

        self.featuresList = list(featuresObj.keys())
        self.n_features = len(self.featuresList)
        self.features = featuresObj
        self.thresholds = {feature: featuresObj[feature]['thresholds'] for feature in featuresObj.keys()}


class Optimization(Design):
    """
    Class encoding a single-target design (optimization)
    """

    def __init__(self, featuresObj, target):
        Design.__init__(self, featuresObj)
        self.listDesigns = [target]
        self.nDesigns = self.listDesigns.__len__()


class RandomSampling(Design):
    """
    Class encoding a random sampling design
    """

    def __init__(self, featuresObj, sample_size=1000):
        Design.__init__(self, featuresObj)
        self.listDesigns = []
        self.nDesigns = sample_size


class FullFactorial(Design):
    """
    Class encoding a multi factorial design
    """

    def __init__(self, featuresObj):
        Design.__init__(self, featuresObj)
        self.listDesigns = self.computeCombinations()
        self.nDesigns = self.listDesigns.__len__()

    def computeCombinations(self):
        pass
        return list(map(".".join, eval("product('" + "','".join(
                [''.join([str(x) for x in self.thresholds[feat]]) for feat in self.featuresList]) + "')")))


class CustomDesign(Design):
    """
    Class encoding a custom design (as many targets as you want)
    """

    def __init__(self, featuresObj, targets=[]):
        Design.__init__(self, featuresObj)

        self.thresholds = OrderedDict()
        for feature in featuresObj.keys():
            self.thresholds[feature + featuresObj[feature]['feattype']] = featuresObj[feature]['thresholds']

        self.listDesigns = targets
        self.nDesigns = len(targets)

        logger.debug('Designs/targets: {}'.format(self.listDesigns))