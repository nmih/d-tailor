import pytest
from dtailor.DesignOfExperiments.Design import Design, CustomDesign


def test_design_init(design_params, pk_mfe_levels, pk_cai_levels, levels_longestrepeat, levels_longesthomopolymer, levels_globalgc, levels_localgc, levels_smallrepeats):
    des = Design(featuresObj=design_params)
    assert des.featuresList == ['mfe', 'ramp', 'rest', 'twist_lrs', 'twist_lh', 'twist_ggc', 'twist_lgc', 'twist_smr']
    assert des.n_features == 8
    assert des.features == design_params
    assert des.thresholds == {
        'mfe' : pk_mfe_levels,
        'ramp': pk_cai_levels,
        'rest': pk_cai_levels,
        'twist_lrs': levels_longestrepeat,
        'twist_lh': levels_longesthomopolymer,
        'twist_ggc': levels_globalgc,
        'twist_lgc': levels_localgc,
        'twist_smr': levels_smallrepeats,
    }


def test_custom_design_init(design_params, targets):
    des = CustomDesign(featuresObj=design_params,
                       targets=targets)
    assert des.listDesigns == targets
    assert des.nDesigns == len(targets)
    assert [k for k in des.thresholds] == ['mfeStructureRNAFoldMFE', 'rampCAI', 'restCAI',
                                           'twist_lrsLongestRepeatedSubstr', 'twist_lhLongestHomopolymer',
                                           'twist_ggcGlobalGC', 'twist_lgcLocalGC', 'twist_smrSmallRepeatPercentage']