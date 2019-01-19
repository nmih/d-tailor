import pytest
from dtailor.DesignOfExperiments.Design import Design, CustomDesign


def test_design_init(design_params, mfe_levels, cai_levels):
    des = Design(featuresObj=design_params)
    assert sorted(des.featuresList) == sorted(['mfe', 'ramp', 'rest'])
    assert des.n_features == 3
    assert des.features == design_params
    assert des.thresholds == {
        'mfe' : mfe_levels,
        'ramp': cai_levels,
        'rest': cai_levels
    }


def test_custom_design_init(design_params, targets):
    des = CustomDesign(featuresObj=design_params,
                       targets=targets)
    assert des.listDesigns == targets
    assert des.nDesigns == len(targets)