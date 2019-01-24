import pytest
import os.path as op
from dtailor.RunningExamples.Designer.GenericDesigner import GenericDesigner
from dtailor.DesignOfExperiments.Design import CustomDesign


@pytest.fixture(scope='module')
def gen_designer(design_params, targets, pk_cai_table, test_sequence):
    des = CustomDesign(featuresObj=design_params,
                       targets=targets)

    name = 'test_name'
    seed = test_sequence
    root_dir = '/tmp/'
    gendes = GenericDesigner(
            name=name,
            seed=seed,
            design=des,
            cai_table=pk_cai_table,
            root_dir=root_dir,
            mutable_region=None,
            cds_region=None,
            keep_aa=True,
            createDB=True)

    return gendes


def test_run(gen_designer):
    gen_designer.run()