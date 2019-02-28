import pytest
import os.path as op
from dtailor.RunningExamples.Designer.GenericDesigner import GenericDesigner
from dtailor.DesignOfExperiments.Design import CustomDesign, RandomSampling


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


@pytest.fixture(scope='module')
def random_designer(design_params, targets, pk_cai_table, test_sequence):
    des = RandomSampling(featuresObj=design_params, sample_size=10)

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


def test_run2(random_designer):
    random_designer.run(selection='neutral')


def highlocalgc(design_params, targets, pk_cai_table, test_sequence):
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