import pytest
import tempfile
from dtailor.RunningExamples.Analyzer.AllAndRampAnalyzer import AllAndRampAnalyzer


@pytest.fixture(scope='module')
def allandramp_analyzer(pk_cds_fasta, pk_cai_table):
    dirpath = tempfile.mkdtemp()
    return AllAndRampAnalyzer(input_file=pk_cds_fasta,
                              input_type="FASTA",
                              cai_table=pk_cai_table,
                              ramp=True,
                              ramp_from_to=(0, 33),
                              root_dir=dirpath)


def test_run(allandramp_analyzer, pk_cai_mfe_precalc_1_33):
    analyzed = allandramp_analyzer.run()
    assert analyzed == pk_cai_mfe_precalc_1_33