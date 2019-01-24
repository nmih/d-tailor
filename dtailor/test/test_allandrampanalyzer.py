import pytest
import tempfile
from dtailor.RunningExamples.Analyzer.AllAndRampAnalyzer import AllAndRampAnalyzer


@pytest.fixture(scope='module')
def allandramp_analyzer(pk_cds_fasta, pk_cai_table):
    dirpath = tempfile.mkdtemp()
    return AllAndRampAnalyzer(input_file=pk_cds_fasta,
                              input_type="FASTA",
                              cai_table=pk_cai_table,
                              ramp_from_to=(0, 39),
                              root_dir=dirpath)


def test_run(allandramp_analyzer):
    analyzed = allandramp_analyzer.run()
