import pytest
from pkg_resources import resource_filename
from collections import OrderedDict
import json

###################
# These levels and CAI table were generated using the 181218-pk_cds_cleaned_ids_TRUNCATED.fna as the CDS base file,
# so that means the table was made directly from that, CAI levels were generated using that table + genes in the file,
# and MFE was generated from genes in the file
###################


_mfe_levels = {'3': (-3.2, -2.4399999999999995),
               '4': (-2.4399999999999995, -1.3599999999999988),
               '5': (-1.3599999999999988, 0.0),
               '2': (-4.48, -3.2),
               '1': (-9999, -4.48)}


_cai_levels = {'3': (0.7244389601788824, 0.7466550646222192),
               '4': (0.7466550646222192, 0.7663049226209018),
               '5': (0.7663049226209018, 1),
               '2': (0.6906052493341361, 0.7244389601788824),
               '1': (0, 0.6906052493341361)}

_targets = ['5.5.5', '4.4.4', '3.3.3', '2.2.2', '1.1.1']

_cai_table = {'aaa': 1.0,
 'aac': 0.610738255033557,
 'aag': 0.8794520547945205,
 'aat': 1.0,
 'aca': 1.0,
 'acc': 0.7087378640776699,
 'acg': 0.3576051779935275,
 'act': 0.9045307443365695,
 'aga': 1.0,
 'agc': 0.38326585695006743,
 'agg': 0.5191256830601093,
 'agt': 0.6707152496626181,
 'ata': 0.49450549450549447,
 'atc': 0.5173288250211326,
 'atg': 1.0,
 'att': 1.0,
 'caa': 1.0,
 'cac': 0.44660194174757284,
 'cag': 0.5353535353535354,
 'cat': 1.0,
 'cca': 1.0,
 'ccc': 0.27874015748031494,
 'ccg': 0.3291338582677165,
 'cct': 0.6740157480314961,
 'cga': 0.1721311475409836,
 'cgc': 0.04644808743169399,
 'cgg': 0.1284153005464481,
 'cgt': 0.26366120218579236,
 'cta': 0.47551020408163264,
 'ctc': 0.2326530612244898,
 'ctg': 0.286734693877551,
 'ctt': 0.4091836734693878,
 'gaa': 1.0,
 'gac': 0.45515846257585973,
 'gag': 0.540119760479042,
 'gat': 1.0,
 'gca': 1.0,
 'gcc': 0.5658914728682171,
 'gcg': 0.262015503875969,
 'gct': 0.7829457364341085,
 'gga': 0.5637583892617449,
 'ggc': 0.39999999999999997,
 'ggg': 0.34630872483221475,
 'ggt': 1.0,
 'gta': 0.32799145299145294,
 'gtc': 0.42841880341880334,
 'gtg': 0.41239316239316237,
 'gtt': 1.0,
 'taa': 1.0,
 'tac': 0.6517482517482518,
 'tag': 0.6071428571428571,
 'tat': 1.0,
 'tca': 0.932523616734143,
 'tcc': 0.6167341430499325,
 'tcg': 0.45209176788124156,
 'tct': 1.0,
 'tga': 0.9642857142857142,
 'tgc': 0.393103448275862,
 'tgg': 1.0,
 'tgt': 1.0,
 'tta': 0.8571428571428572,
 'ttc': 0.6968152866242038,
 'ttg': 1.0,
 'ttt': 1.0}


@pytest.fixture(scope='module')
def pk_mfe_levels():
    return _mfe_levels


@pytest.fixture(scope='module')
def pk_cai_levels():
    return _cai_levels


@pytest.fixture(scope='module')
def targets():
    return _targets


@pytest.fixture(scope='module')
def pk_cai_table():
    return _cai_table


@pytest.fixture(scope='module')
def test_sequence():
    return 'ATGTCTGCAACTTCCGTCACTTTCCCAATTATCAACGAAACTTACCAACAGCCAACCGGGCTTTTCATCAACAATGAATTTGTTAGTGCAAAGTCAGGTAAGACTTTTGATGTTAACACTCCAATTGATGAGTCTCTCATTTGTAAAGTCCAACAGGCCGATGCTGAAGATGTTGAAATTGCCGTTCAAGCAGCATCTAAAGCTTACAAGACTTGGAGATTTACACCGCCAAATGAAAGAGGCAGATACTTGAACAAATTGGCCGATTTGATGGACGAAAAGAGAGACTTACTTGCCAAAATTGAATCCCTTGATAATGGTAAGGCCTTACATTGTGCAAAATTCGATGTCAATCTTGTCATTGAATATTTCAGATACTGTGCAGGTTACTGTGATAAAATCGATGGTAGAACAATTACAACCGATGTCGAACATTTTACCTACACTAGAAAGGAACCTTTAGGTGTCTGTGGTGCAATTACACCTTGGAACTTCCCATTGCTGATGTTTGCTTGGAAAATCGGCCCGGCTTTAGCAACCGGTAATACCATTATTTTGAAGCCTGCCAGTGCAACACCTCTATCAAACCTCTTTACTTGTACCTTGATCAAGGAGGCGGGCATTCCAGCCGGTGTTGTTAATGTTGTTCCAGGTTCCGGTAGAGGCTGTGGTAACTCCATTTTACAACATCCTAAAATTAAGAAGGTTGCGTTTACCGGATCTACAGAAGTTGGTAAAACTGTTATGAAGGAATGTGCTAATTCCATCAAAAAGGTTACTCTCGAATTGGGTGGTAAGTCTCCAAACATTGTTTTCAAAGACTGTAACGTTGAACAAACCATTCAAAATTTGATTACTGGTATTTTCTTCAATGGTGGTGAAGTCTGTTGTGCTGGTTCTAGAATTTACATTGAAGCAACCGATGAGAAATGGTATACTGAATTCTTGACCAAATTCAAGGAGACTGTTGAAAAATTAAAGATTGGTAACCCATTTGAAGAGGGTGTTTTCCAAGGTGCACAAACCACTCCAGATCAATTCCAAACTGTCTTGGACTACATCACCGCTGCTAACGAATCCAGCTTGAAACTATTAACTGGTGGTAAAAGAATTGGCAATAAGGGATACTTTGTTGAGCCAACTATCTTCTACGATGTTCCTCAAAATTCCAAGTTAACTCAAGAAGAAATCTTTGGTCCAGTTGCTGTTGTTTTACCTTTCAAGTCCACTGAAGAATTGATTGAAAAGGCAAATGATTCCGATTTTGGTTTAGGTTCCGGTATTCACACTGAAGATTTCAACAAGGCAATTTGGGTTTCCGAAAGGCTTGAAGCAGGTTCTGTTTGGATCAACACTTACAATGATTTCCACCCAGCTGCTCCATTCGGTGGTTACAAGGAATCCGGTATTGGCAGAGAAATGGGTATTGAAGCTTTCGACAACTATACTCAAACCAAGTTAGTTAGAGCTAGAGTTAACAAGCCAGCTTTTTAG'


@pytest.fixture(scope='module')
def design_params(test_sequence):
    dp = OrderedDict()
    dp['mfe'] = {
            'feattype'      : 'StructureRNAFoldMFE',
            'type'          : 'REAL',
            'mutable_region': (0, 33),
            'thresholds'    : _mfe_levels}
    dp['ramp'] = {
            'feattype'      : 'CAI',
            'type'          : 'REAL',
            'mutable_region': (0, 33),
            'thresholds'    : _cai_levels}
    dp['rest'] = {
            'feattype'      : 'CAI',
            'type'          : 'REAL',
            'mutable_region': (33, len(test_sequence)),
            'thresholds'    : _cai_levels}
    return dp


@pytest.fixture(scope='module')
def pk_cds_fasta():
    return resource_filename(__name__, 'test_files/181218-pk_cds_cleaned_ids_TRUNCATED.fna')


@pytest.fixture(scope='module')
def pk_cai_mfe_precalc_1_33():
    with open(resource_filename(__name__, 'test_files/181218-pk_cds_cleaned_ids_TRUNCATED_CAI_MFE_ramp1-33.json'), 'r') as f:
        return json.load(f)