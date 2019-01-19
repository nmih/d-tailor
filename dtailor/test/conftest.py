import pytest


_mfe_levels = {'1': (-9999, -6.3),
               '2': (-6.3, -4.5),
               '3': (-4.5, -3.1),
               '4': (-3.1, -1.8),
               '5': (-1.8, 0.0)}

_cai_levels = {'1': (0, 0.6937817343645736),
               '2': (0.6937817343645736, 0.7224804374510423),
               '3': (0.7224804374510423, 0.7419150216575061),
               '4': (0.7419150216575061, 0.7666072498092797),
               '5': (0.7666072498092797, 1)}

_targets = ['5.5.5', '4.4.4', '3.3.3', '2.2.2', '1.1.1']

_cai_table = {'aaa': 1.0,
 'aac': 0.6835958956781579,
 'aag': 0.8501819970486965,
 'aat': 1.0,
 'aca': 0.9698236035279295,
 'acc': 0.6997060058798824,
 'acg': 0.3369172616547669,
 'act': 1.0,
 'aga': 1.0,
 'agc': 0.33871546465733765,
 'agg': 0.4327671815167005,
 'agt': 0.6517581628991748,
 'ata': 0.44500644063390016,
 'atc': 0.560287878608823,
 'atg': 1.0,
 'att': 1.0,
 'caa': 1.0,
 'cac': 0.4904431664411366,
 'cag': 0.5018068436200831,
 'cat': 1.0,
 'cca': 1.0,
 'ccc': 0.27423029492297724,
 'ccg': 0.26806429362707185,
 'cct': 0.6827540079008424,
 'cga': 0.13356785576870225,
 'cgc': 0.059216914930650975,
 'cgg': 0.11079352864654904,
 'cgt': 0.2555012000952747,
 'cta': 0.45494302516671425,
 'ctc': 0.2547365847695881,
 'ctg': 0.26136440804574285,
 'ctt': 0.4427739070203317,
 'gaa': 1.0,
 'gac': 0.4646572115247604,
 'gag': 0.5341032773007619,
 'gat': 1.0,
 'gca': 1.0,
 'gcc': 0.4917465037926556,
 'gcg': 0.2468331796194169,
 'gct': 0.8223192606674321,
 'gga': 0.5128889736499158,
 'ggc': 0.3511013635930199,
 'ggg': 0.27271860398588726,
 'ggt': 1.0,
 'gta': 0.2980785532636338,
 'gtc': 0.45507205425261377,
 'gtg': 0.3633088443063012,
 'gtt': 1.0,
 'taa': 1.0,
 'tac': 0.6870831268818843,
 'tag': 0.6097455799913757,
 'tat': 1.0,
 'tca': 0.8974345174022247,
 'tcc': 0.6188733405095085,
 'tcg': 0.4076785073555795,
 'tct': 1.0,
 'tga': 0.6722725312634756,
 'tgc': 0.41045642382875985,
 'tgg': 1.0,
 'tgt': 1.0,
 'tta': 0.8421410042238792,
 'ttc': 0.649768556417164,
 'ttg': 1.0,
 'ttt': 1.0}


@pytest.fixture(scope='module')
def mfe_levels():
    return _mfe_levels


@pytest.fixture(scope='module')
def cai_levels():
    return _cai_levels


@pytest.fixture(scope='module')
def targets():
    return _targets


@pytest.fixture(scope='module')
def cai_table():
    return _cai_table

@pytest.fixture(scope='module')
def test_sequence():
    return 'ATGTCTGCAACTTCCGTCACTTTCCCAATTATCAACGAAACTTACCAACAGCCAACCGGGCTTTTCATCAACAATGAATTTGTTAGTGCAAAGTCAGGTAAGACTTTTGATGTTAACACTCCAATTGATGAGTCTCTCATTTGTAAAGTCCAACAGGCCGATGCTGAAGATGTTGAAATTGCCGTTCAAGCAGCATCTAAAGCTTACAAGACTTGGAGATTTACACCGCCAAATGAAAGAGGCAGATACTTGAACAAATTGGCCGATTTGATGGACGAAAAGAGAGACTTACTTGCCAAAATTGAATCCCTTGATAATGGTAAGGCCTTACATTGTGCAAAATTCGATGTCAATCTTGTCATTGAATATTTCAGATACTGTGCAGGTTACTGTGATAAAATCGATGGTAGAACAATTACAACCGATGTCGAACATTTTACCTACACTAGAAAGGAACCTTTAGGTGTCTGTGGTGCAATTACACCTTGGAACTTCCCATTGCTGATGTTTGCTTGGAAAATCGGCCCGGCTTTAGCAACCGGTAATACCATTATTTTGAAGCCTGCCAGTGCAACACCTCTATCAAACCTCTTTACTTGTACCTTGATCAAGGAGGCGGGCATTCCAGCCGGTGTTGTTAATGTTGTTCCAGGTTCCGGTAGAGGCTGTGGTAACTCCATTTTACAACATCCTAAAATTAAGAAGGTTGCGTTTACCGGATCTACAGAAGTTGGTAAAACTGTTATGAAGGAATGTGCTAATTCCATCAAAAAGGTTACTCTCGAATTGGGTGGTAAGTCTCCAAACATTGTTTTCAAAGACTGTAACGTTGAACAAACCATTCAAAATTTGATTACTGGTATTTTCTTCAATGGTGGTGAAGTCTGTTGTGCTGGTTCTAGAATTTACATTGAAGCAACCGATGAGAAATGGTATACTGAATTCTTGACCAAATTCAAGGAGACTGTTGAAAAATTAAAGATTGGTAACCCATTTGAAGAGGGTGTTTTCCAAGGTGCACAAACCACTCCAGATCAATTCCAAACTGTCTTGGACTACATCACCGCTGCTAACGAATCCAGCTTGAAACTATTAACTGGTGGTAAAAGAATTGGCAATAAGGGATACTTTGTTGAGCCAACTATCTTCTACGATGTTCCTCAAAATTCCAAGTTAACTCAAGAAGAAATCTTTGGTCCAGTTGCTGTTGTTTTACCTTTCAAGTCCACTGAAGAATTGATTGAAAAGGCAAATGATTCCGATTTTGGTTTAGGTTCCGGTATTCACACTGAAGATTTCAACAAGGCAATTTGGGTTTCCGAAAGGCTTGAAGCAGGTTCTGTTTGGATCAACACTTACAATGATTTCCACCCAGCTGCTCCATTCGGTGGTTACAAGGAATCCGGTATTGGCAGAGAAATGGGTATTGAAGCTTTCGACAACTATACTCAAACCAAGTTAGTTAGAGCTAGAGTTAACAAGCCAGCTTTTTAG'


@pytest.fixture(scope='module')
def design_params(test_sequence):
    return {
        "mfe": {
            'feattype'      : 'StructureRNAFoldMFE',
            'type'          : 'REAL',
            'mutable_region': (0, 32),
            'thresholds'    : _mfe_levels},
        "ramp"           : {
            'feattype'      : 'CAI',
            'type'          : 'REAL',
            'mutable_region': (0, 32),
            'thresholds'    : _cai_levels},
        "rest"           : {
            'feattype'      : 'CAI',
            'type'          : 'REAL',
            'mutable_region': (33, len(test_sequence)),
            'thresholds'    : _cai_levels}
    }