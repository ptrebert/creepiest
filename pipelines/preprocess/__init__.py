

ENCODE_LAB_LST = ['Bradley Bernstein, Broad', 'Michael Snyder, Stanford', 'Peggy Farnham, USC',
                  'Kevin White, UChicago', 'John Stamatoyannopoulos, UW', 'Ross Hardison, PennState',
                  'Richard Myers, HAIB', 'Gregory Crawford, Duke', 'Bing Ren, UCSD',
                  'Vishwanath Iyer, UTA', 'Kevin Struhl, HMS', 'Sherman Weissman, Yale',
                  'Barbara Wold, Caltech', 'Thomas Gingeras, CSHL', 'Sheng Zhong UCSD', 'EWHA']

ENCODE_LAB_MAP = {}
for idx, lab in enumerate(ENCODE_LAB_LST):
    ENCODE_LAB_MAP[lab] = 'L' + str(idx).zfill(2)


ENCODE_BIOSAMPLES_MAP = {'HepG2': 'HepG2', 'K562': 'K562', 'GM12878': 'GM12878',
                         'H1-hESC': 'H1hESC', 'MEL cell line': 'MEL',
                         'CH12.LX': 'CH12', 'ES-Bruce4': 'ESB4', 'ES-E14': 'ESE14'}

GENCODE_MAP = {'ENCFF000FVB': 'v7', 'ENCFF000FVR': 'v7',
               'ENCFF000HIY': 'v7', 'ENCFF000EWE': 'v7', 'ENCFF000EYY': 'v7',
               'ENCFF000FEM': 'v7', 'ENCFF000HEY': 'v7', 'ENCFF000FZZ': 'v7',
               'ENCFF000FCA': 'v7', 'ENCFF000CZD': 'v3c', 'ENCFF000CZB': 'v3c',
               'ENCFF000CZF': 'v3c', 'ENCFF000CZI': 'v3c', 'ENCFF000DYF': 'v3c',
               'ENCFF000DYC': 'v3c', 'ENCFF000DHQ': 'v3c', 'ENCFF000DHS': 'v3c',
               'ENCFF000DHU': 'v3c', 'ENCFF000DHW': 'v3c', 'ENCFF000DQM': 'v3c',
               'ENCFF000DQP': 'v3c'}
