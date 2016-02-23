

ENCODE_LAB_LST = ['Bradley Bernstein, Broad', 'Michael Snyder, Stanford', 'Peggy Farnham, USC',
                  'Kevin White, UChicago', 'John Stamatoyannopoulos, UW', 'Ross Hardison, PennState',
                  'Richard Myers, HAIB', 'Gregory Crawford, Duke', 'Bing Ren, UCSD',
                  'Vishwanath Iyer, UTA', 'Kevin Struhl, HMS', 'Sherman Weissman, Yale']

ENCODE_LAB_MAP = {}
for idx, lab in enumerate(ENCODE_LAB_LST):
    ENCODE_LAB_MAP[lab] = 'L' + str(idx).zfill(2)


ENCODE_BIOSAMPLES_MAP = {'HepG2': 'HepG2', 'K562': 'K562', 'GM12878': 'GM12878',
                         'H1-hESC': 'H1hESC', 'MEL cell line': 'MEL',
                         'CH12.LX': 'CH12', 'ES-Bruce4': 'ESB4', 'ES-E14': 'ESE14'}
