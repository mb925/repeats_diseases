import os

absolute = os.path.abspath(os.getcwd())
absolute = absolute + '/../'

data = {
    "trp_diseases": absolute + 'data/trp_diseases',
    "figures": absolute + 'data/figures',
    "backup": absolute + 'data/backup',
    "interactors": absolute + 'data/interactors',
    "PTMs": absolute + 'data/PTMs',
    "modified_residues": absolute + 'data/modified_residues',
    "repeats": absolute + 'data/repeats',
}
