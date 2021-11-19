import os

absolute = os.path.abspath(os.getcwd())
absolute = absolute + '/../'

data = {
    "trp_diseases": absolute + 'data/trp_diseases',
    "figures": absolute + 'data/figures'
}
