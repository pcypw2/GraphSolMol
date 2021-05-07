
##---------------------------------------------------------------------------------------------##
## Code was originally curated by Phyo Phyo Kyaw Zin, the link to the blog containing this is: ##
##                    https://drzinph.com/category/cheminformatics/                            ##
##---------------------------------------------------------------------------------------------##

import pandas as pd
from molvs import standardize_smiles
from Descriptors import ECFP6

def main():
    filename = 'AqSolDB_C.csv'
    df = pd.read_csv(filename)
    smiles = [standardize_smiles(i) for i in df['SMILES'].values]

    graph_descriptor = ECFP6(smiles)
    graph_descriptor.compute_ECFP6(
        filename)


if __name__ == '__main__':
    main()

