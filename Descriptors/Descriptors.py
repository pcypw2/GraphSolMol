
##---------------------------------------------------------------------------------------------##
## Code was originally curated by Phyo Phyo Kyaw Zin, the link to the blog containing this is: ##
##                    https://drzinph.com/category/cheminformatics/                            ##
##---------------------------------------------------------------------------------------------##

import pandas as pd
import numpy as np
from rdkit.Chem import Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.Chem import MACCSkeys
from rdkit import Chem


class MACCS:
    def __init__(self, smiles):
        self.smiles = smiles
        self.mols = [Chem.MolFromSmiles(i) for i in smiles]

    def compute_MACCS(self, name):
        MACCS_list = []
        header = ['bit' + str(i) for i in range(167)]
        for i in range(len(self.mols)):
            ds = list(MACCSkeys.GenMACCSKeys(self.mols[i]).ToBitString())
            MACCS_list.append(ds)
        df = pd.DataFrame(MACCS_list,columns=header)
        df.insert(loc=0, column='smiles', value=self.smiles)
        df.to_csv(name[:-4]+'_MACCS.csv', index=False)

class RDKit_2D:
    def __init__(self, smiles):
        self.mols = [Chem.MolFromSmiles(i) for i in smiles]
        self.smiles = smiles

    def compute_2Drdkit(self, name, smiles):
        rdkit_2d_desc = []
        calc = MoleculeDescriptors.MolecularDescriptorCalculator([x[0] for x in Descriptors._descList])
        header = calc.GetDescriptorNames()
        for i in range(len(self.mols)):
            ds = calc.CalcDescriptors(self.mols[i])
            rdkit_2d_desc.append(ds)
        df = pd.DataFrame(rdkit_2d_desc, columns=header)
        df.insert(loc=0, column='smiles', value=self.smiles)
        df.to_csv(name[:-4] + '_RDKit_2D.csv', index=False)

class ECFP6:
    def __init__(self, smiles):
        self.mols = [Chem.MolFromSmiles(i) for i in smiles]
        self.smiles = smiles

    def mol2fp(self, mol, radius = 3):
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius = radius)
        array = np.zeros((1,))
        DataStructs.ConvertToNumpyArray(fp, array)
        return array

    def compute_ECFP6(self, name):
        bit_headers = ['bit' + str(i) for i in range(2048)]
        arr = np.empty((0,2048), int).astype(int)
        for i in self.mols:
            fp = self.mol2fp(i)
            arr = np.vstack((arr, fp))
        df_ecfp6 = pd.DataFrame(np.asarray(arr).astype(int),columns=bit_headers)
        df_ecfp6.insert(loc=0, column='smiles', value=self.smiles)
        df_ecfp6.to_csv(name[:-4]+'_ECFP6.csv', index=False)

class ECFP4:
    def __init__(self, smiles):
        self.mols = [Chem.MolFromSmiles(i) for i in smiles]
        self.smiles = smiles

    def mol2fp(self, mol, radius = 2):
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius = radius)
        array = np.zeros((1,))
        DataStructs.ConvertToNumpyArray(fp, array)
        return array

    def compute_ECFP4(self, name):
        bit_headers = ['bit' + str(i) for i in range(2048)]
        arr = np.empty((0,2048), int).astype(int)
        for i in self.mols:
            fp = self.mol2fp(i)
            arr = np.vstack((arr, fp))
        df_ecfp4 = pd.DataFrame(np.asarray(arr).astype(int),columns=bit_headers)
        df_ecfp4.insert(loc=0, column='smiles', value=self.smiles)
        df_ecfp4.to_csv(name[:-4]+'_ECFP4.csv', index=False)

class ECFP2:
    def __init__(self, smiles):
        self.mols = [Chem.MolFromSmiles(i) for i in smiles]
        self.smiles = smiles

    def mol2fp(self, mol, radius = 1):
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius = radius)
        array = np.zeros((1,))
        DataStructs.ConvertToNumpyArray(fp, array)
        return array

    def compute_ECFP2(self, name):
        bit_headers = ['bit' + str(i) for i in range(2048)]
        arr = np.empty((0,2048), int).astype(int)
        for i in self.mols:
            fp = self.mol2fp(i)
            arr = np.vstack((arr, fp))
        df_ecfp2 = pd.DataFrame(np.asarray(arr).astype(int),columns=bit_headers)
        df_ecfp2.insert(loc=0, column='smiles', value=self.smiles)
        df_ecfp2.to_csv(name[:-4]+'_ECFP2.csv', index=False)

class test_split:
    def __init__(self, smiles):
        self.mols = [Chem.MolFromSmiles(i) for i in smiles]
        self.smiles = smiles

    def compute_2Drdkit(self, name):
        rdkit_2d_desc = []
        calc = MoleculeDescriptors.MolecularDescriptorCalculator([x[0] for x in Descriptors._descList])
        header = calc.GetDescriptorNames()
        for i in range(len(self.mols)):
            ds = calc.CalcDescriptors(self.mols[i])
            rdkit_2d_desc.append(ds)
        df = pd.DataFrame(rdkit_2d_desc, columns=header)
        df.insert(loc=0, column='smiles', value=self.smiles)
        df.to_csv(name[:-4] + '_RDKit_2D.csv', index=False)


