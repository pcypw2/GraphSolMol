"""
ESOL Model predicts solubility directly from molecular structure with the use of logP, molecular weight (MW),
rotatable bonds (RB) and aromatic proportion (AP).
Delaney JS. J. Chem. Inf. Comput. Sci., 44(3), 1000–1005 (2004).

Method was inspired by The Data Professor https://github.com/dataprofessor
"""
import pandas as pd
import numpy as np
from rdkit.Chem import Descriptors
from rdkit import Chem
from math import sqrt

sol = pd.read_csv('AqSolDB_C.csv')
mol_list = []
for element in sol.SMILES:
    mol = Chem.MolFromSmiles(element)
    mol_list.append(mol)


def generate(smiles, verbose=False):
    moldata = []
    for elem in smiles:
        mol = Chem.MolFromSmiles(elem)
        moldata.append(mol)

    baseData = np.arange(1, 1)
    i = 0
    for mol in moldata:

        desc_MolLogP = Descriptors.MolLogP(mol)
        desc_MolWt = Descriptors.ExactMolWt(mol)
        desc_NumRotatableBonds = Descriptors.NumRotatableBonds(mol)

        row = np.array([desc_MolLogP,
                        desc_MolWt,
                        desc_NumRotatableBonds])

        if i == 0:
            baseData = row
        else:
            baseData = np.vstack([baseData, row])
        i = i + 1

    columnNames = ["MolLogP", "MolWt", "NumRotatableBonds"]
    descriptors = pd.DataFrame(data=baseData, columns=columnNames)

    return descriptors


df = generate(sol.SMILES)


# AP = The proportion of heavy atoms in the molecule that are in an aromatic ring

# Aromatic atoms
def AromaticAtoms(m):
    aromatic_atoms = [m.GetAtomWithIdx(i).GetIsAromatic() for i in range(m.GetNumAtoms())]
    aa_count = []
    for i in aromatic_atoms:
        if i == True:
            aa_count.append(1)
    sum_aa_count = sum(aa_count)
    return sum_aa_count


desc_AromaticAtoms = [AromaticAtoms(element) for element in mol_list]

# Number of heavy atoms
desc_HeavyAtomCount = [Descriptors.HeavyAtomCount(element) for element in mol_list]

# AP
desc_AromaticProportion = [AromaticAtoms(element) / Descriptors.HeavyAtomCount(element) for element in mol_list]

# Combining the dataframes
df_desc_AromaticProportion = pd.DataFrame(desc_AromaticAtoms, columns=['AromaticProportion'])
X = pd.concat([df, df_desc_AromaticProportion], axis=1)
Y = sol['Solubility']
descriptors = pd.concat([sol, X], axis=1)
descriptors.to_csv('AqSolDB_C_ESOL_Desc')

#Standard Scaler
from sklearn.preprocessing import StandardScaler
sc_X = StandardScaler()
sc_y = StandardScaler()
X = sc_X.fit_transform(X)


# Datasplit
from sklearn.model_selection import train_test_split

X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.2)

from sklearn import linear_model
from sklearn.metrics import mean_squared_error, r2_score

model = linear_model.LinearRegression()
model.fit(X_train, Y_train)

# Predicting logs
Y_pred_train = model.predict(X_train)
print('Training set coefficients:', model.coef_)
print('Training set intercept:', model.intercept_)
print('Training set root mean squared error (RMSE): %.2f'
      % sqrt(mean_squared_error(Y_train, Y_pred_train)))
print('Training set coefficient of determination (R^2): %.2f'
      % r2_score(Y_train, Y_pred_train))

Y_pred_test = model.predict(X_test)
print('Test set coefficients:', model.coef_)
print('Test set intercept:', model.intercept_)
print('Test set root mean squared error (RMSE): %.2f'
      % sqrt(mean_squared_error(Y_test, Y_pred_test)))
print('Test set coefficient of determination (R^2): %.2f'
      % r2_score(Y_test, Y_pred_test))

# Deriving the linear regression equation
from sklearn import linear_model
from sklearn.metrics import mean_squared_error, r2_score

model = linear_model.LinearRegression()
model.fit(X_train, Y_train)

yintercept = '%.2f' % model.intercept_
LogP = '%.2f LogP' % model.coef_[0]
MW = '%.4f MW' % model.coef_[1]
RB = '%.4f RB' % model.coef_[2]
AP = '%.2f AP' % model.coef_[3]
print('LogS = ' +
      ' ' +
      yintercept +
      ' ' +
      LogP +
      ' ' +
      MW +
      ' ' +
      RB +
      ' ' +
      AP)

# Full dataset
full = linear_model.LinearRegression()
full.fit(X, Y)

full_pred = model.predict(X)
print('Full dataset coefficients:', full.coef_)
print('Full dataset intercept:', full.intercept_)
print('Full dataset mean squared error (RMSE): %.2f'
      % sqrt(mean_squared_error(Y, full_pred)))
print('Full dataset coefficient of determination (R^2): %.2f'
      % r2_score(Y, full_pred))

full_yintercept = '%.2f' % full.intercept_
full_LogP = '%.2f LogP' % full.coef_[0]
full_MW = '%.4f MW' % full.coef_[1]
full_RB = '+ %.4f RB' % full.coef_[2]
full_AP = '%.2f AP' % full.coef_[3]
print('LogS = ' +
      ' ' +
      full_yintercept +
      ' ' +
      full_LogP +
      ' ' +
      full_MW +
      ' ' +
      full_RB +
      ' ' +
      full_AP)

import matplotlib.pyplot as plt

plt.figure(figsize=(11, 5))

# 1 row, 2 column, plot 1
plt.subplot(1, 2, 1)
plt.scatter(x=Y_train, y=Y_pred_train, c="#7CAE00", alpha=0.3)

z = np.polyfit(Y_train, Y_pred_train, 1)
p = np.poly1d(z)
plt.plot(Y_test, p(Y_test), "#F8766D")

plt.ylabel('Predicted LogS')
plt.xlabel('Experimental LogS')

# 1 row, 2 column, plot 2
plt.subplot(1, 2, 2)
plt.scatter(x=Y_test, y=Y_pred_test, c="#619CFF", alpha=0.3)

z = np.polyfit(Y_test, Y_pred_test, 1)
p = np.poly1d(z)
plt.plot(Y_test, p(Y_test), "#F8766D")

plt.xlabel('Experimental LogS')
plt.ylabel('Predicted LogS')

plt.savefig('plot_horizontal_logS.png')
print(plt.show())


df_predicted = pd.DataFrame(full_pred, columns=['Prediction'])
predictedESOL = pd.concat([sol, df_predicted], axis=1)
predictedESOL.to_csv('Predicted ESOL validation.csv')

