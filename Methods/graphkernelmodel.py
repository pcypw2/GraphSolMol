from rdkit import Chem
from math import sqrt
import grakel
import numpy as np
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import mean_squared_error, r2_score
import pandas as pd
from sklearn.svm import SVR
from sklearn.model_selection import train_test_split
soly = pd.read_csv("AqSolDB_C.csv")
val = pd.read_csv("validation dataset.csv")


def molg_from_smi(smiles):
    mol = Chem.MolFromSmiles(smiles)
    bond_idx = [(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()) for bond in mol.GetBonds()] + \
               [(bond.GetEndAtomIdx(), bond.GetBeginAtomIdx()) for bond in mol.GetBonds()]
    atom_with_idx = {i:atom.GetSymbol() for i, atom in enumerate(mol.GetAtoms())}
    bond_with_idx = {(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()): bond.GetBondTypeAsDouble() for bond in
                     mol.GetBonds()}
    bond_with_idx.update({(bond.GetEndAtomIdx(), bond.GetBeginAtomIdx()): bond.GetBondTypeAsDouble() for bond in
                          mol.GetBonds()})

    adj_m = Chem.GetAdjacencyMatrix(mol).tolist()


    return grakel.Graph(adj_m, edge_labels=bond_with_idx, node_labels=atom_with_idx, graph_format='adjacency')
mols = soly.SMILES
mols_val = val.SMILES

X = [molg_from_smi(mol) for mol in mols]


X_val = [molg_from_smi(mol) for mol in mols_val]
Y = soly.Solubility
y_val = val.Solubility

#Splitting the data
X_train, X_test, y_train, y_test = train_test_split(X,Y, test_size=0.2, random_state=42)

#Choosing the kernel method - grakel has a wide list of kernels to choose from (most
#have been covered in this project, see their github for instructions on how to use them

gk = grakel.WeisfeilerLehmanOptimalAssignment(normalize=True, n_iter=5)


#Calculate the normalized kernel matrices
# Sometimes the kernel matrix would compute NaN, np.nan_to_num would filter this problem out
K_train = gk.fit_transform(X_train)
K_trainr = np.nan_to_num(K_train)

K_test = gk.transform(X_test)
K_testr = np.nan_to_num(K_test)

K_val = gk.transform(X_val)
K_valr = np.nan_to_num(K_val)





#Add SVR method

gsc = GridSearchCV(
    estimator=SVR(kernel='precomputed'),
    param_grid={
       'C': [2**0, 2**1, 2**2, 2**3, 2**4, 2**5, 2**6, 2**7, 2**8],
        'gamma':[2**-15, 2**-13, 2**-11, 2**-9, 2**-7, 2**-5, 2**-3, 2**-1, 2**1, 2**0, 2**1,
                 2**3, 2**5]
      },
    cv=5, scoring='neg_root_mean_squared_error', verbose=2, n_jobs=-1)
sc_arr_y = y_train.ravel()
grid_result = gsc.fit(K_train, y_train)
best_params = grid_result.best_params_
print(best_params)

best_svr = SVR(kernel='precomputed', C= best_params["C"], gamma= best_params['gamma'])

best_svr.fit(K_train, y_train)

#Make predictions
y_pred_train = best_svr.predict(K_trainr)
y_pred_test = best_svr.predict(K_testr)
y_pred_val = best_svr.predict(K_valr)



#Print model performance results
print('Training set RMSE: %.2f' % sqrt(mean_squared_error(y_train,y_pred_train)))
print('Training set R2: %.2f' % r2_score(y_train,y_pred_train))

print('Test set RMSE: %.2f' % sqrt(mean_squared_error(y_test,y_pred_test)))
print('Test set R2: %.2f' % r2_score(y_test,y_pred_test))

print('Validation set RMSE: %.2f' % sqrt(mean_squared_error(y_val,y_pred_val)))
print('Validation set R2: %.2f' % r2_score(y_val,y_pred_val))

import matplotlib.pyplot as plt

plt.figure(figsize=(15, 5))
plt.suptitle("Weisfieler-Lehan Optimal Assignment Graph Kernel Predicted vs Experimental", fontsize=16)
# 1 row, 2 column, plot 1
ax1 = plt.subplot(1, 3, 1)
ax1.scatter(x=y_train, y=y_pred_train, c="#1f77b4", alpha=0.3)
ax1.plot(np.unique(y_train), np.poly1d(np.polyfit(y_train, y_pred_train, 1))(np.unique(y_train)), "#F8766D")
plt.ylabel('Predicted LogS', fontsize=12)
plt.xlabel('Experimental LogS', fontsize=12)


# 1 row, 2 column, plot 2
ax2 = plt.subplot(1, 3, 2)
ax2.scatter(x=y_test, y=y_pred_test, c="#619CFF", alpha=0.3)
#ax2.title.set_text('Test Set')
plt.ylabel('Predicted LogS', fontsize=12)
plt.xlabel('Experimental LogS', fontsize=12)
ax2.plot(np.unique(y_test), np.poly1d(np.polyfit(y_test, y_pred_test, 1))(np.unique(y_test)), "#F8766D")

# 2 row, 1 column, plot 3
ax3 = plt.subplot(1, 3, 3)
ax3.scatter(x=y_val, y=y_pred_val, c='#7CAE00', alpha=0.3)

ax3.plot(np.unique(y_val), np.poly1d(np.polyfit(y_val, y_pred_val, 1))(np.unique(y_val)), "#F8766D")

plt.ylabel('Predicted LogS', fontsize=12)
plt.xlabel('Experimental LogS', fontsize=12)
plt.axis('equal')

plt.subplots_adjust(left=0.125,
                    bottom=0.1,
                    right=0.9,
                    top=0.9,
                    wspace=0.2,
                    hspace=0.35)
plt.savefig('WLOAExpvsPred.png')
print(plt.show())
