import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from math import sqrt
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import GridSearchCV
from sklearn.svm import SVR
from warnings import simplefilter
from sklearn.exceptions import ConvergenceWarning

#Importing AqSol and validation dataset

df = pd.read_csv("AqSolDB_C_ECFP6.csv")
val = pd.read_csv("validation dataset_ECFP6.csv")
X = df.iloc[:, 3:]
y = df.iloc[:, 2]
val_X = val.iloc[:, 3:]
val_y = val.iloc[:, 2]

#Solublity values required as a 2D array
arr_y = np.array(y)
arr_val_y = np.array(val_y)
narr_y = arr_y.reshape(-1, 1)
narr_y_val = arr_val_y.reshape(-1,1)

#Preprocessing the data using standard scaler

sc_X = StandardScaler()
sc_y = StandardScaler()
X_sc = sc_X.fit_transform(X)
y_sc = sc_y.fit_transform(narr_y)


X_val = sc_X.fit_transform(val_X)
y_val = sc_y.fit_transform(narr_y_val)


#Splitting the data into test and train
X_train, X_test, y_train, y_test = train_test_split(X_sc, y_sc, train_size= 0.8,test_size=0.2, random_state=42)

# Hyperparameter optimisation
simplefilter("ignore", category= ConvergenceWarning)
ConvergenceWarning('ignore')
gsc = GridSearchCV(
    estimator=SVR(kernel='rbf'),
    param_grid={
        'C': [2**-5, 2**-4, 2**-3, 2**-2, 2**-1, 2**0, 2**1, 2**2, 2**3, 2**4, 2**5, 2**6],
        'gamma': [2**-15, 2**-13, 2**-11, 2**-9, 2**-7, 2**-5, 2**-3, 2**-1, 2**1, 2**3,
                  2**5]
        },
    cv=5, scoring='neg_root_mean_squared_error', verbose=2, n_jobs=1)
sc_arr_y = y_train.ravel()
grid_result = gsc.fit(X_train, sc_arr_y)
best_params = grid_result.best_params_

best_svr = SVR(kernel='rbf',  C= best_params['C'], gamma=best_params['gamma'])
print(best_svr)


from sklearn.model_selection import cross_validate

rbf_svr = SVR(kernel='rbf', C= best_params['C'], gamma=best_params['gamma'], verbose=0)
scores = cross_validate(rbf_svr, X_train, sc_arr_y, scoring='neg_mean_squared_error', cv=5)
print(scores)

#Fitting the SVR model
rbf_svr.fit(X_train, sc_arr_y)

#Predicting solubility and calculate performance
y_pred_train = rbf_svr.predict(X_train)
y_pred_test = rbf_svr.predict((X_test))
y_pred_val = rbf_svr.predict((X_val))

print('Training set root mean squared error (RMSE): %.2f'% sqrt(mean_squared_error(y_train, y_pred_train)))
print('Training set coefficient of determination (R^2): %.2f'% r2_score(y_train, y_pred_train))
print('Test set root mean squared error (RMSE): %.2f' % sqrt(mean_squared_error(y_test, y_pred_test)))
print('Test set coefficient of determination (R^2): %.2f'% r2_score(y_test, y_pred_test))
print('Validation  root mean squared error (RMSE): %.2f' % sqrt(mean_squared_error(y_val, y_pred_val)))
print('Validation set coefficient of determination (R^2): %.2f'% r2_score(y_val, y_pred_val))

#Convert solubility data from 2D array to 1D for visualisation
y_trainflat = y_train.flatten()
y_testflat = y_test.flatten()
y_valflat = y_val.flatten()

#Visualisation
plt.figure(figsize=(15, 5))
plt.suptitle("ECFP RBF Kernel Predicted vs Experimental", fontsize=16)

# 1 row, 2 column, plot 1
ax1 = plt.subplot(1, 3, 1)
ax1.scatter(x=y_trainflat, y=y_pred_train, c="#1f77b4", alpha=0.3)
ax1.plot(np.unique(y_trainflat), np.poly1d(np.polyfit(y_trainflat, y_pred_train, 1))(np.unique(y_trainflat)), "#F8766D")
plt.ylabel('Predicted LogS', fontsize=12)
plt.xlabel('Experimental LogS', fontsize=12)

# 1 row, 2 column, plot 2
ax2 = plt.subplot(1, 3, 2)
ax2.scatter(x=y_testflat, y=y_pred_test, c="#619CFF", alpha=0.3)
plt.ylabel('Predicted LogS', fontsize=12)
plt.xlabel('Experimental LogS', fontsize=12)
ax2.plot(np.unique(y_testflat), np.poly1d(np.polyfit(y_testflat, y_pred_test, 1))(np.unique(y_testflat)), "#F8766D")

# 2 row, 1 column, plot 3
ax3 = plt.subplot(1, 3, 3)
ax3.scatter(x=y_valflat, y=y_pred_val, c='#7CAE00', alpha=0.3)
ax3.plot(np.unique(y_valflat), np.poly1d(np.polyfit(y_valflat, y_pred_val, 1))(np.unique(y_valflat)), "#F8766D")

plt.ylabel('Predicted LogS', fontsize=12)
plt.xlabel('Experimental LogS', fontsize=12)
plt.axis('equal')

plt.subplots_adjust(left=0.125,
                    bottom=0.1,
                    right=0.9,
                    top=0.9,
                    wspace=0.2,
                    hspace=0.35)
#Save as custom .png title
plt.savefig('ECFPrbfExpvsPred.png')
print(plt.show())

