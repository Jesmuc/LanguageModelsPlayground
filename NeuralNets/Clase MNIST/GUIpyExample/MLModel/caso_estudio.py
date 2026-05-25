
import pandas
import numpy as np
import copy
import matplotlib.pyplot as plt
from random import random

###Reading and plotting data###
longidata=np.array(['Variables_Horno.csv'])
dl = pandas.read_csv(longidata[0])
dl.head()
#print(dlongi.loc[0])
d0 = dl.values

#print('Number of attributes (columns):', len(dl.columns))
ncol = len(dl.columns)
y = d0[:, 0]
#x = d0[:, 1:11]

x = d0[:, 1:(ncol-20)]
x2 = d0[:, 1:(ncol)]
fig, axs = plt.subplots(10)
axs[0].plot(y)
axs[1].plot(x2[:, (ncol-20)])
axs[2].plot(x2[:, (ncol-19)])
axs[3].plot(x2[:, (ncol-18)])
axs[4].plot(x2[:, (ncol-17)])
axs[5].plot(x2[:, (ncol-16)])
axs[6].plot(x2[:, (ncol-15)])
axs[7].plot(x2[:, (ncol-14)])
axs[8].plot(x2[:, (ncol-13)])
axs[9].plot(x2[:, (ncol-12)])

#plt.show()

### Prediction model ###

from sklearn import datasets
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelBinarizer
from sklearn.neural_network import MLPRegressor
from sklearn.linear_model import LinearRegression
from sklearn.metrics import accuracy_score
import numpy as np
import time
import pickle
#Generate a partition in training and test data
validation_size = 0.20
seed = 7
x_train, x_test, y_train, y_test = train_test_split(x, y,
    test_size=validation_size, random_state=seed)

clf = MLPRegressor(solver='adam', alpha=1e-5, tol=1e-50, max_iter=3000, activation='logistic',
                    hidden_layer_sizes=(100,), verbose=False, random_state=1)
#clf = LinearRegression()                   
clf.fit(x_train, y_train)

#Make predictions and report percentual prediction score
from sklearn.metrics import mean_squared_error
#from sklearn.metrics import mean_absolute_percentage_error
#print(clf.score(x_test, y_test)*100)
y_pred = clf.predict(x_test)
print('MSRE:', mean_squared_error(y_pred, y_test, squared=False)/np.average(y_test))
print('Percentage score:', (1.0 - mean_squared_error(y_pred, y_test, squared=False)/np.average(y_test))*100)
#print('MAPE:', mean_absolute_percentage_error(y_pred, y_test))

plt.show()
