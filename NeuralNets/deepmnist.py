import pandas
import numpy as np
import copy
import matplotlib.pyplot as plt
from random import random

'''
data = np.loadtxt('data/visdata.txt', unpack=False)
x = data[:,[0,2]]
y = data[:, [1]]
split = int(0.7*len(x))

xtr = x[0:split,:]
xte = x[split:,:]
ytr = y[0:split]
yte = y[split:]
'''

from keras.utils import to_categorical

data = np.loadtxt('mnist_train.csv', delimiter=',') 
print("Lectura de la base de datos completa")
ncol = data.shape[1]
xtr = data[:,1:ncol]
ytr = data[:,0]
ytr = to_categorical(ytr)

data = np.loadtxt('mnist_test.csv', delimiter=',') 
print("Lectura de la base de datos completa")
ncol = data.shape[1]
xte = data[:,1:ncol]
yte = data[:,0]
yte = to_categorical(yte)

from sklearn.neural_network import MLPRegressor
from sklearn.neural_network import MLPClassifier

import tensorflow as tf

#mlp = MLPRegressor(solver='adam', hidden_layer_sizes=400, max_iter=15000, tol=1e-20 , verbose=True)  
#mlp = MLPClassifier(solver='adam', hidden_layer_sizes=1, max_iter=15000, tol=1e-20 , verbose=True) 
#tmodel = mlp.fit(xtr, ytr)

model = tf.keras.models.Sequential([
  tf.keras.layers.Flatten(),
  tf.keras.layers.Dense(128, activation=tf.nn.relu),
  tf.keras.layers.Dense(128, activation=tf.nn.relu),
  tf.keras.layers.Dense(128, activation=tf.nn.relu),
  tf.keras.layers.Dropout(0.2),
  #tf.keras.layers.Dense(1, activation = 'hard_sigmoid')
  tf.keras.layers.Dense(10, activation = 'softmax')
])
model.compile(optimizer='adam',
              loss='binary_crossentropy',
              metrics=['accuracy'])
#basic_model.compile(loss = 'binary_crossentropy', optimizer = 'sgd', metrics = ['accuracy'])
model.fit(xtr, ytr, epochs=50)

# Calculate the fitness in both train and test datasets

# Evaluate the model on the test data using `evaluate`
print("Evaluate on test data")
results = model.evaluate(xte, yte, batch_size=128)
print("test loss, test acc:", results)

'''
ypred = model.predict(xtr)
#print(int(ypred))
print('Fitness in training region:', fitness(ypred, ytr))
ypred = model.predict(xte)
print('Fitness in test region:', fitness(ypred, yte))
'''
