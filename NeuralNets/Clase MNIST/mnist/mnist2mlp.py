# Este programa entrena una red neuronal de 3 capas para MNIST y guarda el modelo entrenado para posterior uso.

from sklearn import datasets
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelBinarizer
from sklearn.neural_network import MLPClassifier
from sklearn.metrics import accuracy_score
from sklearn.linear_model import LogisticRegression
import numpy as np
import time
import pickle


# leyendo la base de datos MNIST
data = np.loadtxt('mnist_train.csv', delimiter=',') 
print("Lectura de la base de datos completa")
ncol = data.shape[1]
# definiendo entradas y salidas
X = data[:,1:ncol]
y = data[:,0]

#tic = time.clock()
clf = MLPClassifier(solver='adam', alpha=1e-5, tol=1e-50, max_iter=300, activation='logistic',
                    hidden_layer_sizes=(100,), verbose=True, random_state=1)
clf.fit(X, y)
#toc = time.clock()
print("Entrenamiento completo")
#print("Tiempo de procesador para el entrenamiento (seg):")
#print(toc - tic)
        

data = np.loadtxt('mnist_test.csv', delimiter=',') 
#print(data)
ncol = data.shape[1]
# definiendo entradas y salidas
X_test = data[:,1:ncol]
y_test = data[:,0]

# predecir los valores de X_test
predicted = clf.predict(X_test)
	
# para finalizar se calcula el error
error = 1 - accuracy_score(y_test, predicted)
print("Error en el conjunto de prueba:")
print(error)

# el modelo entrenado se salva en disco
filename = 'finalized_modelMLP.sav'
pickle.dump(clf, open(filename, 'wb'))
