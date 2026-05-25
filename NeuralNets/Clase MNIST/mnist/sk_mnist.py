import numpy as np
import sklearn
from sklearn.neural_network import MLPClassifier

# reading MNIST train database
data = np.loadtxt('mnist_train.csv', delimiter=',') 
print("Lecture of the database completed")
ncol = data.shape[1]
# defining inputs and outputs
x = data[:,1:ncol]
y = data[:,0]

model = sklearn.neural_network.MLPClassifier(
    activation='relu', solver='adam', max_iter=500, hidden_layer_sizes=(10), verbose=True, n_iter_no_change=100)
model.fit(x, y.ravel())
ex = 0
y_pred = model.predict(x[ex,:].reshape(1, -1))
print('Prediction:',y_pred)
print('Target:',y[ex])

