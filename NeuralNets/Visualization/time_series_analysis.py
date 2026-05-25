#Este programa carga despliega visualmente un conjunto de datos en formato csv
import pandas
from pandas import Series
import numpy as np
import matplotlib.pyplot as plt

filename=np.array(['renewable_energy.csv', 'petroleum.csv', 'electric_cars.csv'])

i = 1
for f in filename:
    series = pandas.read_csv(f)
    series.plot()

print(series)
plt.show()

