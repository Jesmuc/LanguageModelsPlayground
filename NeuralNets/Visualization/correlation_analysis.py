#Este programa carga despliega un analisis de correlaciones entre atributos de una base de datos
import pandas
import numpy as np
import matplotlib.pyplot as plt
names = [ ' criminalidad ' , ' lotes por zona ' , ' porcentaje de negocios ' , ' rio ' , 
' conc. oxidos ' , ' # habitaciones ' , ' edad de propiedad ' , ' distancia a ctros. empleo ' , 
' acceso a vias ', ' predial ', ' razon maestros por alumno', ' tipo de poblacion ', ' porcentaje pob. bajo ingreso ' , 'valor de propiedad en miles de usd'  ]
names2 = ['1' , '2' , '3' , '4' , '5' , '6' , '7' , '8' , '9' , '10' , '11' , '12' , '13', '14']
data = pandas.read_csv('boston.csv', names=names)
pandas.set_option( 'precision' , 3)
correlations = data.corr(method= 'pearson')
print(correlations)
np.savetxt('correlations.csv', correlations, delimiter=',') 
# plot correlation matrix
fig = plt.figure()
ax = fig.add_subplot(111)
cax = ax.matshow(correlations, vmin=-1, vmax=1)
fig.colorbar(cax)
ticks = np.arange(0,14,1)
ax.set_xticks(ticks)
ax.set_yticks(ticks)
ax.set_xticklabels(names2)
ax.set_yticklabels(names)
plt.suptitle('Correlaciones', fontsize=14)
plt.show()
