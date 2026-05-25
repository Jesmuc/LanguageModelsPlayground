import PySimpleGUI as sg
import pandas as pd
from datetime import date
import numpy as np

# Add some color to the window
sg.theme('DarkTeal9')

EXCEL_FILE = 'prototipo1b.xls'
df = pd.read_excel(EXCEL_FILE)


d1 = df.iloc[:,0]
#d2 = df.iloc[:,2]
d1 = d1.values
#d2 = d2.values

#x = np.vstack((d1,d2))
#x = np.asarray(x)
x = d1
#np.transpose(x)
x = [x]
x = np.array(x).reshape(1,-1)
#print(len(x))
#print(x)


layout = [
    [sg.Text('Antes de proceder, por favor verificar las entradas del archivo prototipo1b.xls')],
    #[sg.Text('Archivo', size=(15,1)), sg.InputText(key='data')],
    
    [sg.Submit('Proceder'), sg.Exit('Salir')],
    
    [sg.Text('Este es un prototipo de muestra. Quedará inactivo el 01/7/2026')],
]

window = sg.Window('Prototipo 1.0', layout)

def clear_input():
    for key in values:
        window[key]('')
    return None


while True:
    delta = date(2022, 7, 1) - date.today()
    event, values = window.read()
    if event == sg.WIN_CLOSED or event == 'Salir' or date.today() > date(2026, 7, 1):
        break
    #if event == 'Borrar':
     #   clear_input()
    if event == 'Proceder':
        import pickle
        with open('model.pickle', 'rb') as f:
             model = pickle.load(f)
             pred = model.predict(x)
        #df = df.append(values, ignore_index=True)
        df.iat[13,4] = pred
        df.to_excel(EXCEL_FILE, index=False)
        sg.popup('Se pronostica un tiempo de PowerOn de:', pred, 'minutos. El archivo prototipo1.xls ha sido actualizado con              esta información. Le recordamos que este prototipo tiene', delta.days , 'días de vigencia.')
        clear_input()
window.close()
