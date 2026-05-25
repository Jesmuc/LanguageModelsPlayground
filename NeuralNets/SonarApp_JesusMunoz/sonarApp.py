import os
import pickle
import joblib
import tkinter as tk
from tkinter import messagebox, ttk
from PIL import Image, ImageTk
import pandas as pd
import numpy as _np

# -------- Configurables --------
DATA_URL = "https://goo.gl/NXoJfR"
MODEL_FILE = 'modelRNA.sav'
IMG_ROCK = 'rock.png'
IMG_MINE = 'mine.png'
# --------------------------------

def load_dataset_from_url_or_local(url):
    df = pd.read_csv(url, header=None)
    # tomar primeras 60 columnas como features y la columna 61 (index 60) como label
    features = df.iloc[:, :60].astype(float).values
    labels = df.iloc[:, 60].values
    return features, labels

def normalize_label_raw(l):
    s = str(l).strip()
    return s.upper()

class SonarApp:
    def __init__(self, master):
        self.master = master
        master.title("Sonar: Mines vs Rock dataset")

        # cargar modelo
        self.model = pickle.load(open('modelRNA.sav', 'rb'))

        # cargar dataset desde URL (o fallback local)
        self.features, self.labels_raw = load_dataset_from_url_or_local(DATA_URL)

        # normalizar labels
        self.labels = [normalize_label_raw(l) for l in self.labels_raw]
        self.n_instances = len(self.features)

        # cargar imágenes
        self.img_rock = Image.open(IMG_ROCK)
        self.img_mine = Image.open(IMG_MINE)

        # --- UI widgets (se construyen tras tener n_instances) ---
        frm = ttk.Frame(master, padding=12)
        frm.grid(row=0, column=0, sticky="nsew")

        ttk.Label(frm, text=f"Número de instancia (1..{self.n_instances}):").grid(row=0, column=0, sticky="w")
        self.spin = tk.Spinbox(frm, from_=1, to=self.n_instances, width=6)
        self.spin.grid(row=0, column=1, sticky="w")

        self.btn_predict = ttk.Button(frm, text="Predecir", command=self.on_predict)
        self.btn_predict.grid(row=0, column=2, padx=8)

        self.btn_reset = ttk.Button(frm, text="Reintentar", command=self.on_reset)
        self.btn_reset.grid(row=0, column=3, padx=8)

        self.msg_var = tk.StringVar(value="Selecciona un número y pulsa Predecir")
        ttk.Label(frm, textvariable=self.msg_var).grid(row=1, column=0, columnspan=4, sticky="w", pady=(8,0))

        # panel de imágenes
        self.panel = ttk.Frame(frm, padding=(0,8))
        self.panel.grid(row=2, column=0, columnspan=4)

        ttk.Label(self.panel, text="Predicho:").grid(row=0, column=0)
        ttk.Label(self.panel, text="Real:").grid(row=0, column=1)

        self.canvas_pred = ttk.Label(self.panel)
        self.canvas_pred.grid(row=1, column=0, padx=10, pady=6)
        self.canvas_true = ttk.Label(self.panel)
        self.canvas_true.grid(row=1, column=1, padx=10, pady=6)

        self.pred_label_var = tk.StringVar(value="")
        self.true_label_var = tk.StringVar(value="")
        ttk.Label(self.panel, textvariable=self.pred_label_var).grid(row=2, column=0)
        ttk.Label(self.panel, textvariable=self.true_label_var).grid(row=2, column=1)

        self.tkimg_pred = None
        self.tkimg_true = None

        # mensaje si hay mismatch de tamaño
        if self.n_instances < 208:
            self.msg_var.set(f"Dataset cargado con {self.n_instances} instancias (esperado 208). Asegúrate que la URL/archivo es correcto.")
        # fin init

    def on_predict(self):
        try:
            idx1 = int(self.spin.get()) - 1
        except Exception:
            messagebox.showwarning("Input", "Número de instancia inválido.")
            return

        if not (0 <= idx1 < len(self.features)):
            messagebox.showwarning("Input", f"El índice debe estar entre 1 y {len(self.features)}.")
            return

        x = self.features[idx1:idx1+1]
        pred0 = self.model.predict(x)

        # Determinar etiqueta predicha:
        prob = float(pred0)
        pred_label = 'R' if prob >= 0.5 else 'M'   # 1 -> rock, 0 -> mine

        # ahora determinar imagen predicha
        pred_text = f"Pred: {pred0}"
        if pred_label == 'M':
            img_pred = self.img_mine
        elif pred_label == 'R':
            img_pred = self.img_rock

        # ground truth
        true_raw = self.labels[idx1]
        true_norm = self._normalize_true(true_raw)
        if true_norm == 'M':
            img_true = self.img_mine; true_text = f"Mina: {true_raw}"
        elif true_norm == 'R':
            img_true = self.img_rock; true_text = f"Roca: {true_raw}"

        def make_tk(img, maxsize=(240,240)):
            im = img.copy()
            im.thumbnail(maxsize, Image.ANTIALIAS)
            return ImageTk.PhotoImage(im)

        self.tkimg_pred = make_tk(img_pred)
        self.tkimg_true = make_tk(img_true)

        self.canvas_pred.config(image=self.tkimg_pred)
        self.canvas_true.config(image=self.tkimg_true)

        self.pred_label_var.set(pred_text)
        self.true_label_var.set(true_text)
        self.msg_var.set(f"Instancia {idx1+1}: predicho vs real mostrados.")

    def _normalize_true(self, v):
        return v.strip().upper()

    def on_reset(self):
        self.canvas_pred.config(image='')
        self.canvas_true.config(image='')
        self.tkimg_pred = None
        self.tkimg_true = None
        self.pred_label_var.set('')
        self.true_label_var.set('')
        self.msg_var.set("Selecciona un número y pulsa Predecir")

if __name__ == '__main__':
    root = tk.Tk()
    app = SonarApp(root)
    root.mainloop()