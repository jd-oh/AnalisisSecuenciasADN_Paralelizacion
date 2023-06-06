import numpy as np
import time
from tqdm import tqdm # esta librería es para mirar el progreso de un for

cant_muestras = 150

# Definimos unas secuencias
Secuencia1 = "ACGTCGTCGAGCTAGCATCGATCAGNNNCATCATCNACTATACNNNNCATCATCATCTACTGCTACGACTACGAGAGAGCTACGACTACG"*cant_muestras
Secuencia2 = "NGCNATCACGATGCATGCACTACGATCGACAGCATCGATCGATGCATCATGCATCGNATGCNTGASCSATCGACGTANGCACTGACNTGA"*cant_muestras
Secuencia2 = Secuencia1

#SecuenciaLTR = "CACTAGACTAGACTAGCNAGCTACGCATGGCTACNCTACGACAGCTAGCTANCTATCNACTACNAGCTACTAGCTANNNACTANCTCGACTACGACTACACTGACCACTAGAC"*cant_muestras
#Secuencia1, Secuencia2 = SecuenciaLTR, SecuenciaLTR

begin = time.time()
dotplot = np.empty([len(Secuencia1),len(Secuencia2)])
print("La matriz de resultado tiene tamaño: ", dotplot.shape)

for i in tqdm(range(dotplot.shape[0])):
  for j in range(dotplot.shape[1]):
    if Secuencia1[i] == Secuencia2[j]:
      dotplot[i,j] = 1
    else:
      dotplot[i,j] = 0

print(f"\n El código se ejecutó en: {time.time() - begin} segundos")

import matplotlib.pyplot as plt
def draw_dotplot(matrix, fig_name='dotplot.svg'):
  plt.figure(figsize=(5,5))
  plt.imshow(matrix, cmap='Greys',aspect='auto')

  plt.ylabel("Secuencia 1")
  plt.xlabel("Secuencia 2")
  plt.savefig(fig_name)

  draw_dotplot(dotplot[:500,:500 ])