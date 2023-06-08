import argparse
from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing as mp
from mpi4py import MPI
import time
from tqdm import tqdm
from multiprocessing import Pool
from scipy.signal import convolve2d

# Función para calcular el dotplot secuencialmente


def calculate_dotplot_sequential(Secuencia1, Secuencia2):
    # Cálculo del dotplot secuencial
    begin = time.time()
    dotplot = np.empty([len(Secuencia1), len(Secuencia2)], dtype=np.int8)
    print("La matriz de resultado tiene tamaño: ", dotplot.shape)

    for i in tqdm(range(dotplot.shape[0])):
        for j in range(dotplot.shape[1]):
            if Secuencia1[i] == Secuencia2[j]:
                dotplot[i, j] = 1
            else:
                dotplot[i, j] = 0

    filtered_dotplot = filter_dotplot(dotplot)

    print(f"\n El código en secuencial se ejecutó en: {time.time() - begin} segundos")
    return filtered_dotplot


def calculate_dotplot_parallel(Secuencia1, Secuencia2, threads=4):
    chunks = np.array_split(range(len(Secuencia1)), threads)
    begin = time.time()

    with Pool(processes = threads) as p:
        dotplot_list = p.starmap(
            dotplot_chunk, [(chunk, Secuencia1, Secuencia2) for chunk in chunks])

    dotplot = np.vstack(dotplot_list)
    end = time.time()



    # Aplicar el filtro al dotplot
    filtered_dotplot = filter_dotplot(dotplot)


    print("La matriz de resultado tiene tamaño: ", filtered_dotplot.shape)
    print(f"\n El código en multiprocessing(porción paralela) se ejecutó en: {end - begin} segundos")
    return filtered_dotplot


def dotplot_chunk(chunk, Secuencia1, Secuencia2):
    dotplot = np.empty([len(chunk), len(Secuencia2)], dtype=np.int8)
    for i in range(len(chunk)):
        for j in range(len(Secuencia2)):
            if Secuencia1[chunk[i]] == Secuencia2[j]:
                dotplot[i, j] = np.int8(1)
            else:
                dotplot[i, j] = np.int8(0)
    return dotplot

# Función para calcular el dotplot utilizando mpi4py


def filter_dotplot(dotplot):
    window_size = 5
    window = np.eye(window_size, dtype=np.float32)
    print("Empecé a haceer el convolve 2d")
    print("Empecé a haceer el convolve 2d")
    print("Empecé a haceer el convolve 2d")
    filtered_dotplot = convolve2d(dotplot, window, mode='same')
    print("Terminé de hacer el convolve 2d")
    print("Terminé de hacer el convolve 2d")
    print("Terminé de hacer el convolve 2d")
    filtered_dotplot = filtered_dotplot.astype(np.float16)
    threshold = np.float16(0.8 * window_size)
    binary_dotplot = (filtered_dotplot >= threshold).astype(np.int8)
    return binary_dotplot
# Función para calcular el dotplot utilizando mpi4py


def calculate_dotplot_mpi(Secuencia1, Secuencia2):
    #beginTotal = time.time()
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    # Dividir la secuencia1 en chunks, uno por cada proceso.
    chunks = np.array_split(range(len(Secuencia1)), size)

    dotplot = np.empty([len(chunks[rank]), len(Secuencia2)], dtype=np.int8)

    begin = time.time()

    for i in range(len(chunks[rank])):
        for j in range(len(Secuencia2)):
            if Secuencia1[chunks[rank][i]] == Secuencia2[j]:
                dotplot[i, j] = np.int8(1)
            else:
                dotplot[i, j] = np.int8(0)

    # gather data from all processes onto the root process
    dotplot = comm.gather(dotplot, root=0)

    # The root process prints the results and generates the plot.
    if rank == 0:
        # merge the gathered data into a single array
        merged_data = np.vstack(dotplot)
        end = time.time()
        # Apply the filter_dotplot function to the merged_data
        filtered_merged_data = filter_dotplot(merged_data)

        #endTotal = time.time()

        #print(f"Tiempo de ejecución total en MPI: {endTotal-beginTotal} segundos")

        
        print(f"Tiempo de ejecución en MPI(porción paralela): {end-begin} segundos")

        return filtered_merged_data

# Función para guardar el dotplot en un archivo de imagen


def save_dotplot(dotplot, output_file):
    # Guardar el dotplot en un archivo de imagen
    plt.figure(figsize=(5, 5))
    plt.imshow(dotplot[:500, :500], cmap='Greys', aspect='auto')

    plt.ylabel("Secuencia 1")
    plt.xlabel("Secuencia 2")
    plt.savefig(output_file)

    # save_dotplot(dotplot[:500,:500 ])


def merge_sequences_from_fasta(file_path):
    sequences = []  # List to store all sequences
    for record in SeqIO.parse(file_path, "fasta"):
        # `record.seq` gives the sequence
        sequences.append(str(record.seq))
    return "".join(sequences)

"""
def calculo_tiempos_multiprocessing(Secuencia1, Secuencia2):
    n_proc = [1,2, 3, 4]
    times = []
    for i in n_proc:
        begin_paralelo = time.time()
        calculate_dotplot_parallel(Secuencia1, Secuencia2,i)
        end_paralelo = time.time()
        times.append(end_paralelo-begin_paralelo)   
        print("Dotplot con ",i," procesadores, tiempo: ",end_paralelo-begin_paralelo," segundos")

    acel = [times[0]/i for i in times]
    efic = [acel[i]/n_proc[i] for i in range(len(n_proc))]
    print("Aceleración: ",acel)
    print("Eficiencia: ",efic)

    
    plt.figure(figsize=(5,5))
    plt.plot(n_proc,times)  
    plt.xlabel("Número de procesadores")
    plt.ylabel("Tiempo de ejecución")

    
    plt.figure(figsize=(10,5))
    plt.subplot(1,2,1)
    plt.plot(n_proc,times)
    plt.xlabel("Número de procesadores")
    plt.ylabel("Tiempo de ejecución")
    plt.subplot(1,2,2)
    plt.plot(n_proc,acel)
    plt.plot(n_proc,efic)
    plt.xlabel("Número de procesadores")
    plt.ylabel("Aceleración y eficiencia")
    plt.legend(["Aceleración","Eficiencia"])

    #plt.show()
    plt.savefig("tiempos2.png")
"""

def calculo_tiempos_mpi4py():
    n_proc = [1,2,3,4]
    times = [352.8144178390503, 225.89974403381348, 187.84728908538818, 171.14893579483032]
    """
    for i in n_proc:
        begin_paralelo_mpi = time.time()
        calculate_dotplot_mpi(Secuencia1, Secuencia2,i)
        end_paralelo_mpi = time.time()
        times.append(end_paralelo_mpi-begin_paralelo_mpi)   
        print("Dotplot con ",i," procesadores, tiempo: ",end_paralelo_mpi-begin_paralelo_mpi," segundos")
    """
    

    acel = [times[0]/i for i in times]
    efic = [acel[i]/n_proc[i] for i in range(len(n_proc))]
    print("Aceleración: ",acel)
    print("Eficiencia: ",efic)

    plt.figure(figsize=(5,5))
    plt.plot(n_proc,times)  
    plt.xlabel("Número de procesadores")
    plt.ylabel("Tiempo de ejecución")

    plt.figure(figsize=(10,5))
    plt.subplot(1,2,1)
    plt.plot(n_proc,times)
    plt.xlabel("Número de procesadores")
    plt.ylabel("Tiempo de ejecución")
    plt.subplot(1,2,2)
    plt.plot(n_proc,acel)
    plt.plot(n_proc,efic)
    plt.xlabel("Número de procesadores")
    plt.ylabel("Aceleración y eficiencia")
    plt.legend(["Aceleración","Eficiencia"])

    #plt.show()
    plt.savefig("tiemposMPI.png")
    
"""

def calculo_escalamiento_multiprocessing(Secuencia2):
    n_proc = [1, 2, 3, 4]
    strong_times = []  
    weak_times = []

    for i in n_proc:
        begin_paralelo = time.time()
        calculate_dotplot_parallel(Secuencia2, Secuencia2,i)
        end_paralelo = time.time()
        strong_times.append(end_paralelo - begin_paralelo)   
        print("Dotplot con ", i, " procesadores, tiempo: ", end_paralelo - begin_paralelo)

    for i in n_proc:
        Secuencia = Secuencia2[:len(Secuencia2)*i]  # Incrementa el tamaño de la secuencia.
        begin_paralelo = time.time()
        calculate_dotplot_parallel(Secuencia, Secuencia, i)
        end_paralelo = time.time()
        weak_times.append(end_paralelo - begin_paralelo)   
        print("Dotplot con ", i, " procesadores, tiempo: ", end_paralelo - begin_paralelo)

    plt.figure(figsize=(10, 5))

    plt.subplot(1, 2, 1)
    plt.plot(n_proc, strong_times, marker='o')
    plt.title("Strong Scaling")
    plt.xlabel("Number of processors")
    plt.ylabel("Time (s)")

    plt.subplot(1, 2, 2)
    plt.plot(n_proc, weak_times, marker='o')
    plt.title("Weak Scaling")
    plt.xlabel("Number of processors")
    plt.ylabel("Time (s)")

    plt.tight_layout()
    plt.savefig("EscalamientoMultiprocessing.png")

"""

def calculate_strong_scalability(Secuencia1, Secuencia2):
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    strong_times = []

    for i in range(rank + 1):
        begin_parallel = time.time()
        calculate_dotplot_mpi(Secuencia1, Secuencia2)
        end_parallel = time.time()
        strong_times.append(end_parallel - begin_parallel)
        print("Dotplot con ", i+1, " procesadores, tiempo: ", end_parallel - begin_parallel)

    return strong_times


if __name__ == "__main__":
    # Configurar la línea de comandos
    parser = argparse.ArgumentParser(description='Dotplot analysis')
    parser.add_argument('--input1', required=True, help='Input fasta file 1')
    parser.add_argument('--input2', required=True, help='Input fasta file 2')
    parser.add_argument('--output', required=True, help='Output image file')
    parser.add_argument('--sequential', action='store_true',
                        help='Run sequential dotplot')
    parser.add_argument('--multiprocessing', action='store_true',
                        help='Run dotplot using multiprocessing')
    parser.add_argument('--mpi', action='store_true',
                        help='Run dotplot using mpi4py')

    args = parser.parse_args()

    # Leer las secuencias desde los archivos fasta
    beginDatos = time.time()
    seq1 = str(SeqIO.read(args.input1, "fasta").seq)
    seq2 = str(SeqIO.read(args.input2, "fasta").seq)
    endDatos = time.time()

    print("El tiempo de carga de datos es", endDatos - beginDatos, "segundos")


    # Calcular el dotplot según la opción seleccionada
    if args.sequential:
        dotplot = calculate_dotplot_sequential(seq1, seq2)
    elif args.multiprocessing:
        dotplot = calculate_dotplot_parallel(seq1, seq2)
        #calculo_tiempos_multiprocessing(seq1, seq2)
        #calculo_escalamiento_multiprocessing(seq2)
    elif args.mpi:
        dotplot = calculate_dotplot_mpi(seq1, seq2)
        #calculo_tiempos_mpi4py()
    else:
        print("Debes seleccionar al menos una opción para calcular el dotplot.")
        exit()

    
     # Guardar el dotplot en un archivo de imagen
    beginImagen = time.time()
    plt.figure(figsize=(10, 10))
    plt.imshow(dotplot[:500, :500], cmap='Greys', aspect='auto')
    plt.savefig(args.output)
    endImagen = time.time()
    
    print("El tiempo de carga de imagenes es", endImagen - beginImagen, "segundos")
    
   


