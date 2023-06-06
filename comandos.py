import argparse
from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing as mp
from mpi4py import MPI
import time
from tqdm import tqdm

# Función para calcular el dotplot secuencialmente
def calculate_dotplot_sequential(Secuencia1, Secuencia2):
    # Cálculo del dotplot secuencial
    begin = time.time()
    dotplot = np.empty([len(Secuencia1),len(Secuencia2)], dtype=np.int8)
    print("La matriz de resultado tiene tamaño: ", dotplot.shape)

    for i in tqdm(range(dotplot.shape[0])):
        for j in range(dotplot.shape[1]):
            if Secuencia1[i] == Secuencia2[j]:
                dotplot[i,j] = 1
            else:
                dotplot[i,j] = 0

    print(f"\n El código se ejecutó en: {time.time() - begin} segundos")
    return dotplot

# Función para calcular el dotplot utilizando multiprocessing
def calculate_dotplot_parallel(Secuencia1, Secuencia2):
    # Cálculo del dotplot utilizando multiprocessing
    begin = time.time()
    dotplot = np.array(parallel_dotplot(Secuencia1, Secuencia2,4))

    print("La matriz de resultado tiene tamaño: ", dotplot.shape)
    print(f"\n El código se ejecutó en: {time.time() - begin} segundos")
    return dotplot

def worker(args):
    i, Secuencia1, Secuencia2 = args
    return [Secuencia1[i] == Secuencia2[j] for j in range(len(Secuencia2))]

def parallel_dotplot(Secuencia1, Secuencia2, threads=mp.cpu_count()):
    with mp.Pool(processes=threads) as pool:
        result = pool.map(worker, [(i, Secuencia1, Secuencia2) for i in range(len(Secuencia1))])
    return result

# Función para calcular el dotplot utilizando mpi4py
def calculate_dotplot_mpi(Secuencia1, Secuencia2):
    begin = time.time()

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    # Dividir la secuencia1 en chunks, uno por cada proceso.
    chunks = np.array_split(range(len(Secuencia1)), size)

    # Cálculo del dotplot utilizando mpi4py
    dotplot = np.empty([len(chunks[rank]),len(Secuencia2)],dtype=np.int32)


    for i in range(len(chunks[rank])):
        for j in range(len(Secuencia2)):
            if Secuencia1[chunks[rank][i]] == Secuencia2[j]:
                dotplot[i,j] = np.int32(1)
            else:
                dotplot[i,j] = np.int32(0)

    # Recopilar los resultados de todos los procesos
    dotplot = comm.gather(dotplot, root=0)

    if rank == 0:
        # Combinar los resultados en un solo dotplot
        merged_data = np.vstack(dotplot)
        end = time.time()
        print(f"Tiempo total de ejecución: {end-begin} segundos")

        return merged_data

# Función para guardar el dotplot en un archivo de imagen
def save_dotplot(dotplot, output_file):
    # Guardar el dotplot en un archivo de imagen
    plt.figure(figsize=(5,5))
    plt.imshow(dotplot[:500,:500], cmap='Greys',aspect='auto')

    plt.ylabel("Secuencia 1")
    plt.xlabel("Secuencia 2")
    plt.savefig(output_file)

    #save_dotplot(dotplot[:500,:500 ])

def merge_sequences_from_fasta(file_path):
    sequences = []  # List to store all sequences
    for record in SeqIO.parse(file_path, "fasta"):
        # `record.seq` gives the sequence
        sequences.append(str(record.seq))
    return "".join(sequences)


if __name__ == "__main__":
    # Configurar la línea de comandos
    parser = argparse.ArgumentParser(description='Dotplot analysis')
    parser.add_argument('sequence1', help='Path to the first sequence in FASTA format')
    parser.add_argument('sequence2', help='Path to the second sequence in FASTA format')
    parser.add_argument('--sequential', action='store_true', help='Calculate dotplot sequentially')
    parser.add_argument('--multiprocessing', action='store_true', help='Calculate dotplot using multiprocessing')
    parser.add_argument('--mpi', action='store_true', help='Calculate dotplot using mpi4py')
    parser.add_argument('--output', default='dotplot.png', help='Path to the output dotplot image file')
    args = parser.parse_args()

    # Leer las secuencias desde los archivos FASTA
    
    seq1 = str(SeqIO.read(args.sequence1, 'fasta').seq)
    seq2 = str(SeqIO.read(args.sequence2, 'fasta').seq)
    
    """
    file_path_1 = "Ecoli2.fna"
    file_path_2 = "Salmonella1.fna"

    seq1 = merge_sequences_from_fasta(file_path_1)
    seq2 = merge_sequences_from_fasta(file_path_2)
    """
    

    # Calcular el dotplot según la opción seleccionada
    if args.sequential:
        dotplot = calculate_dotplot_sequential(seq1, seq2)
    elif args.multiprocessing:
        dotplot = calculate_dotplot_parallel(seq1, seq2)
    elif args.mpi:
        dotplot = calculate_dotplot_mpi(seq1, seq2)
    else:
        print("Debes seleccionar al menos una opción para calcular el dotplot.")
        exit()

    # Guardar el dotplot en un archivo de imagen
    save_dotplot(dotplot, args.output)
