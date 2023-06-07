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
import os
import gc

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

    print(f"\n El código se ejecutó en: {time.time() - begin} segundos")
    return dotplot


def dotplot_chunk_and_filter(args):
    chunk, Secuencia1, Secuencia2 = args
    dotplot = np.empty([len(chunk), len(Secuencia2)], dtype=np.int8)
    for i in range(len(chunk)):
        for j in range(len(Secuencia2)):
            if Secuencia1[chunk[i]] == Secuencia2[j]:
                dotplot[i, j] = np.int8(1)
            else:
                dotplot[i, j] = np.int8(0)

    filtered_dotplot = filter_dotplot(dotplot)
    return filtered_dotplot


def dotplot_generator(num_processes, chunks, Secuencia1, Secuencia2):
    with mp.Pool(num_processes) as p:
        for filtered_dotplot in p.imap(dotplot_chunk_and_filter, [(chunk, Secuencia1, Secuencia2) for chunk in chunks]):
            yield filtered_dotplot


def calculate_dotplot_parallel(Secuencia1, Secuencia2, output_file):
    begin = time.time()
    num_processes = 4
    chunks = np.array_split(range(len(Secuencia1)), num_processes)

    # Create a temporary file to store the partial dotplot
    fp = np.memmap(output_file, dtype='int8', mode='w+',
                   shape=(len(Secuencia1), len(Secuencia2)))

    row_start = 0
    for chunk_dotplot in dotplot_generator(num_processes, chunks, Secuencia1, Secuencia2):
        row_end = row_start + chunk_dotplot.shape[0]
        fp[row_start:row_end, :] = chunk_dotplot
        row_start = row_end

    print(f"\n El código se ejecutó en: {time.time() - begin} segundos")
    return fp


# Función para calcular el dotplot utilizando mpi4py


def filter_dotplot(dotplot, window_size=5, threshold=0.7):
    rows, cols = dotplot.shape
    binary_dotplot = np.zeros_like(dotplot, dtype=np.int8)

    diagonal_sums = np.zeros_like(dotplot, dtype=np.float16)
    for idx in range(window_size):
        diagonal_sums += np.roll(np.roll(dotplot,
                                 shift=idx, axis=0), shift=idx, axis=1)
    diagonal_sums /= window_size

    binary_dotplot = np.where(diagonal_sums > threshold, 1, 0)

    return binary_dotplot
# Función para calcular el dotplot utilizando mpi4py


def calculate_dotplot_mpi(Secuencia1, Secuencia2):
    begin = time.time()

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    # Dividir la secuencia1 en chunks, uno por cada proceso.
    chunks = np.array_split(range(len(Secuencia1)), size)

    dotplot = np.empty([len(chunks[rank]), len(Secuencia2)], dtype=np.int8)

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

        # Apply the filter_dotplot function to the merged_data
        filtered_merged_data = filter_dotplot(merged_data)

        end = time.time()
        print(f"Tiempo total de ejecución: {end-begin} segundos")

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
    seq1 = str(SeqIO.read(args.input1, "fasta").seq)
    seq2 = str(SeqIO.read(args.input2, "fasta").seq)

# Calcular el dotplot según la opción seleccionada
    output_file = os.path.splitext(args.output)[0] + "_temp.npy"
    if args.sequential:
        dotplot = calculate_dotplot_sequential(seq1, seq2)
    elif args.multiprocessing:
        # Add the missing output_file argument
        dotplot = calculate_dotplot_parallel(seq1, seq2, output_file)
    elif args.mpi:
        dotplot = calculate_dotplot_mpi(seq1, seq2)
    else:
        print("Debes seleccionar al menos una opción para calcular el dotplot.")
        exit()

    # Guardar el dotplot en un archivo de imagen
    plt.figure(figsize=(10, 10))
    plt.imshow(dotplot[:500, :500], cmap='Greys', aspect='auto')
    plt.savefig(args.output)

    # Close and release memory for dotplot
    del dotplot
    gc.collect()

    # Remove the temporary file
    if os.path.exists(output_file):
        os.remove(output_file)
