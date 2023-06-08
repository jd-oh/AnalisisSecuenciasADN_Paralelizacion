<h2>Instalación y requerimientos</h2>

• Para la ejecución del programa con MPI4py es necesario tener instalado MPI en el computador. Se describirán los pasos a seguir:

    Para verificar si ya se tiene instalado MPI
      mpiexec --version
    
• Si no está instalado, leer el siguiente informe: [Python MPI: Setup](https://nyu-cds.github.io/python-mpi/setup/#:~:text=Go%20to%20the%20installation%20page,executable%20and%20follow%20the%20instructions.&text=If%20you%20want%20to%20set,the%20output%20is%20as%20expected.](url))

<b>• Los archivos de secuencia deben estar en formato FASTA (.fna)</b>
        
<h2>Comandos:</h2>
  
      Los siguientes comandos se utilizan para la ejecución del programa
        --input1 rutaSecuencia1
        --input2 rutaSecuencia2
        --output nombreDeLaImagenResultante
        --sequential
        --multiprocessing
        --processes num_procesadores
        --mpi
   
<h2>Ejemplos de comando: </h2>
 <h3>Si los archivos FASTA están en la carpeta raíz, se ejecutaría así dependiendo de la forma en que se va a calcular:</h3>
 
  • Para MPI4py se necesita agregar al principio el mpiexec -n [num procesos] para que MPI lo reconozca. Se ejecuta así:
   
    mpiexec -n 4 python Aplicacion.py --input1 Salmonella150.fna --input2 Ecoli150.fna --output dotplot11.png --mpi
    
  • Para Multiprocessing se necesita especificar el número de procesos. Se ejecuta de la siguiente manera:
  
    python Aplicacion.py --processes 4 --input1 Salmonella150.fna --input2 Ecoli150.fna --output dotplot11.png --multiprocessing
    
  • Secuencial se ejecuta de la siguiente manera:
  
    python Aplicacion.py --input1 Salmonella150.fna --input2 Ecoli150.fna --output dotplot11.png --sequential
  
  
<h2>Salida</h2>
En la salida veremos una imagen filtrada con las coincidencias que se encontraron entre las dos secucuencias (lineas diagonales). Esta se guardará en la carpeta raíz del proyecto


