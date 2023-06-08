<h1>Cómo instalar y ejecutar</h1>

Para la ejecución del programa con MPI4py es necesario tener instalado MPI en el computador. Se describirán los pasos a seguir:

    Para verificar si ya se tiene instalado MPI
      mpiexec --version
    
Si no está instalado, leer el siguiente informe: [Python MPI: Setup](https://nyu-cds.github.io/python-mpi/setup/#:~:text=Go%20to%20the%20installation%20page,executable%20and%20follow%20the%20instructions.&text=If%20you%20want%20to%20set,the%20output%20is%20as%20expected.](url))
        
<h2>Comandos:</h2>
  
      Los siguientes comandos se utilizan para la ejecución del programa
        --input1
        --input2
        --output
        --sequential
        --multiprocessing
        --processes
        --mpi
   
<h2>Ejemplos de comando: </h2>
   Para MPI4py se necesita agregar al principio el mpiexec -n [num procesos] para que MPI lo reconozca. Se ejecuta así:
   
    mpiexec -n 4 python Aplicacion.py --input1 Salmonella150.fna --input2 Ecoli150.fna --output dotplot11.png --mpi
    
  Multiprocessing se ejecuta de la siguiente manera:
  
    python Aplicacion.py --processes 4 --input1 Salmonella150.fna --input2 Ecoli150.fna --output dotplot11.png --multiprocessing
    
  Secuencial se ejecuta de la siguiente manera:
  
    python Aplicacion.py --input1 Salmonella150.fna --input2 Ecoli150.fna --output dotplot11.png --sequential
  
      
