Installation
=================

Requirements
---------------------

* 64 bit Linux OS
* CMake version 3.1 or above
* C++ compiler with C++14 support (GCC 5.0+)
* GNU make
* zlib
* python3 with the following packages installed: networkx, llist, numpy, edlib, joblib, Bio, pydot


Downloading and compiling from source code
-------------------------------------

To build from source code run the following commands from code directory.


``` bash

    cmake .
    make 
    make install
```

Binary files will be stored in the bin subdirectory.
We will further assume that you added bin directory to your PATH variable.
Otherwise, you need to use full path to the executable instead of lja in the command line.

Running La Jolla Assembler
=================

Input
-------------------------------------
Input reads can be in fasta or fastq format and can be compressed using gzip.
Make dure that file name extensions correspond to contents of files since they are used to determine file format.
E.g. uncompressed fasta reads are expected to be stored in files with name extensions ".fasta" or ".fa" while compressed fastq reads are expected to be stored in files with name extensions ".fq.gz" or "fastq.gz".

Command line
-------------------------------------
To run LJA use the following command line

``` bash

    lja [options] -o <output_dir> --reads <reads_file> [--reads <reads_file2> ...]
```

## Basic options

`-o <file_name>` (or `--output-dir <file_name>`)
    Name of output folder. Resulting graph will be stored there.

`--reads <file_name>`
    Name of file that contains reads in fasta or fastq format. This option can be used any number of times in the same command line resulting in collecting reads from multiple files.

`--help`
    Print help message.
## Advanced options
`-t <int>` (or `--threads <int>`)
    Number of threads. The default value is 16.

`-k <int>`
Value of k (vertex size) to be used for the initial error correction. k should be odd (otherwise k + 1 is used instead).

`-K <int>`
Value of k (vertex size) to be used for the final error correction and initialization of multiDBG. K should be odd (otherwise K + 1 is used instead).

`--diploid`
Use this option for diploid genome. By default, LJA assumes that the genome is haploid or inbred.

Output of de Bruijn graph construction
=================

All output files are stored in <output_dir> `, which is set by the user.

-   `<output_dir>/graph.fasta` sequences of all edges of de Bruijn graph in fasta format
-   `<output_dir>/graph.gfa` sequences of all edges of de Bruijn graph in gfa format. Note that gfa format is ill-suited for de Bruijn graph storage since it does not represent some vertices such as vertices with indegree 2 and outdegree 0. 
-   `<output_dir>/graph.dot` sequences of all edges of de Bruijn graph
-   `<output_dir>/dbg.log` log file for the run

Feedback and bug reports
=================

Your comments, bug reports, and suggestions are very welcomed.
They will help us to further improve La Jolla Assembler.
If you have any troubles  LJA, please send us `dbg.log` file from the directory `<output_dir>`.

You can send your comments and bug reports via e-mail: [anton.bankevich@gmail.com](mailto:anton.bankevich@gmail.com).

