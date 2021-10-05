Installation
=================

For installation guide please refer to the installation section of the [LJA manual](lja_manual.md).

Running jumboDBG
=================

Input
-------------------------------------
Input reads can be in fasta or fastq format and can be compressed using gzip.
Make sure that file name extensions correspond to contents of files since they are used to determine file format.
E.g. uncompressed fasta reads are expected to be stored in files with name extensions ".fasta" or ".fa" while compressed fastq reads are expected to be stored in files with name extensions ".fq.gz" or "fastq.gz".

Command line
-------------------------------------
Here we assume that you added bin subdirectory from installation directory to your PATH variable.
Otherwise, you need to use full path to the executable instead of jumboDBG in the command line.

``` bash

    jumboDBG [options] -o <output_dir> --reads <reads_file> [--reads <reads_file2> ...] -k <int>
```

## Basic options

`-o <file_name>` (or `--output-dir <file_name>`)
Name of output folder. Resulting graph will be stored there.

`-k <int>`
Value of k (vertex size) to be used for de Bruijn graph construction. k should be odd (otherwise k + 1 is used instead).

`--reads <file_name>`
Name of file that contains reads in fasta or fastq format. This option can be used any number of times in the same command line resulting in collecting reads from multiple files.

`--help`
Print help message.

## Advanced options
`-t <int>` (or `--threads <int>`)
Number of threads. The default value is 16.

`-w <int>` (or `--window <int>`)
The window size to be used for sparse de Bruijn graph construction. The default value is 2000. Note that all reads of length less than k + w are ignored during graph construction.

`--compress` Compress all homopolimers in reads.

`--coverage`
Calculate edge coverage of edges in the constructed de Bruijn graph.

Output of de Bruijn graph construction
=================

All output files are stored in <output_dir> `, which is set by the user.

-   `<output_dir>/graph.fasta` sequences of all edges of de Bruijn graph in fasta format
-   `<output_dir>/graph.gfa` sequences of all edges of de Bruijn graph in gfa format. Note that gfa format is ill-suited for de Bruijn graph storage since it does not represent some vertices such as vertices with indegree 2 and out degree 0.
-   `<output_dir>/graph.dot` sequences of all edges of de Bruijn graph
-   `<output_dir>/dbg.log` log file for the run

Feedback and bug reports
=================

Your comments, bug reports, and suggestions are very welcomed.
They will help us to further improve de Bruijn graph construction tool.
If you have any trouble running graph construction, please send us `dbg.log` file from the directory `<output_dir>`.

You can send your comments and bug reports via e-mail: [anton.bankevich@gmail.com](mailto:anton.bankevich@gmail.com).

