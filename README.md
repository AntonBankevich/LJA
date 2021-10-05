Tool for de Bruijn graph construction
==============

### Version: 0.1

La Jolla Assembler(LJA) is a tool for genome assembly from HiFI reads based on de Bruijn graphs.
LJA uses very high values of k for de Bruijn graph construction thus automatically resolving almost all repeats even in mammalian genomes.
LJA consists of three modules each of which addresses a computational challenge that us unique for HiFI reads:
(i) jumboDBG module for de Bruijn graph construction for arbitrarily large values of k;
(ii) mowerDBG module for almost perfect error correction and
(iii) multiDBG module for using variable values of k in different parts of the genome and using the full length of HiFi reads for repeat resolution.
In addition, since LJA compresses homopolymers in reads to avoid frequent errors in homopolymers, LJA uses LJApolisher tool for uncompressing and polishing final contigs.
For more details please refer to our [paper](https://www.biorxiv.org/content/10.1101/2020.12.10.420448).

Please note that LJA software is still a work in progress.
In current version diploid assembly is an experimental feature and for now we can not combine HiFI reads with other technologies.
In addition, LJA uses more memory and works slower than intended but still better than most other tools.
We constantly work on improving LJA performance and results.
Please check for the new versions regularly and send any questions and bug reports to [anton.bankevich@gmail.com](mailto:anton.bankevich@gmail.com).  


For LJA installation and running instructions please refer to [LJA manual](docs/lja_manual.md).
We also provide jumboDBG module for de Bruijn graph construction as a separate script.
For jumboDB running instructions please refer to [jumboDBG manual](docs/jumbodbg_manual.md).

License
-------

This tool is distributed under a BSD license. See the [LICENSE file](LICENSE) for details.


Credits
-------
If you use this software in your research please cite our [paper](https://www.biorxiv.org/content/10.1101/2020.12.10.420448).
LJA is developed by Anton Bankevich and Andrey Bzikadze in [Pavel Pevzner's lab at UCSD](http://cseweb.ucsd.edu/~ppevzner/)
and Dmitry Antipov in [Center for Algorithmic Biology at SPbSU](https://cab.spbu.ru/).
