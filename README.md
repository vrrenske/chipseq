# chipseq analysis pipeline

This analysis pipeline is based on [Jasper Stinenbosch' R pipeline](https://github.com/JStinenbosch). To increase computing speed, I rewrote the pipeline in Snakemake (python-based).
I added some basic visualization using R (ggplot2 and plotly), which gives an overview of the relative sequencing depth per 200 bps. I plan to further automate the pipeline (add dependencies per rule) & add more automated visualization,
but for now, there are some things which should be installed before running.

## Operating system
I recommend running the whole pipeline in a linux environment, since some programs (for instance Macs2) don't run in Windows. I think OSx might work, but I didn't test this!

If you want to run this pipeline on a Windows computer, install the [Ubuntu app](https://www.microsoft.com/en-us/p/ubuntu/9nblggh4msv6?activetab=pivot:overviewtab) on your computer. From there, you can acces your C drive through:

```shell
cd /mnt/c
```

## Dependencies

This pipeline uses the following programs. I installed all of them using [Conda](https://bioconda.github.io/), which makes life a lot easier. To install conda, follow the instructions [here](https://bioconda.github.io/user/install.html#install-conda).

To install a program, use the following command:

```shell
conda install theprogramyouwanttoinstall
```
### Programs

* [Snakemake](https://snakemake.readthedocs.io/en/stable/)
* [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
* [samtools](http://www.htslib.org/)
* [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)

The pipeline also uses [Macs2](https://github.com/taoliu/MACS) to call for peaks, but the snakemake pipeline automatically uses conda to create a seperate environment to use Macs2. This is because Macs2 needs Python 2.7, which is incompatible with the Python 3.7 environment the rest of the program was written in.

### Programming languages

* [Python 3.7](https://www.python.org/downloads/release/python-370/)
* [R 3.5](https://cran.rstudio.com/)


## How to run the pipeline

The easiest is to copy-paste the whole set of files in this repository to the folder where you want to analyze your files. In the same folder, create a folder named "Test_folder" where you put all your "fastq.gz" files. To run the whole pipeline, run:

```shell

snakemake -s chip_pipeline.smk --use-Conda

```

To run only one of the rules in the pipeline, put the name of the rule behind the command:

```shell

snakemake -s chip_pipeline.smk map_compress

```

## Specifics of the pipeline you might want to change

1. The target folder is named Test_folder. If you want to change this, do so in both "chip_pipeline.smk" and "python_orderfiles_snake.py".
2. Check whether the trimmomatic parameters are fit for your organism/experiment
3. The genome size, mapping & plotting are all based on the D39V *Streptococcus pneumoniae* genome. Change the size, FASTA & gb files,
and referrals in "chip_pipeline.smk" & "convertBAMs.R" to change this.
4. The mapping, trimming & peak calling are counting on *paired end* sequencing.
5. Check your expected fragment size, library type etc.
6. The MACS2 commands are set to find *narrow peaks*. Add `-broad` to the commands to change.
