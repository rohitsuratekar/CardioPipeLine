# CardioPipeLine

This is simple utility to analyse the RNA-seq data with variety of popular
 tools. It is written in `snakemake` which allows you to sclae it
  to any system including clusters and cloud computing architectures.
   Currently it provides analysis with following tools

* [STAR](https://github.com/alexdobin/STAR) (genome based mapping)
* [Salmon](https://salmon.readthedocs.io/en/latest/salmon.html) (transcript based mapping)
* [Kallisto](https://github.com/pachterlab/kallisto) (transcript based  mapping)
* [SortMeRNA](https://github.com/biocore/sortmerna) (rRNA filtering)
* [StringTie](https://ccb.jhu.edu/software/stringtie/) (quantification of BAM files)
* [NCBI-tools](https://github.com/ncbi/sra-tools) (downloading and conversion from NCBI server)
* [DEseq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) (differential expression)

In short, given public run ID from SRA archives (SRRXXXXXXX) it can download
 it from NCBI server, converts to fastq, filters the read, creates indices, maps to
  genome/transcript and generates quantification. In addition it can perform
  differential expression analysis.

However, you can use your own FASTA or SRA files for any of the analysis.
In order to change the default options, you will need to change appropriate
 rules in respective `rules/filename.smk` files.
 
### Dependencies
* `python3+` (tested on Python3.7 and 3.8)
* [snakemake](https://snakemake.readthedocs.io/en/stable/) (worlflow
 management)
* [pyyaml](https://pyyaml.org/) (To read the configuration file)
* [pandas](https://pandas.pydata.org/) (Was needed only for one function but
 I am too lazy to implement that function with in-build libraries.)
* `R 4+` (in case you want DESeq2 analysis. Tested on R 4.0.2)
* `tximport`, `DESeq2`, `AnnotationDbi`, `Rsubread`, `yaml` for `R` related analysis

All above dependencies have been successfully tested on Ubuntu 18.04.5
, Gentoo 5.4.48, Fedora 32.

### Installation
On installation part, it is good idea to make a new virtual environment and
 install `snakemake` on it. Then download all the tool binaries in your
  desired place. After that, just clone this repository.

```
git clone https://github.com/rohitsuratekar/CardioPipeLine.git
```
If you do not have `git`, just download the repository and unzip it in your
 desired folder.

After that only thing you need to do is change `config/config.yaml` and
 `samples.csv` file according to your system. and then run following in the
  repository folder

`
snakemake --core all
`
Above command will essentially perform everything for you, from downloading
 file from NCBI to making quantification :)
 
### Tasks
This pipeline performs tasks in following order. However, you can always
 select which tasks you want :) 

    
| Task No | Details | Tool Used | Input <sup>#</sup> | Main Output<sup> #, * </sup> |
| :------ | :----- | :------ | :----| :--- |
| 0 | All tasks | -- | -- | -- |
| 1 | Download `.sra` file| `prefetch` | samples.csv | srr.sra |
| 2 | Convert to `fastq` | `fasterq-dump` | srr.sra | srr.fastq|
| 3 | rRNA filtering | `sortmerna` | srr.fastq | srr.filtered.fastq |
| 4 | STAR indexing | `star` | genome.fa, anno.gtf | star-idx |
| 5 | STAR mapping  | `star` | star-idx, srr.filtered.fastq / srr.fastq | srr.bam |
| 6 | Quantification | `stringtie` | srr.bam, anno.gtf | srr.tsv |
| 7 | Salmon indexing | `salmon` | anno.gtf | salmon-idx |
| 8 | Salmon mapping | `salmon` | salmon-idx, srr.filtered.fastq / srr.fastq | quants.tsv |
| 9 | Kallisto indexing | `kallisto` | anno.gtf | kallisto-idx |
| 10 | Kallisto mapping | `kallisto` | kallisto-idx, srr.filtered.fastq / srr.fastq | mapping.tsv |
| 11 | Count Matrix generation | `Rsubread` | ssr.bam | counts.csv |
| 12 | Count Matrix generation | `tximport` | quants.tsv / mapping.tsv | counts.csv |
| 13 | Differential Analysis | `DESeq2` | counts.csv | exp.csv |
| 14 | Clean up | `shell` | -- | -- |


<sup>* There might be many other related outputs. </sup> 
<sup># Names are for representative purpose. Actual names will be different
 </sup>

#### Selection of tasks
You have full control of which tasks to run. Your can run desired tasks by
 editing `tasks` parameter in the `config.yaml` file. If you want to perform
  specific task, you need to provide its input file in appropriate location.  

But do not worry, thanks to `snakemake`, if you did not provide appropriate
 input files, it will  search far task which outputs these files and will
  automatically run for you. 
 
By default `tasks: 15` is set. This will run all available tasks.

### Setting up the workflow 
Editing `config.yaml` and `samples.csv` is probably the most important
 aspect of this pipeline. These are files which will change for each user
  according to their system, working 

#### config.yaml
 If you are fine  with default options of tools
  used in this workflow, you can just edit this file and run the pipeline
   with `snakemake`. Just read the comment in  `config/config.yaml` file
    and you will know what to do :) Comments are self explanatory.

Config file has parameter called `base` which provides base folder for all
 the analysis. Every output file crated with this pipeline can be found in
  this base folder. Let us assume we are analysing `SRR0000001` (paired end
  ) and `SRR0000002` (single end). We will get following folder structure
   after full analysis
   
```
base
├── sra
│   ├── SRR0000001.sra
│   └── SRR0000002.sra
├── fastq
│   ├── SRR0000001.sra_1.fastq
│   ├── SRR0000001.sra_2.fastq
|   └── SRR0000002.sra.fastq
├── filtered
│   ├── SRR0000001.sra.filtered_1.fastq
│   ├── SRR0000001.sra.filtered_2.fastq
|   └── SRR0000002.sra.filtered.fastq
├── bams
│   ├── SRR0000001.sra.Aligned.sortedByCoord.out.bam
|   └── SRR0000002.sra.Aligned.sortedByCoord.out.bam
├── index
│   ├── star / .. 
|   ├── sortmerna / ..
|   ├── salmon / ..
|   └── kallisto / ..
├── mappings
│   ├── star / .. 
|   ├── stringtie / ..
|   ├── salmon / ..
|   └── kallisto / ..
├── deseq2
│   ├── star / .. 
|   ├── stringtie / ..
|   ├── salmon / ..
|   └── kallisto / ..
└── logs
```

#### samples.csv
This file will contain names of all the samples which you want to analyse
 with this pipeline. You should provide SRA **Runs** IDs (usually in
  *SRRXXXXXXX* format ). You can find these IDs for publically available
   datasets from [GEO Dataset](https://www.ncbi.nlm.nih.gov/gds) website. If
    you have your own *fastq* files, you can check naming convention and
     rename those files and keep in appropriate folder. and start the
      pipeline.
 
 This `.csv` (comma delimited) file SHOULD have header row and SHOULD
  contains at least following headers.

* `run` : Run ID
* `is_paired` : true (if you are handling paired end sample), false (if it
 is single end sample)

Example file content
```
run,is_paired
SRR0000001,true
SRR0000002,false
```
