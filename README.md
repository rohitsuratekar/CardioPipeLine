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
* [MultiQc](https://multiqc.info/) (Quality Control report)

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
* [pandas](https://pandas.pydata.org/) 
* `R 4+` (in case you want DESeq2 analysis. Tested on R 4.0.2)
* `tximport`, `DESeq2`, `AnnotationDbi`, `Rsubread`, for `R` related analysis

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

    
| Task No | Details | Tool Used | Input <sup>#</sup> | Output<sup> #, * </sup> |
| :------ | :----- | :------ | :----| :--- |
| 0 | All tasks | -- | -- | -- |
| 1 | Download `.sra` file| `prefetch` | samples.csv | srr.sra |
| 2 | Convert to `fastq` | `fasterq-dump` | srr.sra | srr.fastq|
| 3 | rRNA filtering | `sortmerna` | srr.fastq | srr.filtered.fastq |
| 4 | STAR indexing | `star` | genome.fa, anno.gtf | star-idx |
| 5 | STAR mapping  | `star` | star-idx, srr.filtered.fastq / srr.fastq | srr.bam |
| 6 | Quantification | `stringtie` | srr.bam, anno.gtf | st.tsv |
| 7 | Salmon indexing | `salmon` | anno.gtf | salmon-idx |
| 8 | Salmon mapping | `salmon` | salmon-idx, srr.filtered.fastq / srr.fastq | quants.sf |
| 9 | Kallisto indexing | `kallisto` | anno.gtf | kallisto-idx |
| 10 | Kallisto mapping | `kallisto` | kallisto-idx, srr.filtered.fastq / srr.fastq |abundance.tsv |
| 11 | Count Matrix StringTie | `prepDE.py` | st.tsv | stringtie.counts |
| 12 | Count Matrix Star | `Rsubread` | ssr.bam | star.counts |
| 13 | Count Matrix Salmon | `tximport` | quants.sf | salmon.counts |
| 14 | Count Matrix Kallisto | `tximport` | abundance.tsv | kallisto.counts |
| 15 | Differential Analysis | `DESeq2` | *.counts | con1_vs_con2.csv |
| 16 | Quality Control Report | `multiqc` | -- | quality_report.html |


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
 
By default `tasks: 0` is set. This will run all available tasks.

### Setting up the workflow 
Editing `config.yaml` and `samples.csv` is probably the most important
 aspect of this pipeline. These are files which will change for each user
  according to their system and work.

#### config.yaml
 If you are fine  with default options of tools
  used in this workflow, you can just edit this file and run the pipeline
   with `snakemake`. Just read the comment in  `config/config.yaml` file
    and you will know what to do :) Comments are self explanatory.

Config file has parameter called `base` which provides base folder for all
 the analysis. Every output file crated with this pipeline can be found in
  this base folder. Overall output file structure is shown below. More
   detailed file structure can be found [here](https://github.com/rohitsuratekar/CardioPipeLine/wiki/Base-Folder-Structure). 
   
```
base
├── sra
├── fastq
├── filtered
├── bams
├── index
├── mappings
├── deseq2
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

If you are performing any DESeq2 analysis, third column specifying
 conditions is mandatory
 
```
run,is_paired,condition
SRR0000001,true,wt
SRR0000002,false,mt
```
Here third column `condition` will be used in DESeq2 analysis. You should
 update appropriate parameters in `config.yaml` file. For example, Let us
  assume following sample file
  
```
run,is_paired,time
SRR0000001,true,48
SRR0000002,false,24
```
In above file, as third column name is `time` and `24` is your reference
 condition then you should update
 following in the `config.yaml`
```
deseq2:
    design_column: "time"
    reference: "24"
```
In addition, you should also update `design` parameter from `config.yaml
` file. For example, in above situation, it can be `design: "~ time"`

