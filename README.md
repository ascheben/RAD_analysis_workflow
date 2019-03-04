# RAD_analysis_protocol
A bioinformatics protocol for analysing genotypes using restriction site associated DNA sequencing data

## Table of contents
I. [Introduction](#Introduction)  
II. [Software requirements](#Software-requirements)  
III. [SNP calling protocol](#SNP-calling-protocol)  
IV. [Diversity analysis protocol](#Diversity-analysis-protocol)  


## Introduction <a name="Introduction"></a>  
XXXX Cite PBJ paper and talk about diff.

## Software requirements <a name="Software-requirements"></a>
* Stacks : http://catchenlab.life.illinois.edu/stacks/
* trimmomatic : http://www.usadellab.org/cms/?page=trimmomatic
* fastqc : https://www.bioinformatics.babraham.ac.uk/projects/download.html
* multiqc : https://multiqc.info/
* mash : https://github.com/marbl/Mash
* BWA : https://sourceforge.net/projects/bio-bwa/files/    
* SAMtools : http://www.htslib.org/download/        
* BCFtools : http://www.htslib.org/download/  
* VCFtools : https://github.com/vcftools/vcftools
* SNPRelate : https://github.com/zhengxwen/SNPRelate
* plink : http://zzz.bwh.harvard.edu/plink/download.shtml
* fastStructure : https://rajanil.github.io/fastStructure/
* pophelper : https://github.com/royfrancis/pophelper
* vcf2phy : https://github.com/edgardomortiz/vcf2phylip
* RAxML : https://cme.h-its.org/exelixis/software.html
* ggtree : https://github.com/GuangchuangYu/ggtree

I recommend installing these tools using conda (XXX). In my experience some incompatibilites between tools cannot be resolved, therefore it may be necessary to create multiple conda environments.
`conda create`

## SNP calling protocol <a name="SNP-calling-protocol"></a>  
XXX
### Demultiplexing and adapter trimming
``bcl2fastq``
``process_radtags``
``trimmomatic``

### Quality control of reads
``fastqc``
``multiqc``
``mash``

### De novo and reference based SNP calling
XXX

#### De novo assembly and SNP calling
Link to stacks homepage. Refer to manual pipeline for large datasets. Note that cstacks is most computationally intensive step and can be run in batches and results combined.

Assemble and call
``denovo_map.pl``
Export vcf
``populations``

Discuss crucial parameters for denovo_map.pl and link to Paris paper to show how to best explore parameter space for optimized denovo assembly of stacks.

#### Mapping reads to a reference and SNP calling

Include references that suggest bwa-mem + samtools/bcftools is most effective

``bwa-mem`` 
``gstacks``
``samtools``
``bcftools``

### Filtering SNPs

``vcftools`` 

## Diversity analysis protocol <a name="Diversity-analysis-protocol"></a>  
XXX

### Population genetic summary statistics
``populations`` 

### Principal component analysis

``Rscript``

### Structure analysis
``plink`` 
``structure.py`` 

### Inferring a maximum likelihood phylogenetic tree
``filtHet.py``
``vcf2phylip.py``
``RAxML``
``Rscript``
