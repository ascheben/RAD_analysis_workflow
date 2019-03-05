# RAD_analysis_protocol
A bioinformatics protocol for analysing genotypes using restriction site associated DNA sequencing data

## Table of contents
I. [Introduction](#Introduction)  
II. [Software requirements](#Software-requirements)  
III. [SNP calling protocol](#SNP-calling-protocol)  
IV. [Diversity analysis protocol](#Diversity-analysis-protocol)  
V. [References](#References)

## Introduction <a name="Introduction"></a>  
Reduced representation sequencing (RRS) encompasses a suite of methods for sampling a limited number of regions from a genome for DNA sequencing (Scheben et al., 2017). The main advantage of RRS over whole genome sequencing is that it substantially reduces the cost of generating sequencing data with adequate coverage for accurately calling single nucleotide polymorphisms (SNPs). Several RRS approaches including RAD-seq (Miller et al., 2007; Baird et al., 2008) and ddRAD-seq (Peterson et al., 2012) use restriction enzymes to digest the genome. DNA fragments within a suitable size range are then sequenced. In this bioinformatics protocol we will go through the most common steps in the ddRAD-seq analysis pipeline, starting with raw reads and concluding with an exploration of genetic structure and diversity. 

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

I recommend installing these tools using the package manager [conda](https://conda.io/en/latest/) in a Linux environment. On Windows systems, it is possible to use the [Windows Subsystem for Linux](https://docs.microsoft.com/en-us/windows/wsl/about), rather than using a virtual machine. In my experience some incompatibilites between tools cannot be resolved, therefore it may be necessary to create multiple conda environments for different tasks. For example, to create an evironment for demultiplexing and quality control, one could use:  

`conda create -n myenv -c bioconda python=3.5 stacks trimmomatic fastqc multiqc mash`  

The environment can then be loaded at any time with:  

`conda activate myenv`

## SNP calling protocol <a name="SNP-calling-protocol"></a>  
SNP calling is simple from a users perspective when a high-quality reference genome is available, as is the case for most model organisms. In non-model organisms without reference genomes, reads are generally assembled de novo in order to call SNPs. De novo assembly remains an error-prone task and therefore, as a general rule, reference-based SNP calling is preferred.
### Demultiplexing and adapter trimming
Illumina sequencing providers often upload reads in fastq format, although they may also provide reads in the raw Illumina base call (bcl) format. For converting bcl to fastq, the [bcl2fastq](http://sapac.support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software.html) pipeline can be used. As its main input, the pipeline takes the path to the bcl run and the sample sheet made available by the sequencing provider. To convert the bcl data to fastq using 8 threads and default parameters, we can use:  

``bcl2fastq --runfolder-dir /path/to/run/dir --sample-sheet=./mysamplesheet.csv -o /path/to/out/dir -r 8 -p 8 -w 8``  

Below is an example of a real sample sheet. The most important information is the lane number and the IDs. If no sample sheet is provided with your run, you can create a mock sample sheet with this information and use that as input for ``bcl2fastq``.  

``[Data]
FCID,Lane,SampleID,SampleName,Index,Description,Control,Recipe,Operator,SampleProject
HH3LTCCXY,7,1_FD01,Other,,external_id:ddRAD1;family_id:;pool_name:CST Batch 1_FD01-1;pool_limsid:2-1;guess_library_type:None;guess_library_name:None;guess_library_limsid:None;l_0:lib#180625_001_P_DNF lanes 1-2 & 7-8_G1|limsid#2-1|proc#CreateProductionCSTBatch|prot#CSTCreation;l_1:lib#LP7-NTP_D2|limsid#2-1|proc#AUTOMATED-MakeNTP|prot#SubmittedLibraryQC;l_2:lib#LP02-DCT_D2|limsid#2-1|proc#AUTOMATED-MakeDCT|prot#SubmittedLibraryQC;l_3:lib#LP8-TSP1_D2|limsid#2-1|proc#MakeTSP1|prot#SubmittedLibraryQC;l_4:lib#lib_B1|limsid#2-1|proc#CreateProductionSeqLabBatch|prot#SampleAccessioning&InitialQC;l_5:lib#FD01_11|limsid#XXX|proc#|prot#,N,,SIX,R_XXX_M002``  

The fastq files will almost always contain multiple pooled (i.e., multiplexed) samples. There are many tools that can separate out (i.e., demultiplex) the individual samples. Here, for our RAD data we use ``process_radtags``, which is part of the stacks package and was developed especially for demultiplexing RRS read libraries. We will need to provide a text file containing the barcodes for each of our samples. Here is an example of the format expected by ``process_radtags``:  

``AACCA	sample_01
AAGGA	sample_02
AATTA	sample_03
``

An advantage of ``process_radtags`` is that unlike most demultiplexing tools, it can use information on restriction cut sites to quality control the reads using the option ``-e`` or , for ddRAD-seq, ``--renz_1`` and ``--renz_2``. Depending on the library preparation technique, reads may include the restriction site(s) targeted by the enzyme(s) used. A read missing the restriction site may be the result of technical errors. These reads should be discarded, particularly when you plan to de novo assemble your reads, because the assembly method we will use requires uniform reads. The rescue option ``-r`` will attempt to rescue restriction sites and barcodes if they have a minor mismatch with the expected sequence. In addition, we can discard reads containing unknown nucleotides (Ns) and low per base quality scores using the ``-c`` and ``-q`` options respectively. Additional information on options and a list of supported enzymes can be found [here](http://catchenlab.life.illinois.edu/stacks/comp/process_radtags.php).  To demuliplex our paired-end sequences we can run the program like so:  

``process_radtags -i gzfastq -P -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz -b ./barcodes.txt -o . -r -c -q --renz_1 pstI --renz_2 nlaIII``  

Demultiplexed reads may contain adapter contamination, which can hinder read alignment and assembly. Reads containing adapter sequences should therefore be discarded for de novo assembly using stacks. This can be achieved in a range of ways; here will do it using ``trimmomatic``, by using the provided paired-end Illumina adapters and setting the ``MINLEN`` option to the read length.  

``trimmomatic sample_01_R1.fq sample_01_R2.fq sample_01_pe_R1.fq sample_01_se_R1.fq sample_01_pe_R2.fq sample_01_se_R2.fq ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 MINLEN:120``  

If you are aligning reads to a reference genome, then you can remove the above ``MINLEN`` option, and adapters will simply be trimmed off the read. Tools like ``trimmomatic`` also allow trimming based on quality, however, de novo assembly with stacks relies on uniform read lengths and aggressive quality trimming can also reduce read alignment to a reference genome. Therefore quality trimming is not recommended.

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

## References <a name="References"></a>  
