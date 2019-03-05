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

``
[Data]  
FCID,Lane,SampleID,SampleName,Index,Description,Control,Recipe,Operator,SampleProject
HH3LTCCXY,7,1_FD01,Other,,external_id:ddRAD1,N,,SIX,R_XXX_M002``  

The fastq files will almost always contain multiple pooled (i.e., multiplexed) samples. There are many tools that can separate out (i.e., demultiplex) the individual samples. Here, for our RAD data we use ``process_radtags``, which is part of the stacks package and was developed especially for demultiplexing RRS read libraries. We will need to provide a text file containing the barcodes for each of our samples. Here is an example of the format expected by ``process_radtags``:  

``
AACCA	sample_01  
AAGGA	sample_02  
AATTA	sample_03  
``

An advantage of ``process_radtags`` is that unlike most demultiplexing tools, it can use information on restriction cut sites to quality control the reads using the option ``-e`` or , for ddRAD-seq, ``--renz_1`` and ``--renz_2``. Depending on the library preparation technique, reads may include the restriction site(s) targeted by the enzyme(s) used. A read missing the restriction site may be the result of technical errors. These reads should be discarded, particularly when you plan to de novo assemble your reads, because the assembly method we will use requires uniform reads. The rescue option ``-r`` will attempt to rescue restriction sites and barcodes if they have a minor mismatch with the expected sequence. In addition, we can discard reads containing unknown nucleotides (Ns) and low per base quality scores using the ``-c`` and ``-q`` options respectively. Additional information on options and a list of supported enzymes can be found [here](http://catchenlab.life.illinois.edu/stacks/comp/process_radtags.php).  To demuliplex our paired-end sequences we can run the program like so:  

``process_radtags -i gzfastq -P -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz -b ./barcodes.txt -o . -r -c -q --renz_1 pstI --renz_2 nlaIII``  

Demultiplexed reads may contain adapter contamination, which can hinder read alignment and assembly. Reads containing adapter sequences should therefore be discarded for de novo assembly using stacks. This can be achieved in a range of ways; here will do it using ``trimmomatic``, by using the provided paired-end Illumina adapters and setting the ``MINLEN`` option to the read length.  

``trimmomatic sample_01_R1.fq sample_01_R2.fq sample_01_pe_R1.fq sample_01_se_R1.fq sample_01_pe_R2.fq sample_01_se_R2.fq ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 MINLEN:120``  

If you are aligning reads to a reference genome, then you can remove the above ``MINLEN`` option, and adapters will simply be trimmed off the read. Tools like ``trimmomatic`` also allow trimming based on quality, however, de novo assembly with stacks relies on uniform read lengths and aggressive quality trimming can also reduce read alignment to a reference genome. Therefore quality trimming is not recommended.

### Quality control of reads
Now that we have demultiplexed and filtered our samples, we should review read quality statistics and estimate the genetic similarity between samples. This will allow us to spot issues with the data immediately, before proceeding to more computationally intensive steps of the analysis. Using the tools ``fastqc`` and ``multiqc`` we can generate a single html quality report for our samples. To generate an individual quality report per sample, we first run ``fastqc`` in each fastq (or gzipped fastq) file in your directory of demultiplexed, adapter-trimmed samples.

``for sample in *.fq; do fastqc ${sample};done``  

<div class="alert alert-notice">
Note that as with most analyses in this protocol ``fastqc`` runtime scales linearly with sample number and sample size. For large sample sizes (>100) and large read numbers (>1 Gigabases / sample), parallel execution of commands on a high-performance computing resource is suggested.
</div>

Inspecting each quality report per sample is difficult for large sample numbers, therefore we use ``multiqc`` to generate a single master report from the individual ``fastqc`` reports.

``multiqc /path/to/fastqc/outputs/``

The most important quality parameters for our data are the quality scores, GC content and the adapter content. If samples fail these tests in multiqc, additional filtering may be required. Depending on the study organism, deviation from the expected GC content may indictae bacterial contamination, as bacteria have a much wider range of genomic GC content than most plants and animals. In contrast, it is alright for RAD samples to fail 'Per Base Sequence Content' tests, if the restriction site has not been removed from the samples. Because R1, R2 or both R1 and R2 reads are expected to begin with a restriction site, this is not an issue, and we can expect a pattern like the one from real paired end ddRAD-seq samples digested with pstI and nlaIII shown below.

![alt text](https://github.com/ascheben/RAD_analysis_protocol/raw/master/src/images/mqc_pbsq.png "Multiqc Per Base Sequence Content Report for ddRAD samples")  

Duplicate reads are also expected in ddRAD-seq by design, as both the end and the start of the sequence is determined by a restriction site.  
A final quality control step consists of a k-mer based analysis of genetic distances between samples. This can help identify contamination and mislabeling. Mash provides a rapid algorithm for estimating genetic distances between samples using sets of reads. First a sketch is produced for each sample (this can take a while depending on sample size).  

``for sample in *.fq; do mash sketch -r -m 4 -s 1000 ${sample};done``  

The minimum copies of each k-mer (``-m``) can be increased to ensure rare and therefore possibly spurious k-mers are excluded. Increasing the sketch size (``-s``) can increase the accuracy of the distance estimation.  

Next the distance between each pair of sketches is estimated and written to an output file (this is fast).

``
for sample_x in *.msh; do
    for sample_y in *.msh; do
        mash  dist ${sample_x} ${sample_y}
    done
done > All_Distances.txt``

We can check the p-values in the output file to ensure they are significant, and then we can remove the superfluous columns:

``cut -f1-3 All_Distances.txt > Distances.txt``

Finally we can use R to plot a tree based on the mash distance matrix.

``
library("phangorn") # may have to install phangorn first  
a <- read.table("Distances.txt", stringsAsFactors=F, sep="\t")  
matrix <- reshape(a, direction="wide", idvar="V2", timevar="V1")  
distance <- as.dist(matrix[,-1], upper=F, diag=F)  
attr(distance, "Labels") <- matrix[,1]  
plot(upgma(distance),cex = 0.5)  
add.scale.bar(ask = TRUE)  
``

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
