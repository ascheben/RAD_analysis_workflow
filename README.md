# RAD_analysis_protocol
A bioinformatics protocol for analysing genotypes using restriction site associated DNA sequencing data

## Table of contents
I. [Introduction](#Introduction)  
II. [Software requirements](#Software-requirements)  
III. [SNP calling protocol](#SNP-calling-protocol)  
- [Demultiplexing and adapter trimming](#Demultiplexing-and-adapter-trimming)  
- [Quality control of reads](#Quality-control-of-reads)  
- [_De novo_ assembly and SNP calling](#De-novo-assembly-and-SNP-calling)  
- [Mapping reads to a reference and SNP calling](#Mapping-reads-to-a-reference-and-SNP-calling)  
- [Filtering SNPs](#Filtering-SNPs)  

IV. [Diversity analysis protocol](#Diversity-analysis-protocol)  
- [Population genetic summary statistics](#Population-genetic-summary-statistics)
- [Principal component analysis](#PCA)  
- [Structure analysis](#Structure-analysis)  
- [Inferring a maximum likelihood phylogenetic tree](#Tree-building)  

V. [References](#References)

## Introduction <a name="Introduction"></a>  
Reduced representation sequencing (RRS) encompasses a suite of methods for sampling a limited number of regions from a genome for DNA sequencing ([Scheben et al., 2017](https://onlinelibrary.wiley.com/doi/full/10.1111/pbi.12645)). The main advantage of RRS over whole genome sequencing is that it substantially reduces the cost of generating sequencing data with adequate coverage for accurately calling single nucleotide polymorphisms (SNPs). Several RRS approaches including RAD-seq ([Miller et al., 2007](https://genome.cshlp.org/content/17/2/240.short); [Baird et al., 2008](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0003376)) and ddRAD-seq ([Peterson et al., 2012](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0037135)) use restriction enzymes to digest the genome. DNA fragments within a suitable size range are then sequenced. In this bioinformatics protocol we will go through the most common steps in the ddRAD-seq analysis pipeline, starting with raw reads and concluding with an exploration of genetic structure and diversity. 

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
SNP calling is simple from a users perspective when a high-quality reference genome is available, as is the case for most model organisms. In non-model organisms without reference genomes, reads are generally assembled _de novo_ in order to call SNPs. _De novo_ assembly remains an error-prone task and therefore, as a general rule, reference-based SNP calling is preferred.
### Demultiplexing and adapter trimming <a name="Demultiplexing-and-adapter-trimming"></a> 
Illumina sequencing providers often upload reads in fastq format, although they may also provide reads in the raw Illumina base call (bcl) format. For converting bcl to fastq, the [bcl2fastq](http://sapac.support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software.html) pipeline can be used. As its main input, the pipeline takes the path to the bcl run and the sample sheet made available by the sequencing provider. To convert the bcl data to fastq using 8 threads and default parameters, we can use:  

``bcl2fastq --runfolder-dir /path/to/run/dir --sample-sheet=./mysamplesheet.csv -o /path/to/out/dir -r 8 -p 8 -w 8``  

Below is an example of a real sample sheet. The most important information is the lane number and the IDs. If no sample sheet is provided with your run, you can create a mock sample sheet with this information and use that as input for ``bcl2fastq``.  

``
[Data]  
FCID,Lane,SampleID,SampleName,Index,Description,Control,Recipe,Operator,SampleProject
HH3LTCCXY,7,1_FD01,Other,,external_id:ddRAD1,N,,SIX,R_XXX_M002``  

The fastq files will almost always contain multiple pooled (i.e., multiplexed) samples. There are many tools that can separate out (i.e., demultiplex) the individual samples. Here, for our RAD data we use ``process_radtags``, which is part of the stacks package and was developed especially for demultiplexing RRS read libraries. We will need to provide a text file containing the barcodes for each of our samples. Here is an example of the format expected by ``process_radtags``:  

```
AACCA	sample_01  
AAGGA	sample_02  
AATTA	sample_03  
```

An advantage of ``process_radtags`` is that unlike most demultiplexing tools, it can use information on restriction cut sites to quality control the reads using the option ``-e`` or , for ddRAD-seq, ``--renz_1`` and ``--renz_2``. Depending on the library preparation technique, reads may include the restriction site(s) targeted by the enzyme(s) used. A read missing the restriction site may be the result of technical errors. These reads should be discarded, particularly when you plan to _de novo_ assemble your reads, because the assembly method we will use requires uniform reads. The rescue option ``-r`` will attempt to rescue restriction sites and barcodes if they have a minor mismatch with the expected sequence. In addition, we can discard reads containing unknown nucleotides (Ns) and low per base quality scores using the ``-c`` and ``-q`` options respectively. Additional information on options and a list of supported enzymes can be found [here](http://catchenlab.life.illinois.edu/stacks/comp/process_radtags.php).  To demuliplex our paired-end sequences we can run the program like so:  

``process_radtags -i gzfastq -P -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz -b ./barcodes.txt -o . -r -c -q --renz_1 pstI --renz_2 nlaIII``  

Demultiplexed reads may contain adapter contamination, which can hinder read alignment and assembly. Reads containing adapter sequences should therefore be discarded for _de novo_ assembly using stacks. This can be achieved in a range of ways; here will do it using ``trimmomatic``, by using the provided paired-end Illumina adapters and setting the ``MINLEN`` option to the read length.  

``trimmomatic sample_01.1.fq sample_01.2.fq sample_01_pe_R1.fq sample_01_se_R1.fq sample_01_pe_R2.fq sample_01_se_R2.fq ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 MINLEN:120``  

For further analysis using ``stacks``, samples should be renamed using the expected naming scheme:  

```
# R1 read
<sample_name>.1.fq 
# R2 read
<sample_name>.2.fq 
```

If you are aligning reads to a reference genome, then you can remove the above ``MINLEN`` option, and adapters will simply be trimmed off the read. Tools like ``trimmomatic`` also allow trimming based on quality, however, _de novo_ assembly with stacks relies on uniform read lengths and aggressive quality trimming can also reduce read alignment to a reference genome. Therefore quality trimming is not recommended.

When PCR is used to amplify DNA before sequencing, PCR duplicates can occur in reads and inflate coverage estimation. As genotype confidence is usually based on read depth, removing PCR duplicates helps increase the accuracy of genotype confidence estimation. Clones in RAD-seq data can be removed using a ``stacks`` tool by identifying identical reads in a sample.

``clone_filter -1 ./sample_01.1.fq.gz -2 ./sample_01.2.fq.gz -i gzfastq -o ./clones_removed/``

However, because all sequence reads from a ddRAD-seq locus are identical, clones cannot (!) be removed from ddRAD-seq if no random oligo tags were introduced into reads during molecular library construction (e.g., [Hoffberg et al., 2016](https://onlinelibrary.wiley.com/doi/full/10.1111/1755-0998.12566). When ddRAD-seq library construction uses a technique to mark duplicates, ``clone_filter`` can be used to remove these.

### Quality control of reads <a name="Quality-control-of-reads"></a> 
Now that we have demultiplexed and filtered our samples, we should review read quality statistics and estimate the genetic similarity between samples. This will allow us to spot issues with the data immediately, before proceeding to more computationally intensive steps of the analysis. Using the tools ``fastqc`` and ``multiqc`` we can generate a single html quality report for our samples. To generate an individual quality report per sample, we first run ``fastqc`` in each fastq (or gzipped fastq) file in your directory of demultiplexed, adapter-trimmed samples.

``for sample in *.fq; do fastqc ${sample};done``  

Note that as with most analyses in this protocol ``fastqc`` runtime scales linearly with sample number and sample size. For large sample sizes (>100) and large read numbers (>1 Gigabases / sample), parallel execution of commands on a high-performance computing resource is suggested.

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

```
library("phangorn") # may have to install phangorn first  
a <- read.table("Distances.txt", stringsAsFactors=F, sep="\t")  
matrix <- reshape(a, direction="wide", idvar="V2", timevar="V1")  
distance <- as.dist(matrix[,-1], upper=F, diag=F)  
attr(distance, "Labels") <- matrix[,1]  
plot(upgma(distance),cex = 0.5)  
add.scale.bar(ask = TRUE)  
```

### _De novo_ assembly and SNP calling <a name="De-novo-assembly-and-SNP-calling"></a> 
Reads are assmebled and SNPs are called using ``stacks`` in six consecutive steps. Firstly reads per samples are clustered into unique stacks with ``ustacks``. Although we are using paired-end reads, the R2 reads will only be integrated at the tsv2bam stage, so the firest three step of the analysis are carried out without the R2 reads. Below an example is shown for a single file (note that the value for ``-i`` should be incremented for each additional file).

``
ustacks -o . -m 3 -M 1 -p 1 -t gzfastq -f sample_01.1.fq.gz --name sample_01 -i 1
``
A popmap file must be provided for the further analysis steps. This tab-separated file contains one row per sample in the dataset,  with on column for sample names and one for the corresponding population, for example:

```
sample_01   pop1
sample_02   pop1
sample_03   pop2
sample_04   pop2
```

The number of populations in the popmap file can range from one to many. In the next four steps, a catalogue of stacks (i.e., RAD loci) for the complete data set is constructed and genotypes are called. R2 reads in the data directory are automatically integrated into the corresponding R1 RAD locus based on sample names.  

```
popmap="/path/to/popmap.txt"
datadir="/path/to/data/dir"
threads="24"

cstacks -n 1 -P ${datadir} -M ${popmap} -p ${threads}  
sstacks -P ${datadir} -M ${popmap} -p ${threads}  
tsv2bam -P ${datadir} --pe-reads-dir ${datadir} -M ${popmap} -t ${threads}  
gstacks -P ${datadir} -M ${popmap} -t ${threads}  
```

The ``cstacks`` step is generally the computationally intensive step. Newer version of stacks can use existing catalogues as input with the ``-c`` option, which allows the catalogue construction step to be split into multiple consecutive steps if users are limited by walltime.  

The key parameters for _de novo_ assembly of RAD reads are ``-M`` (maximum distance in nucleotides allowed between stacks during constuction of unique stacks per sample) and ``-n`` (number of mismatches allowed between sample loci during catalogue construction). It can be useful to explore the parameter space for optimized denovo assembly of RAD loci as discussed in detail in [Paris et al., 2017](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12775). As a general rule, it is recommended to use the same value for ``-n`` and ``-M``. The larger the value for these options, the more SNPs are generally identified. The optimal values for these key parameters depend largely on the genetic diversity of the sampled population, with more diverse populations requiring larger values. In less diverse populations larger distances allowed between stacks will lead to spurious SNP calls.

Genotype data can be exported in various standard formats including vcf format. To export a vcf file, we can use the following command.

``populations -P ${datadir} -M ${popmap} -t ${threads} --vcf``  

Further details can be found in the [stacks manual](http://catchenlab.life.illinois.edu/stacks/manual/#phand).

### Mapping reads to a reference and SNP calling <a name="Mapping-reads-to-a-reference-and-SNP-calling"></a> 

[A recent comparison](https://link.springer.com/article/10.1186/s12859-017-2000-6) of standard SNP calling methods suggested that bwa-mem alignment followed by SNP calling with samtools/bcftools is the most broadly accurate method. Alternative methods include [freebayes](https://github.com/ekg/freebayes), [GATK](https://software.broadinstitute.org/gatk/) and stacks. Although GATK offers a more complex overall algorithm than samtools, including local realignment of reads, it may not always provide more accurate results, and it is geared toward human data. A major advantage of samtools over GATK is its higher speed and ease of use. 

Reads are aligned to a reference genome using ``bwa mem``, sorted and converted from SAM to BAM format before indexing.

```
bwa mem reference.fa -t 1 -M -R '@RG\tID:sample_01\tPL:illumina\tPU:sample_01\tSM:sample_01' sample_01.1.fq sample_01.2.fq | samtools sort | samtools view -bS -h  > sample_01.bam
samtools index sample_01
``` 
Simple alignment statistics per sample can be calculated using ``samtools flagstat``.  

``samtools flagstat sample_01.bam``

Following read alignment, bam files of all samples are merged into a single bam file.

``samtools merge merged.bam *.bam``

Large bam files can be split by chromosome (or scaffold, as the case may be) using [bamtools](https://github.com/pezmaster31/bamtools).  

``bamtools split -in merged.bam -reference``

If the genome consists of many small scaffolds, these can then be merged again into groups of scaffolds using ``samtools merge``. This then allows parallel processing of each chromosome/scaffold to accelerate SNP calling.

For SNP calling, a pileup file is first created using ``samtools mpileup``. Next SNPs are called using bcftools, which can output genotypes in vcf format.
``` 
samtools mpileup -d 1000 -I -go merged.bcf -ugf reference.fa -t DP merged.bam 
bcftools call --threads 4 -mv -O u -o merged.vcf merged.bcf 
```   

### Filtering SNPs <a name="Filtering-SNPs"></a> 

We employ an iterative filtering approach as discussed in [O'Leary et al.](https://onlinelibrary.wiley.com/doi/full/10.1111/mec.14792). By first removing individuals with extremely low genotyping rates, the number of SNPs passing quality filters is increased. Indels are removed and only biallelic SNPs are retained as both are generally rare and more difficult to analyse with standard software. Minor allele frequency (MAF) is further used to filter rare alleles with a threshold generally set between 0.01 and 0.05. The filters that will remove the most SNPs are usually the depth (``--minDP``) and the missingness (``max-missing``) filters. For heterozygous samples, read depths >=5 have been suggested to avoid undercalling of heterozygous genotypes ([Maruki & Lynch, 2017]( http://www.g3journal.org/content/7/5/1393)). Depth and missingness filters should be fine-tuned for each data set to find a balance between the quantity and quality of SNPs. As a final step, individuals with low genotyping rates based on the filtered SNPs are removed. A simple bash script can carry out all of these steps using ``vcftools``:

```
# Take input vcf file name from user input
myvcf=$1  
# Calculate individual missingness
vcftools --vcf $myvcf --missing-indv --out ${myvcf%%.vcf}  
# Identify samples with over 90% genotypes missing
tail -n +2 ${myvcf%%.vcf}.imiss | awk '$5>0.9' | cut -f1 > ${myvcf%%.vcf}.rm  
# Remove samples with high missingness and apply standard filters
vcftools --vcf $myvcf --remove ${myvcf%%.vcf}.rm --maf 0.05 --max-missing 0.8 --remove-indels --max-alleles 2 --min-alleles 2 --minDP 5 --recode --out ${myvcf%%.vcf}_mim09_biallelic_minDP5_mm08_maf005  
# Randomly keep only a single SNP from each RAD locus (always <500bp in size)
vcftools --vcf ${myvcf%%.vcf}_mim09_biallelic_minDP5_mm08_maf005.recode.vcf --thin 500 --recode --out ${myvcf%%.vcf}_mim09_biallelic_minDP5_mm08_maf005_thin500  
# Recalculate individual missingness
vcftools --vcf ${myvcf%%.vcf}_mim09_biallelic_minDP5_mm08_maf005_thin500.recode.vcf --missing-indv --out ${myvcf%%.vcf}_mim09_biallelic_minDP5_mm08_maf005_thin500  
# Identify samples with over 50% genotypes missing
tail -n +2 ${myvcf%%.vcf}_mim09_biallelic_minDP5_mm08_maf005_thin500.imiss | awk '$5>0.5' | cut -f1 >${myvcf%%.vcf}_mim09_biallelic_minDP5_mm08_maf005.rm  
# Remove samples with over 50% missingness
vcftools --vcf ${myvcf%%.vcf}_mim09_biallelic_minDP5_mm08_maf005_thin500.recode.vcf --remove ${myvcf%%.vcf}_mim09_biallelic_minDP5_mm08_maf005.rm --recode --out ${myvcf%%.vcf}_mim09_biallelic_minDP5_mm08_maf005_thin500_mim05  
```   

The script can be run on a vcf file with the command below.

``./filter_vcf.sh my.vcf``

## Diversity analysis protocol <a name="Diversity-analysis-protocol"></a>  

### Population genetic summary statistics <a name="Population-genetic-summary-statistics"></a>  
A range of useful summary statistics, incorporating population information, can be calculated using ``stacks populations``. A popmap file (described above) must be provided. Statistics can be calculated on nucleotide diversity (pi), heterozygosity and common population genetic statistics such as Fst.

``populations -V merged_filtered.vcf -M popmap.txt --fst --fst_correction p_value`` 

### Principal component analysis <a name="PCA"></a>  
Principal component analysis (PCA) provides a straightforward approach to obtain and overview of genetic diversity within the data. Using ``SNPRelate`` and the ``ggrepel`` package, we can carry out the PCA and plot the results for our samples. The popmap information is used to colour samples by population. 

``Rscript scripts/pca.R <input.vcf> <popmap.txt> <out_prefix>``

### Structure analysis <a name="Structure-analysis"></a>  
An approach complementary to PCA is structure analysis using a Bayesian framework. The software ``fastStructure``  can be used to carry out structure analysis on large SNP datasets. Before analysis, the vcf file must be converted to plink format. The ``plink`` tools was developed for human data and is not able to deal with diverse chromosome or scaffold names in the vcf file. Therefore, for the purpose of the structure analysis, the scaffold names in the vcf can be substituted with 'chrUn' using ``vim``. Depending on the chromosome or scaffold names in a specific dataset, this step may be superfluous.

``vim -c 'g/^[^#]/s/.\{{-}}\t/chrUn\t/' -c 'wq' {input}``

The vcf file can then be converted to plink format.

``plink --vcf snps.vcf --double-id --allow-extra-chr --recode --out snps --make-bed`` 

Only the bed, bim and fam output files are required by ``fastStructure``. The structure analysis is carried out for multiple values of ``-K``, which reflects the number of populations in the dataset. The range of values tested depends on prior knowledge of the dataset. A common approach is to test values of 1-10 and to determine the best fitting population size using a likelihood approach implemented in ``fastStructure chooseK``.

```
for l in {1..10};do   
    structure.py -K $l --input=snps --output=snps_structure  
done  
chooseK.py --input=snps_structure
``` 
Although ``fastStructure`` provides rudimentary plotting tools, the R package ``pophelper`` offers a more user-friendly set of functions for creating plots from ``fastStructure`` output. A simple R script can be used to generate plots for each value of K tested, and for all values of K. 

``Rscript scripts/plotStructure.R /path/to/structure/outputs names.txt groups.txt 10 `` 

The input file path must point to the directory containg the ``meanQ`` files generated by the structure analysis. The ``names.txt`` file contains all sample names (one per row) in the same order as they occured in the vcf file converted to plink format. The ``groups.txt`` is a tab-separated popmap file with a column of sample names and a column of corresponding group names, in any order. The value of 10 in the command is the total number of values tested for K.

### Inferring a maximum likelihood phylogenetic tree <a name="Tree-building"></a>  

A phylogenetic tree is a further approach to analyse diversity and relationships within the sampled population. RAxML is a widely used tool for inferring trees using maximum likelihood (ML), offering fast and multithreaded algorithms. Alternative tools are [IQ-Tree](http://www.iqtree.org/) and [MrBayes](http://nbisweden.github.io/MrBayes/). RAxML not originally designed for analysing SNP data, therefore appropriate data preparation and careful model selection are required.

Invariant sites with only heterozygous calls in all called individuals must be filtered. A simple python script can be used to remove all SNPs without an alternative alleles call and with more than 90% heterozygous calls. 

``./scripts/filterHets.py <input.vcf> 0.9 1 > <output.vcf> ``  

Next the vcf file is converted to phylip alignment format using [vcf2phylip](https://github.com/edgardomortiz/vcf2phylip).

``./vcf2phylip.py -i <input.vcf>``  

A phylogeny can then by inferred using the ASC_GTRCAT model and an ascertainment bias correction (``--asc-corr lewis``). These options help handle the SNP input, which leads to a matrix of sites that are all variable. The number of bootstrap replicates (``-#``) is generally set to at least 100. Multithreaded RAxML can be run in rapid bootstrap mode using the command below.  

```
raxmlHPC-PTHREADS-SSE3 -f a -V -T 12 -m ASC_GTRCAT --asc-corr lewis -p 12345 -x 12345 -# 100 -s snps.phy -n mysnps -o outgroup
```
The main output of the above RAxML command is the best-scoring ML tree with support values, which is a file named RAxML_bipartitions.mysnps, where "mysnps" is the run name set with ``-n``. This tree can be visualised using ``ggtree``. An alternative tool with a graphic user interface is [FigTree](http://tree.bio.ed.ac.uk/software/figtree/). An advantage of ``ggtree`` is that it allows greater automation of the visualisation process, including colouring branches based on groups. To plot the best ML tree and colour branches by group, we can execute a simple R script, using the RAxML output file and the popmap file as input.

``Rscript scripts/ggtree.R <raxml.tre> <popmap.txt>``  

## Reference <a name="References"></a>  


