# IDAM #
## Brief Description ##
IDAM provides a framework for the identification of microbial disease-associated gene modules based on metagenomic and metatranscriptomic data. The raw sequencing data (matched metagenomic and metatranscriptomic data) or the expression matrix is required. The output is gene modules consisting of gene and sample subsets.
## Environment ##
IDAM is an integrated C package that requires a basic UNIX/Linux environment. The gcc compiler with version 4.8.5 or higher is required to be prior installed. More details can be found [here](https://gcc.gnu.org/wiki/InstallingGCC). Currently, IDAM does not support Mac or Windows systems.
## Usage ##
### 1. Dependencies ###
Dependencies of IDAM are listed below. You can click the software name to navigate to its website. Note that after installing each dependency, you should [add it to path](https://unix.stackexchange.com/questions/26047/how-to-correctly-add-a-path-to-path) to ensure a successful run.
#### 1) [HUMAnN2](https://huttenhower.sph.harvard.edu/humann2) (>=2.8.1)
HUMAnN2 is used to process the raw sequencing file and obtain the relative expression of each sample.
#### 2) [BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/) (>=2.10.0) 
BLAST is used for sequence alignment between the sequences from UniRef90 gene families and the gene sequences of uber-operons.


### 2. Installation ###
The source code of IDAM is freely available at https://github.com/OSU-BMBL/IDAM. To install IDAM, you can download the zip file manually from GitHub, or use the code below in Linux.
   	 
	cd /your/working/path/ 
	wget https://cloud.osubmi.com/downloadFiles/idam/idam_master.zip

Then, unzip the file:

	unzip idam_master.zip && rm idam_master.zip

Compile codes:

	cd ./idam_master
	make


### 3. Data preparation ###
To run the tutorial, first download the data [here](https://cloud.osubmi.com/downloadFiles/idam/) (1.0G) into the folder ./idam_master, or use wget code in the Linux system.	
	
	wget https://cloud.osubmi.com/downloadFiles/idam/meta_data.tar.gz

Then unzip the file.
	
	tar -xvzf meta_data.tar.gz && rm meta_data.tar.gz

In this way, you can obtain three subfolders including the data prepared for input, the data used for an example of HUMAnN2, and the sequence data used for an example of alignment.
#### 3.1 The example data prepared for input:  ./meta_data/input/
Two files are included in this folder, matrix.tsv and tblastx.tsv. The first is an expression matrix obtained from 734 metagenomic and metatranscriptomic datasets in this paper, which consists of  917,127 genes as rows and 734 samples as columns. The second is the alignment result between genes within the matrix mentioned above and genes in uber-operons.
#### 3.2 The example data used for HUMAnN2: ./meta_data/humann2/
Two files are included in this folder, metagenome.fastq.gz and metatranscriptome.fastq.gz. They correspond to the metagenomic and metatranscriptomic data from one sample, respectively. 
#### 3.3 The example data used for alignment: ./meta_data/sequence/
The file, uniref.fasta,  consists of all sequences corresponding to genes (UniRef90 identifiers). You can extract the sequences of genes from your matrix by NCBI (for gene names) or the UniProt database (for UniRef identifiers).

**Note 1:** You can run the code in 4.3 using the two files in ./meta_data/input/ and get the results in this paper (**This can be used for testing successful run**). Or, change these two files into your own data.   

**Note 2:** If your data is raw sequencing data, you have to run 4.1 for the expression matrix. Then, do an alignment between the sequences of genes and the sequence files (sequence_output.fa) via 4.2.  

**Note 3:** If you have a well-formatted expression matrix, you can skip 4.1 and do an alignment between the sequences of genes and the sequence files (sequence_output.fa, used as the aligned database) via 4.2.  

### 4. Running ###
#### 4.1 The expression matrix construction ####
This step is to obtain an expression matrix from raw sequencing reads. If you have matched metagenomic and metatranscriptomic data of multiple samples, you can first run HUMAnN2 for gene relative abundance of each sample.  Here, we use the data in *./meta_data/sequence/* as an example.

	humann2 --input  ./meta_data/humann2/metagenome.fastq.gz --output ./meta_data/humann2/metagenome_result
	humann2 --input  ./meta_data/humann2/metatranscriptome.fastq.gz --taxonomic-profile ./meta_data/humann2/metagenome_result/metagenome_humann2_temp/metagenome_metaphlan_bugs_list.tsv --output ./meta_data/humann2/metatranscriptome_result
	humann2_rna_dna_norm --input_dna ./meta_data/humann2/metagenome_result/metagenome_genefamilies.tsv --input_rna ./meta_data/humann2/metatranscriptome_result/metatranscriptome_genefamilies.tsv --example 

Then, you can merge the community-level gene expression of each sample based on gene family file and obtain a matrix, in which each row represents a gene and each column is a sample.

 
#### 4.2 The alignment with uber-operon sequences ####
This step is to obtain the sequence alignment result for the relationships between genes and uber-operons. 

First, the gene (UniRef identifiers) sequence can be extracted from the UniProt database via the [Retrieve/ID mapping tool](https://www.uniprot.org/uploadlists/). When you input the ID within the matrix from the previous step, you can get all sequences. Put them into one file, as shown in the file ./meta_data/sequence/uniref.fasta. Then, run the command.

	makeblastdb -in ./sequence_output.fa -dbtype nucl -out gene -parse_seqids
	tblastn -query ./meta_data/sequence/uniref.fasta -db gene -soft_masking true -outfmt 6 -out ./meta_data/input/tblastx.tsv -evalue 0.001 -num_threads 8 -max_hsps 1

You can obtain a result file named tblastx.tsv in the folder ./meta_data/input/.



#### 4.3 Module generation ####
This step is to output gene modules. The data in the folder ./meta_data/input/ can be run as follows. 

	./idam -i ./meta_data/input/matrix.tsv -n ./meta_data/input/tblastx.tsv -k 20
	

**Note:** This will output three files: 1) The file ./matrix.tsv.rules records the discretization rules; 2) The file ./matrix.tsv.chars records the discretization data; 3) The file **./matrix.tsv.modules** records the output modules.  


### 5. Other ###
#### Parameter description ####
-i : the input matrix must be one of two tab-delimited formats;

-n : the blast result file must be of -outfmt 6 format (BLAST+);

-w : the parameter in the combined function, default: 0.3;

-q : use quantile discretization for continuous data, default: 0.01;

-r : the number of ranks used to discrete data, default: 1;

-o : number of modules to report, default: 300;

-f : filtering overlapping modules, default: 0.1;

-k : minimum column width of the module,default: 5% of columns, minimum 2 columns.



## Contact ##
Any questions, problems, or bugs are welcome and should be dumped to [Qin Ma](Qin.Ma@osumc.edu).










