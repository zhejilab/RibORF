RibORF: Identifying genome-wide translated open reading frames using ribosome profiling. 

Contact: Zhe Ji (zhe.ji@northwestern.edu)

Here I present the RibORF software originally published in (Ji et al., eLife, 2015). I made several significant changes of the RibORF algorithm to make the software more powerful and user-friendly. First, I used the logistic regression model to calculate the translated probability, instead of support vector machine (SVM), because logistic Regression function was included by default R installation, while SVM was not. SVM seems to perform slightly better to distinguish non-ribosomal protein-RNA complex, but the differences are pretty minor. Second, this version of RibORF can train the prediction parameters based on the read distribution pattern of each individual dataset, using canonical ORFs in mRNA as the positive training set, and internal off-frame ORFs as the negative training set. As the 3-nt periodicity patterns of different ribosome profiling datasets are quite variable, the modification can make the prediction more accurate. Users will obtain a statistical summary of the algorithm performance, such as false discovery rates and negative discovery rates using different cutoffs. The prediction accuracy is correlated with ribosome profiling data quality. Third, I introduced a new learning parameter to identify translated ORFs, i.e. fraction of codons with 1st nucleotides containing the maximum number of reads. This parameter helps to reduce false positives resulting from enriched reads in a small subset of codons.

Materials:
Ribosome profiling datasets in Fastq format;
Genome assembly file in Fasta format;
Ribosomal RNA (rRNA) sequence file in Fasta format;
Transcript annotation file in genePred format;
Linux high performance computing cluster;
Perl program installation;
R program installation in the PATH;
Read mapping software (such as Bowtie(Langmead and Salzberg, 2012) and Tophat(Kim et al., 2013))

Protocol steps 

1. Download RibORF package from https://github.com/zhejilab/RibORF/. 

2. Obtain genome annotation files, including the genome assembly file in Fasta format and reference transcriptome annotation file in genePred format. Run “ORFannotate.pl” to get candidate ORFs in transcripts. The program allows users to pick start codons types and select ORF length cutoff. The default setting considers 5 types of start codons “ATG/CTG/GTG/TTG/ACG”, which are most frequently used ones. The program also annotate candidate ORF types based on the transcript types and ORF locations in the transcripts, such as canonical ORFs, uORFs in 5’UTRs and internal off-frame ORFs in coding regions. 

Usage: perl ORFcandidate.pl -g genomeSequenceFile -t transcriptomeFile -o outputDir  [-s startCodon] [-l orfLengthCutoff]
   -g genomeSequenceFile: the genome assembly file in fasta format;
   -t transcriptomeFile: the reference transcriptome annotation file in genePred format;
   -o outputDir: output directory;
   -s startCodon [optional]: start codon types to be considered separated by “/”, default: ATG/CTG/GTG/TTG/ACG;
   -l orfLengthCutoff [optional]: cutoff of minimum candidate ORF length, default: 6nt.

Example commands: 
2a. Download human genome and transcriptome annotation files from GENCODE(Harrow et al., 2012). 
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/gencode.v28.annotation.gtf.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/GRCh38.primary_assembly.genome.fa.gz
2b. Get gtfToGenePred command line from UCSC Genome Browser, and use the tool to convert GTF file to genePred format. 
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred
gtfToGenePred gencode.v28.annotation.gtf gencode.v28.annotation.genePred.txt
2c. Run “ORFannotate.pl” and generate candidate ORFs in genePred format. 
perl ORFannotate.pl -g GRCh38.primary_assembly.genome.fa -t gencode.v28.annotation.genePred.txt -o outputDir

There will be 2 files generated in the output directory, including “candidateORF.genepred.txt” with candidate ORFs in genePred format, and “candidateORF.fa” with candidate ORF sequences in Fasta format. We generated the candidate ORF IDs with the following format: “TranscriptID:chromatin:strand|RankNumber|transcriptLength:startCodonPosition: stopCodonPosition|candidateORFType|startCodonType” (as an example: ENST00000420190.6:chr1:+|1|1578:87:357|uORF|TTG). 

For the genePred format, each row contains the following information of a transcript (take “ENST00000377898.3” as an example): 
“Name of transcript”: ENST00000377898.3
“Chromosome”: chr1
“Strand”: +
“Transcription start site”: 6244191
“Transcription end site”: 6245578
“Start codon position”: 6244365
“Stop codon position”: 6245507
“Number of exons”: 4
“Exon start positions”: 6244191,6244350,6244547,6245109,
“Exon end positions”: 6244241,6244446,6244629,6245578,

3. Obtain the ribosome profiling dataset with cycloheximide treatment or without drug treatment, and run “removeAdapter.pl” to trim 3’ adapters of sequencing reads. For some datasets, 3’ adapters were sequenced. This helps to precisely define the length of inserted RNA fragments. 

Usage: perl removeAdapter.pl -f fastqFile -a adapterSequence -o outputFile [-l readLengthCutoff]
   -f fastqFile: raw sequencing reads in fastq format;
   -a adapterSequence: sequence of 3’ adapters, first 10nt is recommended;
   -o outputFile: output file;
   -l readLengthCutoff [optional]: minimal read length after trimming 3’ adapters, default is 15nt.

Example commands: 
3a. Download an example ribosome profiling dataset from GEO database using the fastq-dump command line (available from the NIH software sratoolkit, https://www.ncbi.nlm.nih.gov/sra/docs/toolkitsoft/).  
fastq-dump -Z SRR1802146 > SRR1802146.fastq 
3b. Remove 3’ adapters of sequencing reads.
perl removeAdapter.pl -f SRR1802146.fastq -a CTGTAGGCAC -o adapter.SRR1802146.fastq

4. Map trimmed reads to rRNAs, and then map non-rRNA reads to the reference transcriptome and genome.

Example commands:
4a. Get human rRNA sequences from NCBI database, including 5S rRNA (NR_023363), 5.8S rRNA (NR_003285), 18S rRNA (NR_003286) and 28S rRNA (NR_003287). Put the rRNA sequences in the file “human.ribosomal.rna.fa” with fastq format.
4b. Use Bowtie to index rRNA sequences. 
bowtie2-build human.ribosomal.rna.fa hg.ribosome
4c. Align trimmed ribosome profiling reads from step 3b to rRNAs, and obtain non-rRNA reads.
bowtie2 -x hg.ribosome -U adapter.SRR1802146.fastq --un norrna.adapter.SRR1802146.fastq -S ribosome.adapter.SRR1802146.fastq
4d. Align non-rRNA reads to the reference transcriptome and genome, and obtain the alignment file in SAM format.
tophat --GTF gencode.v28.annotation.gtf --no-convert-bam -o outputDir GRCh38genome.index  norrna.adapter.SRR1802146.fastq

5. Run “readDist.pl” to group reads based on fragment length, and check the distance between 5' ends of the reads around the start and stop codons of canonical ORFs of mRNAs. The length of ribosome protected fragments is ~30nt. Users can specify the fragment length (i.e. 28,29,30), and examine read distribution. 

During ribosome profiling, RNase I cannot always completely digest RNA regions unprotected by protein complexes. As a result, some reads do not show clear 3-nt periodicity across ORFs, and these reads cannot be used to identify genome-wide translated ORFs. It is important to examine the read distribution across canonical ORFs of mRNAs, and ensure that reads show clear 3-nt periodicity. 

Usage: perl readDist.pl -f readFile -g geneFile -o outputDir [-d readLength] [-l leftNum] [-r rightNum]
   -f readFile: read alignments to the reference transcriptome and genome in SAM format;
   -g geneFile: canonical protein-coding ORF annotation file in genePred format;
   -o outputDir: output directory;
   -d readLength [optional]: specified RPF length (nt), default: 25,26,27,28,29,30,31,32,33,34,;
   -l leftNum [optional]: N nucleotides upstream start codon and downstream stop codon, default: 30;
   -r rightNum [optional]: N nucleotides downstream start codon and upstream stop codon, default: 50.

Example command: 
perl readDist.pl -f SRR1802146.mapping.sam -g gencode.v28.annotation.genePred.txt -o outputDir -d 28,29,30, -l 40 -r 70

Several files will be generated in the output directory. “plot.readDist.*.pdf” shows the plot of read distribution around the start and stop codons of canonical ORFs. “read.dist.sample.*.txt” shows numeric read densities around codons, indicated by read per million (RPM) values. “sta.read.dist.*.txt” contains read number statistics, and fractions of reads in each nucleotide of codons (1st, 2nd and 3rd). High quality reads show clear 3-nt periodicity, with high percentage of reads in 1st nucleotides of codons (>50% is recommended). Low quality reads do not show obvious 3-nt periodicity, and cannot be used for further analyses. 

6. Run “offsetCorrect.pl” and correct read locations based on offset distances between 5’ ends and ribosomal A-sites. The offset distances can be inferred from the plots from step 5. Based on the read distribution around the start and stop codon of canonical ORFs, users can manually check the offset distances between 5’ ends of the reads and ribosomal A-site. Users can manually check the read distribution plots and ensure the data quality. Put correction parameters in a file, i.g. "offset.corretion.parameters.txt", with 2 columns. The first column shows the read fragment length, and the second column shows the offset distance. Ribosomal profiling experiments can have different offset correction parameters. 

Usage: perl offsetCorrect.pl -r readFile -p offsetParameterFile -o readCorrectedFile
   -r readFile: read mapping file before offset correction in SAM format;
   -p offsetParameterFile: parameters for offset correction, 1st column: read length, 2nd column: offset distance;
   -o readCorrectedFile: output file after offset correction in SAM format. 

Example commands: 
6a. Generate the file “offset.corretion.parameters.txt”, with the content as following.
28	15
29	16
30	16
6b. Run “offsetCorrect.pl”. 
perl offsetCorrect.pl -r SRR1802146.mapping.sam -p offset.corretion.parameters.txt -o corrected. SRR1802146.mapping.sam

7. Run “readDist.pl” and check the corrected read locations around the start and stop codons of canonical ORFs. This step is to check whether read locations after offset correction show clear in-frame 3-nt periodicity. As the read length after offset correction is 1, put the parameter “-d” as “1”. 

Usage: perl readDist.pl -f readFile -g geneFile -o outputDir -d 1 [-l leftNum] [-r rightNum]
   -f readFile: read alignments to the reference transcriptome and genome, in SAM format;
   -g geneFile: canonical protein-coding ORF annotation file, in genePred format;
   -o outputDir: output directory;
   -d readLength:1;
   -l leftNum [optional]: N nucleotides upstream start codon and downstream stop codon, default: 30;
   -r rightNum [optional]: N nucleotides downstream start codon and upstream stop codon, default: 50.

Example command:
perl readDist.pl -f corrected.SRR1802146.mapping.sam -g gencode.v28.annotation.genePred.txt -o outputDir -d 1

8. Run “ribORF.pl” to identify translated ORFs using corrected ribosome profiling read alignment file from step 6 and the candidate ORF file from step 2. Users can pick cutoff parameters, including ORF length, supported read number and predicted translated P-value. 

Usage: perl ribORF.pl -f readCorrectedFile -c candidateORFFile -o outputDir [-l orfLengthCutoff] [-r orfReadCutoff] [-p predictPvalueCutoff]
   -f readCorrectedFile: input read mapping file after offset correction in SAM format;
   -c candidateORFFile: candidate ORFs in genePred format;
   -o outputDir: output directory, with files reporting learning parameters and predicted translating probability;
   -l orfLengthCutoff [optional]: cutoff of ORF length (nt), default: 6;
   -r orfReadCutoff [optional]: cutoff of supported read number, default: 11.
   -p predictPvalueCutoff [optional]: cutoff used to select predicted translated ORF, default: 0.7.

Example: perl ribORF.pl -f corrected.SRR1802146.mapping.sam -c candidateORF.genepred.txt -o outputDir

A few output files will be generated. “pred.pvalue.parameters.txt” contains training parameters for candidate ORFs and predicted P-values. The columns include the following: candidate ORF ID (orfID), chromosome (chrom), strand, start codon location (codon5), stop codon location (codon3), ORF length, supporting read number (readNum), fraction of reads in 1st nucleotides of codons (f1), fraction of reads in 2nd nucleotides of codons (f2), fraction of reads in 3rd nucleotides of codons (f3), entropy value of read distribution (entropy), maximum entropy value of randomized distribution (MAXentropy), percentage of maximum entropy value (PME), number of codons with sequencing reads (codonNum), fraction of codons with 1st nucleotides containing more reads than 2nd and 3rd (f1max), and predicted translated probability (pred.pvalue). The predicted translated P-value for each candidate ORFs was calculated using the logistic regression, and was based on 3 input parameters, including f1, PME and f1max. 

The predicted translated P-values represent the binary classification of the candidate ORFs. The values can be smaller than 0 or be great than 1. P-values closer to 1 represents that the ORF is likely to be translated. Users can select the cutoff P-value to define translated ORFs. “stat.cutoff.txt” file contains the cutoffs and associated statistics of true positive, false positive, true negative and false negative estimations. Users can take the numbers as the reference to pick appropriate cutoff. “plot.ROC.curve.pdf” contains the ROC curve plot based on “false positive rates” and “true positive rates” from the “stat.cutoff.txt” file. Users can use the following R script to estimate Area Under ROC Curve (AUC) value. Users can run the following R script to get the ROC curve and the AUC value, and the R script requires “MESS” package installation.

A <- read.table ("stat.cutoff.txt", sep="\t", header=T)
fpr <- A[,6]
tpr <- A[,7]
plot(fpr, tpr, col=0)
lines(fpr, tpr,col=1, lwd=3)
library("MESS")
auc(fpr,tpr, type = 'spline')

Many candidate ORFs can overlap with each other with the same stop codon and different start codons. In this case, a representative ORF is selected based on the following criteria: we first pick AUG as start codons if present, and we then choose 5’ most start codon as the representative one. But if there is no read between the picked one and the next downstream candidate, we choose the next one as the representative start codon. The output files “repre.valid.pred.pvalue.parameters.txt” and “repre.valid.ORF.genepred.txt” contain the information of representative ORFs. 

