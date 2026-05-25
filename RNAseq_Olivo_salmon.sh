## Script for RNA-seq analysis in rmblazquez@ORDENADORANALISIS
## Make sure to include the metadata in the analysis directory
## Besides, install miniconda first, since the script will need a conda environment 
## To ensure conda bins are in the path, execute this command

source ~/miniconda3/bin/activate

## Create a conda environment

conda create -n 4rnaseq -c bioconda -c conda-forge -c defaults -c r fastqc multiqc fastp salmon

conda activate 4rnaseq

## Create directories for the analysis (introduce a custom path)

projectDir="/home/rmblazquez/Documentos/Resultados/RNAseq_acebuche"

mkdir -p $projectDir/FastQC_raw $projectDir/Fastp $projectDir/FastQC_trimm $projectDir/Salmon

## Generate FASTQ names list (since samples are named after 1 to 49)

readsPath="/home/rmblazquez/Documentos/Data/reads/HN00199498/RawFASTQ"

for i in {1..49}; do echo $i; done > $projectDir/lista.txt

## Check raw sequences quality with FastQC

cat $projectDir/lista.txt | while read line; do
  mkdir $projectDir/FastQC_raw/$line.1.fastq
  fastqc -o $projectDir/FastQC_raw/$line.1.fastq $readsPath/$line.1.fastq.gz 
  mkdir $projectDir/FastQC_raw/$line.2.fastq
  fastqc -o $projectDir/FastQC_raw/$line.2.fastq $readsPath/$line.2.fastq.gz
done

cd $projectDir/FastQC_raw

multiqc .

## Remove adapters and low quality positions and sequences with Fastp 

cd $projectDir

cat $projectDir/lista.txt | while read line; do
  fastp --in1 $readsPath/$line.1.fastq.gz \
  --in2 $readsPath/$line.2.fastq.gz \
  --out1 $projectDir/Fastp/$line.1.trimm.fastq.gz \
  --out2 $projectDir/Fastp/$line.2.trimm.fastq.gz \
  --thread 12
  --detect_adapter_for_pe \
  --cut_front \
  --cut_tail \
  --cut_window_size 12 \
  --cut_mean_quality 30 \
  --length_required 75 \
  --json $projectDir/Fastp/$line.json \
  --html $projectDir/Fastp/$line.html
done

## Check sequence quality again, this time from the trimmed sequences

cat $projectDir/lista.txt | while read line; do
  mkdir $projectDir/FastQC_trimm/$line.1.trimm.fastq
  fastqc -o $projectDir/FastQC_trimm/$line.1.trimm.fastq $projectDir/Fastp/$line.1.trimm.fastq.gz 
  mkdir $projectDir/FastQC_trimm/$line.2.trimm.fastq
  fastqc -o $projectDir/FastQC_trimm/$line.2.trimm.fastq $projectDir/Fastp/$line.2.trimm.fastq.gz
done

cd $projectDir/FastQC_trimm

multiqc .

## Download and index the Olea europaea var. Picual transcriptome (cDNA) with Salmon

cd $projectDir/Salmon

wget "https://genomaolivar.dipujaen.es/downloads/Sequences/Olea_europaea_cDNA_v061.fasta.gz" -O Olea_europaea_cDNA_v061.fasta.gz

salmon index -t $projectDir/Salmon/Olea_europaea_cDNA_v061.fasta.gz -i $projectDir/Salmon/Oleur_salmon_cdna_index -k 31

## Quantification of gene expression per sample 

cd $projectDir

cat $projectDir/lista.txt | while read line; do
  salmon quant -i $projectDir/Salmon/Oleur_salmon_cdna_index \
  -l IU \
  -1 $projectDir/Fastp/$line.1.trimm.fastq.gz \
  -2 $projectDir/Fastp/$line.2.trimm.fastq.gz \
  -p 12 \
  --validateMappings \
  -o $projectDir/Salmon/$line.quant
done

## Generate expression matrix (optional)

cat $projectDir/lista.txt | while read line; do
  cut -f1 $projectDir/Salmon/$line.quant/quant.sf > $projectDir/$line.referencia.sf
  cut -f5 $projectDir/Salmon/$line.quant/quant.sf > $projectDir/Salmon/$line.quant.sf
  paste $projectDir/Salmon/$line.referencia.sf $projectDir/Salmon/$line.quant.sf > $projectDir/Salmon/$line.quant.csv
  rm $projectDir/Salmon/$line.referencia.sf $projectDir/Salmon/$line.quant.sf
done

## Salmon output is to be analyzed with R library 'tximport'
