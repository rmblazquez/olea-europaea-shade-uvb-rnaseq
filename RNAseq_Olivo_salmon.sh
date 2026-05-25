## Script para análisis RNA-seq en rmblazquez@ubuntupc 
## Las lecturas están en el directorio /home/rmblazquez/Documentos/Data/reads/HN00199498/RawFASTQ/
## Además de los archivos FASTQ (son PE, tendrán extensión .1.fastq.gz y .2.fastq.gz por muestra), hay un CSV con los metadatos del experimento
## Los análisis de este script están realizados desde el home de rmblazquez

## Trabajaremos en entornos de conda para evitar contaminar de programas el servidor, así que hay que instalar miniconda lo primero
## Cada vez de que abre la sesión, hay que refrescar los bins de conda, así que primero se ejecuta el comando

source ~/miniconda3/bin/activate

## Configurar sesión de conda 

conda create -n 4rnaseq -c bioconda -c conda-forge -c defaults -c r fastqc multiqc fastp salmon

conda activate 4rnaseq

## Generar árbol de directorios

projectDir="/home/rmblazquez/Documentos/Resultados/RNAseq_acebuche"

mkdir -p $projectDir/FastQC_raw $projectDir/Fastp $projectDir/FastQC_trimm $projectDir/Salmon

## Generar lista de archivos FASTQ para el análisis

readsPath="/home/rmblazquez/Documentos/Data/reads/HN00199498/RawFASTQ"

#cut -f2 $readsPath/RNAseq_olivo_design.csv | tail -n +2 > $projectDir/lista.txt
for i in {1..49}; do echo $i; done > $projectDir/lista.txt

## Comprobar calidad de lecturas con FastQC sin filtrar

cat $projectDir/lista.txt | while read line; do
  mkdir $projectDir/FastQC_raw/$line.1.fastq
  fastqc -o $projectDir/FastQC_raw/$line.1.fastq $readsPath/$line.1.fastq.gz 
  mkdir $projectDir/FastQC_raw/$line.2.fastq
  fastqc -o $projectDir/FastQC_raw/$line.2.fastq $readsPath/$line.2.fastq.gz
done

cd $projectDir/FastQC_raw

multiqc .

## Eliminar adaptadores y secuencias y posiciones de baja calidad

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

## Comprobar calidad de lecturas con FastQC filtradas

cat $projectDir/lista.txt | while read line; do
  mkdir $projectDir/FastQC_trimm/$line.1.trimm.fastq
  fastqc -o $projectDir/FastQC_trimm/$line.1.trimm.fastq $projectDir/Fastp/$line.1.trimm.fastq.gz 
  mkdir $projectDir/FastQC_trimm/$line.2.trimm.fastq
  fastqc -o $projectDir/FastQC_trimm/$line.2.trimm.fastq $projectDir/Fastp/$line.2.trimm.fastq.gz
done

cd $projectDir/FastQC_trimm

multiqc .

## Indexar transcriptoma (CDS) de olivo var. picual

cd $projectDir/Salmon

wget "https://genomaolivar.dipujaen.es/downloads/Sequences/Olea_europaea_cDNA_v061.fasta.gz" -O Olea_europaea_cDNA_v061.fasta.gz

salmon index -t $projectDir/Salmon/Olea_europaea_cDNA_v061.fasta.gz -i $projectDir/Salmon/Oleur_salmon_cdna_index -k 31

## Cuantificar transcriptoma de las muestras

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

## Generar tabla de cuantificación (matriz de expresión)

cat $projectDir/lista.txt | while read line; do
  cut -f1 $projectDir/Salmon/$line.quant/quant.sf > $projectDir/$line.referencia.sf
  cut -f5 $projectDir/Salmon/$line.quant/quant.sf > $projectDir/Salmon/$line.quant.sf
  paste $projectDir/Salmon/$line.referencia.sf $projectDir/Salmon/$line.quant.sf > $projectDir/Salmon/$line.quant.csv
  rm $projectDir/Salmon/$line.referencia.sf $projectDir/Salmon/$line.quant.sf
done

## Hay dos opciones a partir de la cuantificación por salmon (elige una):

# I) Generar "manualmente" la matriz de expresión con UNIX 

# tr '\n' '\t' < lista.txt > rnaseqHeader
# paste *.csv | cut -f1,2,4,6,8,10 | tail -n +2 | cat rnaseqHeader - > OleurQuantSalmon.csv
# rm rnaseqHeader *_quant.csv

# II) Descargarse los directorios de 'salmon' y usar la librería de R 'tximport' para generarla en R
