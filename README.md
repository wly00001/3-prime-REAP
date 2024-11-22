
# 3'REAP (3’ Reads Enrichment using Annotated PolyA sites) pipeline, using reverse reads data of any type of 3' RNA-Seq methods as input. 

#### By Luyang Wang, lwang@wistar.org / wly00001@gmail.com, Bin Tian Lab @ The Wistar Institute


-Before running the code below, please create a data name list, DataNameList.txt, which contains only one column of each raw sequencing data name. For example, if the raw sequencing data name is 12345-02-01-01_S168_L008_R1_001.fastq, the name should be put into the txt file should be "12345-02-01-01_S168_L008_R1_001".

-Please also have "TrimGalore", "cutadapt", "STAR", "bedtools", "genomeCoverageBed", UCSC "norm_bedgraph.pl", UCSC "bedGraphToBigWig", UCSC "chrom.sizes" file, and "R" ready.

### Set your WORK DIR
```
WORK_DIR=/full/path/of/your/project
cd $WORK_DIR
```


### generate PAS STAR index
Convert PAS_table into PAS_table.bed
PAS_table_full_path is the full path of a tab-separated PAS table, in which rows are PASs and columns must contain "Chromosome","Strand","Position" for each PAS.
```
Rscript PAStbl2bed.R -i full_path_of_PAS_table -o full_path_of_PAS_table.bed
```
Convert PAS_table.bed into PAS_table.fasta
```
bedtools getfasta -fi Your_genome.fasta -bed full_path_of_PAS_table.bed -name > full_path_of_PAS_table.fasta
```
Generate PAS STAR index
```
STAR --runThreadN 2 --runMode genomeGenerate --genomeDir full_path_of_PAS_table_STAR_index/ --genomeFastaFiles full_path_of_PAS_table.fasta
```



### Use TrimGalore to trim adapters
```
mkdir TrimGalore
while IFS=$' \t\r\n' read -r sample; do
trim_galore /full/path/of/your/raw/data/${sample}.fastq -o $WORK_DIR/TrimGalore/
done < DataNameList.txt
```



### Some reverse reads contain 5'Ts, remove remaining 5'Ts before the alignment
```
python trim_5T.py --rawfastq_dir $WORK_DIR/TrimGalore --project_dir $WORK_DIR
```



### Mapping
```
mkdir star_out
while IFS=$' \t\r\n' read -r sample; do
STAR --runThreadN 2 --genomeDir full_path_of_PAS_table_STAR_index/ --readFilesIn ./TrimGalore/${sample}.5Ttrimmed.fastq --outFilterType BySJout --outFilterMultimapNmax 1 --outFilterScoreMinOverLread 0.2 --outFilterMatchNminOverLread 0.2 --outSAMattributes NH HI NM MD --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $WORK_DIR/star_out/${sample}
done < DataNameList.txt
```



### Define LAP (last aligned position)
```
mkdir $WORK_DIR/LAP
cd $WORK_DIR/LAP
while IFS=$' \t\r\n' read -r sample; do
bedtools bamtobed -cigar -i $WORK_DIR/star_out/${sample}Aligned.sortedByCoord.out.bam > ${sample}.bed
sort -k 1,1 ${sample}.bed > ${sample}.bed.sorted
done < DataNameList.txt
wc -l *.bed.sorted
cd $WORK_DIR
```



### Match LAP with PAS, and generate PAS count table, only reads whose LAP is within +/-24nt PAS will be kept.
```
mkdir result
mkdir result/csv
while IFS=$' \t\r\n' read -r sample; do
Rscript PAS_match.R -bedLAP $WORK_DIR/LAP/${sample}.bed.sorted -out $WORK_DIR/result/csv/${sample}
done < DataNameList.txt
```



### Generate bigwig for LAP
```
chromsizes=your.chrom.sizes # Chromosome sizes dir/file
mkdir $WORK_DIR/bigwig_LAP
cd $WORK_DIR/bigwig_LAP
while IFS=$' \t\r\n' read -r sample; do
  echo "Working on $sample..."
  ## Count total read number
  totalReadNum=`wc -l $WORK_DIR/result/csv/${sample}_PASS_bw.bed | sed s/[[:blank:]].*//`
  echo "for file ${sample}_PASS_bw.bed, TotalReadNum=$totalReadNum"

  sort -k 1,1 $WORK_DIR/result/csv/${sample}_PASS_bw.bed > ./${sample}_PASS_bw.bed.sorted

  ## it is strand-specific 
  ## if reverse data: + on genomeCoverageBed is "minus", and - is "plus"
  ## if forward data: + on genomeCoverageBed is "plus", and - is "minus"
  ## Generate bedgraph file
  echo "Generate bedgraph files for + and - strands..."
  genomeCoverageBed -bg -split -i ./${sample}_PASS_bw.bed.sorted -strand '-' -g $chromsizes > $sample.plus.bedgraph
  genomeCoverageBed -bg -split -i ./${sample}_PASS_bw.bed.sorted -strand '+' -g $chromsizes > $sample.minus.bedgraph

  ## Normalize bedgraph counts
  echo "Normalize bedgraph counts..."
  norm_bedgraph.pl -t $totalReadNum -i "$sample.plus.bedgraph"
  norm_bedgraph.pl -t $totalReadNum -i "$sample.minus.bedgraph"

	## give minus strand negative value
  awk -v FS="\t" -v OFS="\t" 'NR>1 {print $1, $2, $3, -$4}' $sample.minus.bedgraph.normolized > $sample.minus.bedgraph.normolized1 

  ## Convert to bigwig file
  echo "Convert to bigwig file..."
  bedGraphToBigWig $sample.plus.bedgraph.normolized  $chromsizes $sample.plus.bw
  bedGraphToBigWig $sample.minus.bedgraph.normolized1  $chromsizes $sample.minus.bw

done < DataNameList.txt

chmod 775 *.bw #change file permission
rm *.bedgraph
rm *.bedgraph.normolized
rm *.bedgraph.normolized1
cd $WORK_DIR
```



### Combine all PAS count tables into one csv table.
```
Rscript combine_all_sample_PAS_count_tables.R -csv ./result/csv -out ./result/cluster.all.reads.csv
```



### The cluster.all.reads.csv is the PAS count table of all samples, in which rows are PASs and columns are samples.



