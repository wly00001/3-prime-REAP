
if (!require('tidyverse')) install.packages('tidyverse'); library('tidyverse')
if (!require('GenomicFeatures')) install.packages('GenomicFeatures'); library('GenomicFeatures')
if (!require('GenomicRanges')) install.packages('GenomicRanges'); library('GenomicRanges')
if (!require('GenomicAlignments')) install.packages('GenomicAlignments'); library('GenomicAlignments')
if (!require('plyr')) install.packages('plyr'); library('plyr')



args<- commandArgs(TRUE)
args_v<-NA
if(length(args)>1){
   args_v <- args[c(FALSE, TRUE)]
   names(args_v) <- sub("\\-+","", args[c(TRUE, FALSE)])
   print(args_v)
}
args<-NULL

bed_dir <- ifelse(is.na(args_v["bedLAP"]), "_bedLAP_", args_v["bedLAP"])
out_dir <- ifelse(is.na(args_v["out"]), "_out_", args_v["out"])

sample <- gsub(".bed.sorted","",tail(strsplit(bed_dir,"/")[[1]],1))
print(paste0("sample name, ", sample))

bed <- read.table(bed_dir, header = F, sep = "\t", quote = "", stringsAsFactors = FALSE)
names(bed) <- c("chr", "start", "end", "readID", "MAPQ", "strand", "CIGAR")

bed$from <- ""
bed$to <- ""
bed$from[bed$strand == "-"] <- bed$start[bed$strand == "-"]
bed$to[bed$strand == "-"] <- bed$end[bed$strand == "-"]
bed$from[bed$strand == "+"] <- 125-bed$end[bed$strand == "+"]
bed$to[bed$strand == "+"] <- 125-bed$start[bed$strand == "+"]
bed$from <- as.numeric(bed$from)
bed$to <- as.numeric(bed$to)
print(paste0("mappable read #, ", length(unique(bed$readID))))
bed$strand_rev <- ""
bed$strand_rev[bed$strand == "-"] <- "+"
bed$strand_rev[bed$strand == "+"] <- "-"

df <- bed
df <- df %>% separate(chr, sep = "::", into = c("PasID","miniREFid"), remove = T)
df <- df %>% separate(PasID, sep = ":", into = c("Chromosome","Strand","Position"), remove = F)

df2 <- subset(df, Strand==strand_rev)
df2$d_LAP2PAS <- abs(df2$to-100)
df3_0 <- df2 %>% 
  group_by(readID) %>% 
  dplyr::slice(which.min(d_LAP2PAS))   ### deal with multimappers, useless for unique reads
df3_0$CIGAR_head_symbol <- sapply(explodeCigarOps(df3_0$CIGAR), head, 1)
df3_0$CIGAR_tail_symbol <- sapply(explodeCigarOps(df3_0$CIGAR), tail, 1)
df3_0$CIGAR_head_value <- sapply(explodeCigarOpLengths(df3_0$CIGAR), head, 1)
df3_0$CIGAR_tail_value <- sapply(explodeCigarOpLengths(df3_0$CIGAR), tail, 1)
df3_0$read_from <-""
df3_0$read_from[df3_0$Strand == "+" & df3_0$CIGAR_tail_symbol=="S"] <- df3_0$CIGAR_tail_value[df3_0$Strand == "+" & df3_0$CIGAR_tail_symbol=="S"]+1
df3_0$read_from[df3_0$Strand == "+" & df3_0$CIGAR_tail_symbol!="S"] <- 1
df3_0$read_from[df3_0$Strand == "-" & df3_0$CIGAR_head_symbol=="S"] <- df3_0$CIGAR_head_value[df3_0$Strand == "-" & df3_0$CIGAR_head_symbol=="S"]+1
df3_0$read_from[df3_0$Strand == "-" & df3_0$CIGAR_head_symbol!="S"] <- 1
df3_0$read_from <- as.numeric(df3_0$read_from)
df3 <- subset(df3_0, d_LAP2PAS<=24)
print(paste0("PAS matched (PASS) # (<=24), ", nrow(df3)))
df3$Position <- as.numeric(df3$Position)


### generate bed file which will be used for creating UCSC genome browser bigwig file.
PASS_bw_bed <- df3[c("Chromosome", "Position", "start", "end", "readID", "MAPQ", "strand")]
PASS_bw_bed$start_GENOME <- ""
PASS_bw_bed$end_GENOME <- ""
PASS_bw_bed$start_GENOME <- as.numeric(PASS_bw_bed$start_GENOME)
PASS_bw_bed$end_GENOME <- as.numeric(PASS_bw_bed$end_GENOME)
PASS_bw_bed$end_GENOME[PASS_bw_bed$strand == "-"] <- PASS_bw_bed$end[PASS_bw_bed$strand == "-"]-100+PASS_bw_bed$Position[PASS_bw_bed$strand == "-"]
PASS_bw_bed$start_GENOME[PASS_bw_bed$strand == "-"] <- PASS_bw_bed$end_GENOME[PASS_bw_bed$strand == "-"]-1
PASS_bw_bed$start_GENOME[PASS_bw_bed$strand == "+"] <- PASS_bw_bed$start[PASS_bw_bed$strand == "+"]-25+PASS_bw_bed$Position[PASS_bw_bed$strand == "+"]
PASS_bw_bed$end_GENOME[PASS_bw_bed$strand == "+"] <- PASS_bw_bed$start_GENOME[PASS_bw_bed$strand == "+"]+1
PASS_bw_bed <- PASS_bw_bed[c("Chromosome", "start_GENOME", "end_GENOME", "readID", "MAPQ", "strand")]
options(scipen = 999)
write.table(PASS_bw_bed, paste0(out_dir,"_PASS_bw.bed"), col.names = F, row.names = F, sep = "\t", quote = F)

cluster.all.reads <- data.frame(table(df3$PasID))
print(paste0("PASS# in csv, ", sum(cluster.all.reads$Freq)))
names(cluster.all.reads) <- c("PasID",sample)
write.csv(cluster.all.reads, paste0(out_dir,"_cluster.all.reads.csv"), col.names = T, row.names = F, quote = F)





