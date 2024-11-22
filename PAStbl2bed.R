

args<- commandArgs(TRUE)
args_v<-NA
if(length(args)>1){
  args_v <- args[c(FALSE, TRUE)]
  names(args_v) <- sub("\\-+","", args[c(TRUE, FALSE)])
  print(args_v)
}
args<-NULL

in_dir <- ifelse(is.na(args_v["i"]), "_i_", args_v["i"])
out_dir <- ifelse(is.na(args_v["o"]), "_o_", args_v["o"])


if (!require('tidyverse')) install.packages('tidyverse'); library('tidyverse')

PAStbl <- read.table(in_dir, header = TRUE, sep = "\t")[, c("Chromosome","Strand","Position")]
PAStbl$PasID <- paste(PAStbl$Chromosome,PAStbl$Strand,PAStbl$Position,sep=":")

### get genomic range for each PAS, (-100nt, PAS, +25nt)
PAStbl$from <- ifelse(PAStbl$Strand == "-", PAStbl$Position-25-1, PAStbl$Position-100)
PAStbl$to <- ifelse(PAStbl$Strand == "-", PAStbl$Position+100-1, PAStbl$Position+25)
PAStbl.bed <- PAStbl[c("Chromosome", "from", "to", "PasID")]
options(scipen = 999)
write.table(PAStbl.bed, out_dir, col.names = F, row.names = F, sep = "\t", quote = F)

