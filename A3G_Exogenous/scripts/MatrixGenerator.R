# Script to obtain a 192 unstranded tri-nucleotide matrix from variant positions
# Author: Carlos Martinez-Ruiz

#Load libraries-------------------------------------------------------
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(GenomicFeatures)
library(tidyr)
library(dplyr)
library(Biostrings)
library(fst)

#Functions-------------------------------------------------------
# Function to reverse strand in any string containing nucleotides
reverse_seq <- function(string){
  string_split <- unlist(strsplit(string, split = ""))
  string_split_new <- string_split
  if (any(grepl("A", string_split))){
    string_split_new[string_split == "A"] <- "T"
  }
  
  if (any(grepl("C", string_split))){
    string_split_new[string_split == "C"] <- "G"
  }
  
  if (any(grepl("G", string_split))){
    string_split_new[string_split == "G"] <- "C"
  }
  
  if (any(grepl("T", string_split))){
    string_split_new[string_split == "T"] <- "A"
  }
  new_string <- paste(string_split_new, collapse = "")
  return(new_string)
}

#Function to turn variant dataframe into tri nucleotide matrix
trinuc_mat_from_table <- function(rna_editing_sites){
  
  #Keep only columns of interest for the trinucleotide matrix, taking only changes by strandedness
  to_trinuc <- unique(rna_editing_sites[, c("sample_id", "var_pos", "muttype_nostrand", "flank_nostrand")])
  #Remove events without strand info
  to_trinuc <- to_trinuc[to_trinuc$muttype_nostrand != "*", ]
  
  #Keep only the tri-nucleotide context
  to_trinuc$trinuc <- gsub(pattern = "([TCGA]+)([tcga])([TCGA]+)",
                           replacement = "\\1\\U\\2\\3",
                           x = to_trinuc$flank_nostrand,
                           perl = TRUE)
  
  #Generate the trinucleotide matrix
  trinuc_mat <- to_trinuc %>% group_by(sample_id, trinuc, muttype_nostrand) %>% tally()
  #Ensure all patients have all contexts-----------------------------------------------------------------------------------------------------
  #Get all combinations possible of contexts
  all_nts <- c("A", "T", "C", "G")
  all_trinuc <- do.call(paste0, expand.grid(all_nts, all_nts, all_nts))  
  all_channels <- expand.grid(all_nts, all_nts)
  #Remove nonsensical channels (e.g. C>C)
  all_channels <- all_channels[all_channels$Var1 != all_channels$Var2, ]
  all_channels <- do.call(paste, c(all_channels, list(sep = ">")))
  
  #Get all contexts
  all_context <- expand.grid(all_channels, all_trinuc)
  #Keep only combinations that make sense (ref == 2nd nt in the trinucleotide)
  all_context$ref <- gsub(pattern = ">.",
                          replacement = "",
                          x = all_context$Var1)
  all_context$mid_nt <- gsub(pattern = "([ATCG])([ATCG])([ATCG])",
                             replacement = "\\2",
                             x = all_context$Var2)
  
  all_context <- all_context[all_context$ref == all_context$mid_nt, ]
  
  all_context <- all_context[, c("Var1", "Var2")]
  #Stop if not all 192 contexts are present 
  if (nrow(unique(all_context)) != 192){
    stop("Not all context are present!")
  }
  
  colnames(all_context) <- c("muttype_nostrand", "trinuc")
  #Add patient
  all_patients <- unique(trinuc_mat$sample_id)
  all_context <- data.frame(muttype_nostrand = rep(all_context$muttype_nostrand,
                                                   length(all_patients)),
                            trinuc = rep(all_context$trinuc,
                                         length(all_patients)),
                            sample_id = rep(all_patients, each = nrow(all_context)),
                            stringsAsFactors = FALSE)
  
  #Merge back with the context df
  trinuc_mat <- full_join(all_context, trinuc_mat)
  #Replace NAs by 0s
  trinuc_mat$n[is.na(trinuc_mat$n)] <- 0
  trinuc_mat$context <- paste(trinuc_mat$muttype_nostrand, trinuc_mat$trinuc, sep = ":")
  trinuc_mat$context <- gsub(pattern = "(.+)(:)([TCGA])([TCGA])([TCGA])",
                             replacement = "\\3[\\1]\\5",
                             x = trinuc_mat$context)
  
  #Sort by context within patient
  trinuc_mat <- trinuc_mat[order(trinuc_mat$sample_id, trinuc_mat$context), ]
  #Get a trinucleotide matrix
  trinuc_mat <- reshape2::dcast(trinuc_mat, sample_id ~ context, value.var = "n")
  rownames(trinuc_mat) <- trinuc_mat$sample_id
  trinuc_mat$sample_id <- NULL
  return(trinuc_mat)
}
#Load the gtf used for variant calling-----------------------------------------------------------------------------
gtf_txdb <- makeTxDbFromGFF("/Users/benfixman/Library/Mobile Documents/com~apple~CloudDocs/Chen_Lab/RNA_SBS/Transient overexpression of exogenous APOBEC3A causes C-to-U RNA editing of thousands of genes/gencode.v43.annotation.gtf", format = "gtf")

# Check and fix sequence names
gtf_file <- "/Users/benfixman/Library/Mobile Documents/com~apple~CloudDocs/Chen_Lab/RNA_SBS/Transient overexpression of exogenous APOBEC3A causes C-to-U RNA editing of thousands of genes/gencode.v43.annotation.gtf"
txdb <- makeTxDbFromGFF(gtf_file, format = "gtf")

# Fix sequence names
seqlevelsStyle(txdb) <- "UCSC" # Ensure UCSC naming style for the GTF file
seqlevelsStyle(BSgenome.Hsapiens.UCSC.hg38) <- "UCSC" # Ensure UCSC naming style for the reference genome

genes_hg38 <- genes(txdb)

genes_hg38 <- genes(gtf_txdb)

rna_editing_table <- read.csv("/Users/benfixman/Library/Mobile Documents/com~apple~CloudDocs/Chen_Lab/RNA_SBS/The double-domain cytidine deaminase APOBEC3G is a cellular site-specific RNA editing enzyme/data/variants/variants.csv")

rna_editing_pos_gr <- GRanges(Rle(rna_editing_table$chr),
                              IRanges(start = rna_editing_table$pos, end = rna_editing_table$pos))

#Get the nucleotides flanking each position, counting either from the start or from the end
flanks_start <- flank(rna_editing_pos_gr, 1, start = TRUE)
flanks_end <- flank(rna_editing_pos_gr, 1, start = FALSE)
rna_editing_table$start_seq <-  as.character(getSeq(BSgenome.Hsapiens.UCSC.hg38, flanks_start))
rna_editing_table$end_seq <-  as.character(getSeq(BSgenome.Hsapiens.UCSC.hg38, flanks_end))
rna_editing_table$ref <- gsub(pattern = ">.+",
                              replacement = "",
                              x = rna_editing_table$muttype)
#Put together into a single string
rna_editing_table$flank <- paste0(rna_editing_table$start_seq, tolower(rna_editing_table$ref),
                                  rna_editing_table$end_seq)

#Remove unnecesary columns
rna_editing_table$start_seq <- NULL
rna_editing_table$end_seq <- NULL
rna_editing_table$ref <- NULL

#Get strandedness of each variant based on gene expression
uniq_sites <- unique(rna_editing_table[, c("chr", "pos", "muttype", "flank")])
#Get the ref and alt alleles
uniq_sites <- uniq_sites %>% separate(muttype, into = c("ref", "alt"), sep = ">")

#Get GRanges object for RNA variants
editing_pos <- GRanges(Rle(uniq_sites$chr),
                       IRanges(start = as.numeric(uniq_sites$pos), width = 1),
                       Rle("*"),
                       REF = uniq_sites$ref, ALT = uniq_sites$alt)

#Get the strandedness of each gene from the same gtf used on variant calling
overlap_genes <- mergeByOverlaps(genes_hg38, editing_pos)

strand_varpos <- data.frame(chr = as.character(seqnames(overlap_genes$editing_pos)),
                            pos = as.integer(start(overlap_genes$editing_pos)),
                            ref = overlap_genes$editing_pos$REF,
                            alt = overlap_genes$editing_pos$ALT,
                            strand = as.character(strand(overlap_genes$genes_hg38)),
                            stringsAsFactors = FALSE)

#Remove duplicated positions (arising from different genes in each strand at the same position)
strand_varpos <- strand_varpos %>% group_by(chr, pos) %>% filter(n_distinct(strand) == 1)

#Ensure only unique overlaps are kept, these arise from the same position overlapping multiple genes with the same strandedness
strand_varpos <- unique(strand_varpos)

uniq_sites <- left_join(uniq_sites, strand_varpos)
uniq_sites$strand <- ifelse(is.na(uniq_sites$strand),
                            "*", uniq_sites$strand)
uniq_sites$muttype <- paste(uniq_sites$ref, uniq_sites$alt, sep = ">")

#Replace the editing type when "untranscribed"
uniq_sites$muttype_nostrand <- "*"
uniq_sites$muttype_nostrand[uniq_sites$strand == "-"] <- sapply(uniq_sites$muttype[uniq_sites$strand == "-"], reverse_seq)
uniq_sites$muttype_nostrand[uniq_sites$strand == "+"] <- uniq_sites$muttype[uniq_sites$strand == "+"]

#Same with flank
#Get flank positions to change
uniq_sites$flank_nostrand <- "*"
uniq_sites$flank_nostrand[uniq_sites$strand == "-"] <- as.character(reverseComplement(DNAStringSet(uniq_sites$flank[uniq_sites$strand == "-"])))
uniq_sites$flank_nostrand[uniq_sites$strand == "+"] <- uniq_sites$flank[uniq_sites$strand == "+"]

uniq_sites$flank_nostrand <- gsub(pattern = "([TCGA])([TCGA])(.+)",
                                  replacement = "\\1\\L\\2\\U\\3", perl = TRUE,
                                  x = uniq_sites$flank_nostrand)

rna_editing_table_strand <- rna_editing_table %>% mutate(muttype_nostrand = NULL,
                                                         flank_nostrand = NULL,
                                                         strand = NULL)

uniq_sites$var_pos <- paste(uniq_sites$chr, uniq_sites$pos , sep = ":")
rna_editing_table_strand$var_pos <- paste(rna_editing_table_strand$chr,
                                          rna_editing_table_strand$pos , sep = ":")
rna_editing_table_strand <- inner_join(rna_editing_table_strand, uniq_sites)

# Generate tri-nucleotide matrix and output
trinuc_mat <- trinuc_mat_from_table(rna_editing_table_strand)

write.table(trinuc_mat, "/Users/benfixman/Library/Mobile Documents/com~apple~CloudDocs/Chen_Lab/RNA_SBS/The double-domain cytidine deaminase APOBEC3G is a cellular site-specific RNA editing enzyme/data/matrices/matrix_192.tsv", quote = FALSE, sep = "\t")
