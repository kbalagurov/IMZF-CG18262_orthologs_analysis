library(data.table)
library(tidyverse)
library(stringr)
library(reutils)
library(Biostrings)
library(ape)

################################################################################
#0. read list of files (seq and ID) and rename fasta
################################################################################

getwd()
setwd("")

file_list <- list.files(path="")

Seq_list_ortho <- file_list[str_detect(file_list, "[:digit:](?=_Metazoa.fa$)")]
#Вектор со всеми ID
Ortho_table_list <- file_list[str_detect(file_list, "[:digit:](?=_Metazoa.txt$)")]



Nog_list <- file_list[!str_detect(file_list, "[:digit:](?=_Metazoa.txt$)|[:digit:](?=_Metazoa.fa$)")]
Nog_table_list <- Nog_list[str_detect(Nog_list, "\\d.txt")]
Nog_table_list  <- Nog_table_list [c(-1, -6)]

Seq_list_NOG <- Nog_list[str_detect(Nog_list, "\\d.fa")]
Seq_list_NOG <- Seq_list_NOG[c(-1, -6)]


ccc <- fread("")
cccc <- as.data.table(ccc)




general_fasta_ortho <- data.table(name_and_seq=character())
for (j in 1:length(Seq_list_ortho)) {
    #j <- 1
    str <- read_lines(Seq_list_ortho[j], skip_empty_rows = TRUE)
    str <- str_replace_all(str, '\\\\','')
    
 
    b <- fread(Ortho_table_list[j])
    b[, ID := gsub('(^\\w+\\s)', "", b[, organism_name])]
    b[, ID := str_sub(b[, ID], 1, 3)]
    b[, ID := str_c(str_sub(b[, organism_name], 1, 1),".", b[,ID])]
    b[, pub_gene_id := str_replace_all(pub_gene_id, '\\\\','')]
    tmp_a <- gsub("(_Metazoa.txt)$", "", Ortho_table_list[j])
    tmp_b <- cccc[flybase_gene_id==tmp_a,][1,1]
    b[, interest_gene := tmp_b]
    b[tab_NOGortho, gene_id := i.gene_id, on = "pub_gene_id"]
    tmp_fasta <-  data.table(name_and_seq=character())
    for (i in 1:(length(str))) {
      #i <- 1
      
      if ((i %% 2) != 0) {
        idid <- str_extract(str[i], ('"pub_gene_id":.*,'))
        tmp_idseq <- data.table(pub_gene_id = idid,  pr_sequence = str[i+1])
        tmp_idseq[, pub_gene_id := str_extract(pub_gene_id, (":.*"))]
        tmp_idseq[, pub_gene_id := str_replace_all(pub_gene_id, '"','')]
        tmp_idseq[, pub_gene_id := str_replace_all(pub_gene_id, ',','')]
        tmp_idseq[, pub_gene_id := str_replace_all(pub_gene_id, ':','')]
        tmp_idseq[, pub_gene_id := str_replace_all(pub_gene_id, '\\\\','')]
        tmp_string <- b[pub_gene_id == tmp_idseq[,pub_gene_id], ]
        tmp_idseq[, ':=' ( ID = tmp_string[1,8], gene_id = tmp_string[1,10], interest_gene = tmp_string[1,9])]
        new_name <- paste0('>',tmp_idseq[,ID],'_',tmp_idseq[,interest_gene],
                           '_',tmp_idseq[,gene_id],'_',tmp_idseq[,pub_gene_id],'\n', tmp_idseq[,pr_sequence] )
        
        
        if ((str_detect(str_extract(new_name, ('>.*\\n')), 'NA')) == TRUE ) { #b| (length(grep(sp_or_sp, new_name))==0)) {
          next} else {
            new_name <- strsplit(new_name, "\n")
            tmp_fasta <- rbind(tmp_fasta,new_name) 
          }
      } 
    }
    
    general_fasta_ortho <- rbind(general_fasta_ortho, tmp_fasta) 
    
}

#EggNOG

general_fasta_NOG <- data.table(name_and_seq=character())
for (j in 1:length(Seq_list_NOG)) {
  #j <- 5
  str <- read_lines(Seq_list_NOG[j], skip_empty_rows = TRUE)
  str <- str_replace_all(str, '\\\\','')
  

  b <- fread(Nog_table_list[j])
  setnames(b, colnames(b), c("pub_gene_id", "ID_2", "organism_name", "organism_taxid", "other"))
  #b <- tab
  b[, ID := gsub('(^\\w+\\s)', "", b[, organism_name])]
  b[, ID := str_sub(b[, ID], 1, 3)]
  b[, ID := str_c(str_sub(b[, organism_name], 1, 1),".", b[,ID])]
  b[, pub_gene_id := str_replace_all(pub_gene_id, '\\\\','')]
  tmp_a <- gsub("(.txt)$", "", Nog_table_list[j])
  tmp_b <- cccc[external_gene_name==tmp_a,][1,1]
  b[, interest_gene := tmp_b]
  b[tab_NOGortho, gene_id := i.gene_id, on = "pub_gene_id"]
  #b[tmp_tab_ortho, gene_id := i.gene_id, on = "pub_gene_id"]
  tmp_fasta <-  data.table(name_and_seq=character())
  for (i in 1:(length(str))) {
    #i <- 1
    
    if ((i %% 2) != 0) {
      idid <- str_extract(str[i],  "(?<=>)[^>]+")
      tmp_idseq <- data.table(pub_gene_id = idid,  pr_sequence = str[i+1])
      tmp_string <- b[pub_gene_id == tmp_idseq[,pub_gene_id], ]
      if (nrow(tmp_string) == 0) {
        print(paste(tmp_idseq[,pub_gene_id], "is not in fastalist"))
        next
      }
      tmp_idseq[, ':=' ( ID = tmp_string[1,6], gene_id = tmp_string[1,8], interest_gene = tmp_string[1,7])]
      new_name <- paste0('>',tmp_idseq[,ID],'_',tmp_idseq[,interest_gene],
                         '_',tmp_idseq[,gene_id],'_',tmp_idseq[,pub_gene_id],'\n', tmp_idseq[,pr_sequence] )
      
      if ((str_detect(str_extract(new_name, ('>.*\\n')), 'NA')) == TRUE ) { #b| (length(grep(sp_or_sp, new_name))==0)) {
        next} else {
          new_name <- strsplit(new_name, "\n")
          tmp_fasta <- rbind(tmp_fasta,new_name) 
        }
    } 
  }

  general_fasta_NOG <- rbind(general_fasta_NOG, tmp_fasta) 
  
}

general_fasta <- rbind(general_fasta_ortho, general_fasta_NOG) 



fwrite( general_fasta_ortho , file = "general_fasta_ortho_sp.fa", append = TRUE)
fwrite( general_fasta_NOG , file = "general_fasta_NOG_sp.fa", append = TRUE)
fwrite( general_fasta , file = "general_fasta_all_sp.fa", append = TRUE)


  
  
typeof(str)
length(str)

###############################################################################

shorts_of_sp<- fread(".txt", sep = "")
typeof(shorts_of_sp)

full_and_short <- data.table(full_name=character(),
                          short_name= character())

shorts_of_sp[, ID_sp := gsub('(^\\w+\\s)', "", shorts_of_sp[, full_name])]
shorts_of_sp[, ID_sp := str_sub(shorts_of_sp[, ID_sp], 1, 3)]
shorts_of_sp[, ID_sp := str_c(str_sub(shorts_of_sp[, full_name], 1, 1),".", shorts_of_sp[,ID_sp])]
fwrite(shorts_of_sp, file = ".txt", sep = "\t")

species_list <- fread(".txt")

sp_or_sp <- paste(species_list[, ID_sp], collapse = "|")
###############################################################################


library(Biostrings)  
library(stringr)    

setwd("")
getwd()

wide_seqs <- readAAStringSet("general_fasta_all_sp_rmdup_align.fa")  
short_seqs <- readAAStringSet("deleted_seqs_file.fa")


names(wide_seqs) <- str_replace(names(wide_seqs), "/.*", "")
names(short_seqs) <- str_replace(names(short_seqs), "/.*", "")

matching_names <- intersect(names(short_seqs), names(wide_seqs))  
cat("Found", length(matching_names), "matching sequences.\n")  

if (length(matching_names) == 0) stop("No matching sequence names found!")  

matched_wide_seqs <- wide_seqs[matching_names]  

short_seqs[matching_names] <- matched_wide_seqs  

writeXStringSet(short_seqs, "short_list_updated.fasta") 
