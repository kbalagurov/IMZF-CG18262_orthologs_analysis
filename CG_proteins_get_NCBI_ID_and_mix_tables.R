library(data.table)
library(tidyverse)
library(stringr)
library(reutils)
library(taxize)
library(ape)
library(ggtree)
library(Biostrings)


getwd()
setwd("")

getwd()


####################################################################################
#EggNOG


file_list <- list.files(path="")
#EGGNOG - GeneName (CG1602)

Nog_list <- file_list[!str_detect(file_list, "[:digit:](?=_Metazoa.txt$)|[:digit:](?=_Metazoa.fa$)")]

Nog_table_list <- Nog_list[str_detect(Nog_list, "\\d.txt")]
Nog_table_list <- Nog_table_list[c(-1, -6)]

tmp_tab_nog <- data.table(ID_1=character(),
                          ID_2= character(),
                          organism_name=character(),
                          organism_taxid= character(),
                          other= character(),
                          ID_sp = character(), 
                          interest_gene = character(),
                          XP_NCBI = character(),
                          gene_id = character())


for (i in 1:length(Nog_table_list)) {
  #i <- 1
  #str <- read_lines(Seq_list[i])
  tmp_nog <- fread(Nog_table_list[i])
  setnames(tmp_nog, colnames(tmp_nog), c("ID_1", "ID_2", "organism_name", "organism_taxid", "other"))
  tmp_nog[, ID_sp := gsub('(^\\w+\\s)', "", tmp_nog[, organism_name])]
  tmp_nog[, ID_sp := str_sub(tmp_nog[, ID_sp], 1, 3)]
  tmp_nog[, ID_sp := str_c(str_sub(tmp_nog[, organism_name], 1, 1),".", tmp_nog[,ID_sp])]
  
  tmp_nog[, XP_NCBI := str_extract(tmp_nog[, other], "XP_\\d+\\.1" )]
  
  a <- gsub("(.txt)$", "", Nog_table_list[i])
  tmp_nog[, interest_gene := a]
  tmp_nog[, gene_id := a] 
  tmp_tab_nog <- rbind(tmp_tab_nog, tmp_nog, fill = TRUE)
  
}
id_list <- tmp_tab_nog$XP_NCBI

for (i in 1:length(id_list)) {
  #i <- 350
  if (!is.na(tmp_tab_nog$XP_NCBI[i])) {
    x <- efetch(tmp_tab_nog$XP_NCBI[i], "protein", rettype = "gp", retmode = "xml")
    write(content(x, "text"), file = "tmp.txt")
    text_id <- paste(fread('tmp.txt', fill = TRUE))
    id_of_gene <- str_extract_all(text_id, "GeneID:\\d+")[[1]]
    
    num_ids <- length(id_of_gene)
    
    if (num_ids > 1) {
      tmp_tab_nog[i, gene_id := ID_2]
    } else if (num_ids == 1) {
      tmp_tab_nog[i, gene_id := str_extract(id_of_gene[1], "\\d+")]
    } else {
      tmp_tab_nog[i, gene_id := NA]
    }
  }
  
  if (is.na(tmp_tab_nog[i, XP_NCBI])) {
    tmp_tab_nog[i, gene_id := ID_2]
  }
}

length(unique(tmp_tab_nog$gene_id))

################################################################################
#OrthoDB



Ortho_table_list <- file_list[str_detect(file_list, "[:digit:](?=_Metazoa.txt$)")]

tmp_tab_ortho <- data.table(interest_gene=character(),
                            organism_name= character(),
                            organism_taxid= character(),
                            pub_gene_id= character(),
                            corrected_pubgene = character(),
                            description= character(),
                            ID_sp= character())

for (i in 1:length(Ortho_table_list)) {
  #i <- 1
  tmp_orthodb <- fread(Ortho_table_list[i])
  corrected_pubgene <- str_replace_all(tmp_orthodb[, pub_gene_id], "\\\\", "_")
  corrected_pubgene <- str_replace_all(corrected_pubgene, "^\\w+\\d+;", "")
  tmp_orthodb <- cbind(tmp_orthodb,corrected_pubgene )
  tmp_a <- gsub("(_Metazoa.txt)$", "", Ortho_table_list[i])
  tmp_orthodb[, interest_gene := tmp_a]
  tmp_orthodb[, ID_sp := gsub('(^\\w+\\s)', "", tmp_orthodb[, organism_name])]
  tmp_orthodb[, ID_sp := str_sub(tmp_orthodb[, ID_sp], 1, 3)]
  tmp_orthodb[, ID_sp := str_c(str_sub(tmp_orthodb[, organism_name], 1, 1),".", tmp_orthodb[,ID_sp])]
  tmp_tab_ortho <- rbind(tmp_tab_ortho, tmp_orthodb[, .(interest_gene,corrected_pubgene, organism_name, organism_taxid,
                                                        pub_gene_id, description, ID_sp)])
}
  
  for (j in 462:length(tmp_tab_ortho$corrected_pubgene)) {
    #j <-  3
    ncbi_from_pubog <- esearch(paste(tmp_tab_ortho[j, 'organism_name'], 'and', tmp_tab_ortho[j, 'corrected_pubgene'], sep = ' '), "gene")
    
    write(content(ncbi_from_pubog, "text"), file = "tmp.txt")
    text_id <- paste(fread('tmp.txt', fill=TRUE))
    if (grepl("No items found\\.", text_id)) {
      tmp_tab_ortho[j, 'corrected_pubgene'] <- tmp_tab_ortho[j, 'pub_gene_id']
      next 
    }
    id_matches <- str_extract_all(text_id, "<Id>\\d+</Id>")[[1]]
    num_ids <- length(id_matches)
    if (num_ids > 1) {
      tmp_tab_ortho[j, 'corrected_pubgene'] <- tmp_tab_ortho[j, 'pub_gene_id']
    } else if (num_ids == 1) {
      id_of_gene <- str_extract(id_matches[1], "\\d+")
      tmp_tab_ortho[j, 'corrected_pubgene'] <- id_of_gene
    } else {
      tmp_tab_ortho[j, 'corrected_pubgene'] <- NA  # Or any default value
      
      if (is.na(tmp_tab_ortho[j, 'corrected_pubgene'])) {
        tmp_tab_ortho[j, 'corrected_pubgene'] <- tmp_tab_ortho[j, 'pub_gene_id']
      }
    }
    }

tmp_tab_ortho[, pub_gene_id := str_replace_all(pub_gene_id, '\\\\','')]
  
colnames(tmp_tab_ortho)[5] <- "gene_id"

flybase_to_cg_ortho <- fread("")#Table of ID and gene names
for (i in 1:nrow(flybase_to_cg_ortho)) {
  #i <- 1
  tmp_tab_ortho[interest_gene==flybase_to_cg_ortho$gene[i], interest_gene:= flybase_to_cg_ortho$ID[i]]
}


tmp_tab_nog[, XP_NCBI := ifelse(!grepl("^XP_", XP_NCBI), ID_1, XP_NCBI)]


#Combine data from OrthoDB and EggNOG

colnames(tmp_tab_nog)
colnames(tmp_tab_ortho)
colnames(tmp_tab_nog)[1] <- "pub_gene_id"
c <- colnames(tmp_tab_nog)[c(7,9,6, 3, 4,1 )]
ct <- colnames(tmp_tab_ortho)[c(1,5,7,2,3,4 )]
tmp_tab_ortho1 <- tmp_tab_ortho[, ..ct]
tmp_tab_nog1 <- tmp_tab_nog[, ..c]
#setnames(tabn, c, colnames(tab))

tab_NOGortho <- data.table(interest_gene=character(),
                           gene_id= character(),
                           ID_sp= character(),
                           organism_name= character(),
                           organism_taxid = character())
tab_NOGortho <- rbind(tmp_tab_ortho1, tmp_tab_nog1)
tab_NOGortho <- rbind(tmp_tab_ortho1)


fwrite( tab_NOGortho , file = "tab_NOGortho.txt", append = TRUE)
getwd()
#fwrite(tab, file="BTB_table_OrthoDB_EggNOG.txt", sep = "\t")


###################################################


file_list_code <- list.files(path="")
getwd()
fasta_true_list <- file_list_code[str_detect(file_list_code, "_statistics\\.fa$")]
fasta_lines <- readLines(fasta_true_list)
header_lines <- fasta_lines[grepl("^>", fasta_lines)]
fasta_ids <- sub("^>", "", header_lines)
fasta_ids_cleaned <- sub("/.*", "", fasta_ids)
fasta_list_cleaned <- data.table(name = fasta_ids_cleaned)

tab_NOGortho_upd <- copy(tab_NOGortho)
tab_NOGortho_upd[, fasta_id := paste0(ID_sp, "_", interest_gene, "_", gene_id, "_", pub_gene_id)]
tab_NOGortho_upd <- tab_NOGortho_upd[fasta_id %in% fasta_list_cleaned$name]


specie_check <- tab_NOGortho_upd[ID_sp == 'M.dom']
length(unique(specie_check$gene_id))
#################################################################################

length(unique(tab_NOGortho_upd$organism_name))


summary(grepl("Drosophila", unique(tab_NOGortho_upd$organism_name), ignore.case = TRUE, .N))
summary(grepl("Drosophila", unique(tab_NOGortho_upd$organism_name), ignore.case = FALSE, .N))


#List of species
species <- (unique(tab_NOGortho_upd$organism_name))
#taxize check specie name
specie_name_check <- gnr_resolve(species)

# retrieve class and order
get_class_and_order <- function(specie_name) {
  tryCatch({
    classification <- taxize::classification(specie_name, db = "ncbi")
    order <- classification[[1]]$name[classification[[1]]$rank == "order"]
    class <- classification[[1]]$name[classification[[1]]$rank == "class"]
    family <- classification[[1]]$name[classification[[1]]$rank == "family"]
    if (length(order) == 0) order <- NA
    if (length(class) == 0) class <- NA
    if (length(family) == 0) family <- NA
    return(list(Order = order, Class = class, Family = family))
  }, error = function(e) {
    return(list(Order = NA, Class = NA, Family = NA)) 
  })
}

classes_and_orders_ncbi <- lapply(species, get_class_and_order)

specie_and_order_ncbi <- data.table(
  Species = species,
  Order = sapply(classes_and_orders_ncbi, `[[`, "Order"),
  Class = sapply(classes_and_orders_ncbi, `[[`, "Class"),
  Family = sapply(classes_and_orders_ncbi, `[[`, "Family")
)

unique_count_orders <- specie_and_order_ncbi[, .(
  Count = .N,      
  Order = unique(Order)  
), by = Family] 

#################################################################################
classification <- taxize::classification('Contarinia nasturtii', db = "ncbi")

length(unique(tab_NOGortho$gene_id))
length(unique(tab_NOGortho$gene_id))


#################################################################################

species_list_simple <- c("Drosophila virilis","Drosophila melanogaster", "Ceratitis capitata","Musca domestica"
                         , "Homo sapiens",
                  "Lucilia sericata", 'Contarinia nasturtii', 'Culex pipiens','Danaus plexippus')

tax_info_simple <- classification(species_list_simple, db = "ncbi")

tree_from_taxonomy <- class2tree(tax_info_simple)
phylo_tree_simple <- tree_from_taxonomy$phylo
plot(phylo_tree_simple, cex = 0.8)
title("Taxonomy-Based Phylogenetic Tree")

ggtree_obj <- ggtree(phylo_tree_simple) + 
  theme_tree2() + 
  geom_tiplab(size = 3, angle = 0, hjust = 0) 
ggsave("phylogenetic_tree.svg", ggtree_obj, width = 20, height = 10, dpi = 300)
plot(ggtree_obj, cex =0.8)
#################################################################################
#Number of each protein in different species

filtered_data1 <- tab_NOGortho_upd[organism_name %in% species_list_simple]
selected_columns1<- c("interest_gene", "organism_name")

filtered_data1 <- filtered_data1[, ..selected_columns1]
unique_count_filtered_data1 <- filtered_data1[, .(
  Protein_Count = .N  
), by = .(organism_name, interest_gene)]


protein_matrix <- dcast(unique_count_filtered_data1, 
                        organism_name ~ interest_gene, value.var = "Protein_Count", fill = 0)
protein_matrix$Total_Proteins <- rowSums(protein_matrix[, -1], na.rm = TRUE)
protein_matrix$cg1603 <- 0




zero_row <- data.table(organism_name = c("Culex pipiens", "Danaus plexippus", "Homo sapiens"),
                       cg1603 = 0, cg1602 = 0, cg1605 = 0, cg2129 = 0, 
                       cg10959 = 0, cg18262 = 0, cg8643 = 0, cg8944 = 0, 
                       Total_Proteins = 0)

protein_matrix <- rbind(protein_matrix, zero_row, fill = T)


desired_column_order <- c("organism_name", "cg1602", "cg1603","cg1605", "cg2129", 
                          "cg10959", "cg18262", 'cg8643', 'cg8944', 'Total_Proteins') 
setcolorder(protein_matrix, desired_column_order)


# Rearrange rows (organism names) 

desired_row_order <- c("Drosophila virilis","Drosophila melanogaster", "Ceratitis capitata","Musca domestica", "Lucilia sericata",
                       'Culex pipiens', 'Contarinia nasturtii', 'Danaus plexippus', 'Homo sapiens')

protein_matrix[, organism_name := factor(organism_name, levels = desired_row_order)]
setorder(protein_matrix, organism_name)

protein_long <- melt(protein_matrix, id.vars = "organism_name", 
                     variable.name = "Protein", value.name = "Count")
plot <- ggplot(protein_long, aes(x = Protein, y = organism_name, fill = Count)) +
  geom_tile(color = "white") +  
  geom_text(aes(label = Count), size = 3) +  # Add counts inside the tiles
  scale_fill_gradient(low = "white", high = "steelblue") + 
  theme_minimal() +  
  labs(title = "Protein Count in Species",
       x = "Protein",
       y = "Species") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, family = "Helvetica"),  # Rotate protein names for readability
        axis.text.y = element_text(size = 10, family = "Helvetica"),
        axis.title = element_text(size = 12),family = "Helvetica",
        plot.title = element_text(family = "Helvetica")) +
  scale_y_discrete(limits = rev(levels(factor(protein_long$organism_name))))


ggsave("protein_count_heatmap.svg", plot = plot, width = 10, height = 7, dpi = 300)


########################################################################################
#Filter fasta sequences 

general_fasta_all_sp <- readAAStringSet('/Users/kostya/Documents/bioinformatics/Phylogenetics_ZNFs/cg_proteins(true_folder)/results/full_list/general_fasta_all_sp.fa')
sequence_ids <- names(general_fasta_all_sp)
fasta_list_cleaned1 <- str_replace_all(fasta_list_cleaned$name,"^>",'')
filtered_sequences <- general_fasta_all_sp[sequence_ids %in% fasta_list_cleaned1]
writeXStringSet(filtered_sequences, 'fasta_all_sp_cleaned.fa')


library(rotl)
library(ape)
library(ggtree)


families_grapf_abstr <- c("Cecidomyiidae", "Glossinidae", "Calliphoridae", 
              "Muscidae", "Tephritidae", "Drosophilidae", "Diopsidae")

matched <- tnrs_match_names(families_grapf_abstr)


tree <- tol_induced_subtree(ott_ids = matched$ott_id)

plot(tree, main = "Phylogenetic Tree of Selected Diptera Families")



