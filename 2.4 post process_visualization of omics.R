# An example to show how we can map the fluxes on to the map
# the xml file obtained based on the automatic layout.
# otherwise the followed code can't be used

# function to define the color based on the fold change
defineGeneFoldColor <- function(omic_fold, up=2, down=0.5) {
  color <- vector()
  for (i in seq_along(omic_fold)) {
    if (omic_fold[i] >= up) {
      color[i] <- "ffff0000" # red
    } else if (omic_fold[i] < down) {
      color[i] <- "ff0eb10e" # green
    } else {
      color[i] <- "ffccff66"
    }
  }
  
  return(color)
}



yeast_map <- readLines(file("result/model_test_check.xml"))
index_rxn <- which(str_detect(yeast_map, "reaction metaid"))

# extract the rxns
rxnid <- yeast_map[index_rxn]
rxnid0 <- str_split(rxnid, " ")
rxnid1 <- vector()
for (i in 1:length(rxnid0)){
  rxnid1[i] <- rxnid0[[i]][3]
}
rxnid1 <- str_replace_all(rxnid1, "id=","") %>%
  str_replace_all(.,"\"","")

flux_data <- data.frame(rxnid = rxnid1, stringsAsFactors = FALSE)
flux_data$gene_id <- str_replace_all(flux_data$rxnid, 'r_','g_')
flux_data$protein_id <- str_replace_all(flux_data$rxnid, 'r_','p_')
flux_data$g_value <- getSingleReactionFormula(flux_map$gene_fold, flux_map$Abbreviation,flux_data$rxnid)
flux_data$p_value <- getSingleReactionFormula(flux_map$protein_fold, flux_map$Abbreviation,flux_data$rxnid)

color_gene <- defineGeneFoldColor(omic_fold = flux_data$g_value, up=1.5, down = 0.6)
color_protein <- defineGeneFoldColor(omic_fold = flux_data$p_value)


#"ff0eb10e\" scheme=\"Color\"/>"   # green
#ffff0000\" scheme=\"Color\"/>"   # red
#ffccff66 # original background



# for gene
# firstly we get the coordinate of gene and protein on the map based on gene_id and protein_id
gene0 <- flux_data$gene_id
gene_index <- vector()
for(i in seq_along(gene0)){
  gene_index[i] <- which(str_detect(yeast_map, "<species") & str_detect(yeast_map, gene0[i]))
}

gene_mapid <- yeast_map[gene_index]
gene_mapid0 <- str_split(gene_mapid, " ")

gene_mapid1 <- vector()
for (i in 1:length(gene_mapid0)){
  gene_mapid1[i] <- gene_mapid0[[i]][3]
}

gene_mapid1 <- str_replace_all(gene_mapid1, "id=","") %>%
  str_replace_all(.,"\"","")

# get the coordinate of gene
gene_position <- vector()
for (i in seq_along(gene_mapid1)){
  gene_position[i] <- which(str_detect(yeast_map,"<celldesigner:speciesAlias") & str_detect(yeast_map, gene_mapid1[i]))
}

# change the color
color_position <- gene_position+9 #fill color

new_color <- paste("<celldesigner:paint color=\"",color_gene,"\" scheme=\"Color\"/>", sep = "")
yeast_map[color_position] <- new_color







# for protein
# firstly we get the coordinate of protein and protein on the map based on protein_id and protein_id
protein0 <- flux_data$protein_id
protein_index <- vector()
for(i in seq_along(protein0)){
  protein_index[i] <- which(str_detect(yeast_map, "<species") & str_detect(yeast_map, protein0[i]))
}

protein_mapid <- yeast_map[protein_index]
protein_mapid0 <- str_split(protein_mapid, " ")

protein_mapid1 <- vector()
for (i in 1:length(protein_mapid0)){
  protein_mapid1[i] <- protein_mapid0[[i]][3]
}

protein_mapid1 <- str_replace_all(protein_mapid1, "id=","") %>%
  str_replace_all(.,"\"","")

# get the coordinate of protein
protein_position <- vector()
for (i in seq_along(protein_mapid1)){
  protein_position[i] <- which(str_detect(yeast_map,"<celldesigner:speciesAlias") & str_detect(yeast_map, protein_mapid1[i]))
}

# change the color
color_position <- protein_position+9 #fill color
new_color <- paste("<celldesigner:paint color=\"",color_protein,"\" scheme=\"Color\"/>", sep = "")
yeast_map[color_position] <- new_color



# change all other background into the same color
yeast_map <- str_replace_all(yeast_map, "ffccff66", "ff33ffff")


writeLines(yeast_map, file("result/model_test_check.xml"))