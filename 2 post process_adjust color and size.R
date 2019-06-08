## change the size for the none-base metabolites
## define a list of none-base metabolites list
## adjust the size based on none-base metabolite list
# define the none-base metabolite
# test code


# step 1 produce the map
source('1.1 preprocess metabolite and reaction for yeastGEM.R')



# step 2 manually adjust the layout of map in CellDesigner




# step 3 adjust the size or color of components in the map
yeast_map <- readLines(file("result/model_with_protein-add rxn connection.xml"))
metabolite_define <- met_annotation
metabolite_define$id_mapping <- paste(metabolite_define$rxnID, metabolite_define$name, sep = "@")
rxn_define <- rxn_core_carbon
rxn_define$id_mapping <- paste(rxn_define$v2, rxn_define$v3, sep = "@")
metabolite_define$name0 <- str_replace_all(metabolite_define$name, "\\[.*?\\]", "")
metabolite_define$currency <- NA
metabolite_define$currency[metabolite_define$name0 %in% currency_metabolites] <-'YES'

# here we can define the currency metabolite and non currency metabolite
# choose the none-base metabolites
metabolite_choose <- filter(metabolite_define, currency == "YES" | type =='PROTEIN' | type=='GENE')


# batch process
for (i in 1:nrow(metabolite_choose)){
  print(i)
  print(metabolite_choose$species[i])
  yeast_map <- changeFontSize(positionID = metabolite_choose$id[i],
                              uniqueID = metabolite_choose$species[i],
                              onemap = yeast_map)
}

writeLines(yeast_map, file("result/model_test_check.xml"))