library(hongR)
library(tidyverse)

## change the size for the none-base metabolites
## define a list of none-base metabolites list
## adjust the size based on none-base metabolite list
# define the none-base metabolite
# test code
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
changeFontSize <- function(positionID, uniqueID, onemap) {
  index <- which(str_detect(onemap, "celldesigner:speciesAlias") & str_detect(onemap, positionID) & str_detect(onemap, uniqueID))
  # bounds line
  pos_bounds <- index + 2
  # size line
  pos_size <- index + 3
  # box size line
  pos_box <- index + 7
  # change the size of metabolite component in the graph
  onemap[pos_bounds] <- str_replace_all(onemap[pos_bounds], "w=\"70.0\" h=\"25.0\"/>", "w=\"35.0\" h=\"12.0\"/>")
  onemap[pos_size] <- str_replace_all(onemap[pos_size], "<celldesigner:font size=\"12\"/>", "<celldesigner:font size=\"8\"/>")
  onemap[pos_box] <- str_replace_all(onemap[pos_box], "<celldesigner:boxSize width=\"70.0\" height=\"25.0\"/>", "<celldesigner:boxSize width=\"35.0\" height=\"12.0\"/>")
  return(onemap)

}

# batch process
for (i in 1:nrow(metabolite_choose)){
  print(i)
  print(metabolite_choose$species[i])
  yeast_map <- changeFontSize(positionID = metabolite_choose$id[i],
                              uniqueID = metabolite_choose$species[i],
                              onemap = yeast_map)
}

writeLines(yeast_map, file("result/model_test_check.xml"))