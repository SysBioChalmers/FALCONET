library(hongR)
library(tidyverse)
## these code can be used the medium process of draft map after automatic layout
yeast_map <- readLines(file("result/model_with_protein-add rxn connection_auto_layout_adjustment.xml"))
yeast_map <- readLines(file("result/model_with_protein-add rxn connection.xml"))

anchor <- readLines(file("data/anchor.xml"))
index_linkAnchor <- which(str_detect(yeast_map,"celldesigner:linkAnchor")==FALSE)
yeast_map_new <- yeast_map[index_linkAnchor]

index_linkAnchor_arrow <- which(str_detect(yeast_map_new,"celldesigner:editPoints")==TRUE)
index_linkAnchor_arrow1 <- which(str_detect(yeast_map_new,"<celldesigner:lineDirection index=\"1\" value=\"unknown\"/>")==TRUE)
index_linkAnchor_arrow2 <- which(str_detect(yeast_map_new,"<celldesigner:lineDirection index=\"2\" value=\"unknown\"/>")==TRUE)
index_arrow0 <- c(index_linkAnchor_arrow,index_linkAnchor_arrow1,index_linkAnchor_arrow2)
yeast_map_new1 <- yeast_map_new[-index_arrow0]

writeLines(yeast_map_new1, file("result/model_test_protein-add rxn connection_using function_auto_layout_remove the linkAnchor.xml")) # save
##others: update the reaction reversibility




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
metabolite_define$mettype <- getMultipleReactionFormula(rxn_define$note,rxn_define$id_mapping,metabolite_define$id_mapping)

# here we can define the currency metabolite and non currency metabolite
# choose the none-base metabolites
metabolite_choose <- filter(metabolite_define, mettype !="base" | type =='PROTEIN' | type=='GENE')
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

writeLines(yeast_map, file("result/model_test_check.xml")) # save



















# test code
yeast_map <- readLines(file("result/model_test.xml"))
# example
# sa26, s26
index <- which(str_detect(yeast_map,"celldesigner:speciesAlias") & str_detect(yeast_map,"sa26") & str_detect(yeast_map,"s26"))
# bounds line
pos_bounds <- index + 2
# size line
pos_size <- index + 3
# box size line
pos_box <- index + 7
# change the size of metabolite component in the graph
yeast_map[pos_bounds] <- str_replace_all(yeast_map[pos_bounds], "w=\"70.0\" h=\"25.0\"/>", "w=\"35.0\" h=\"12.0\"/>")
yeast_map[pos_size] <- str_replace_all(yeast_map[pos_size], "<celldesigner:font size=\"12\"/>", "<celldesigner:font size=\"8\"/>")
yeast_map[pos_box] <- str_replace_all(yeast_map[pos_box], "<celldesigner:boxSize width=\"70.0\" height=\"25.0\"/>", "<celldesigner:boxSize width=\"35.0\" height=\"12.0\"/>")

writeLines(yeast_map, file("result/model_test_check.xml")) # save
##others: update the reaction reversibility

