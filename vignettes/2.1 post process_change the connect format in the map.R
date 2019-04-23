## these code can be used adjust the line and arraw connect in map of cellDesigner4.4

## only change the arraw connection
yeast_map <- readLines(file("result/model_test_check.xml"))
anchor <- readLines(file("data/anchor.xml"))
index_linkAnchor <- which(str_detect(yeast_map,"celldesigner:linkAnchor")==FALSE)
yeast_map_new <- yeast_map[index_linkAnchor]
writeLines(yeast_map_new, file("result/model_test_protein-add rxn connection_using function_auto_layout_remove the linkAnchor.xml")) # save


# only change the line direction between two metabolites
yeast_map <- readLines(file("result/model_test_check.xml"))
index_linkAnchor_arrow <- which(str_detect(yeast_map,"celldesigner:editPoints")==TRUE)
index_linkAnchor_arrow1 <- which(str_detect(yeast_map,"<celldesigner:lineDirection index=\"1\" value=\"unknown\"/>")==TRUE)
index_linkAnchor_arrow2 <- which(str_detect(yeast_map,"<celldesigner:lineDirection index=\"2\" value=\"unknown\"/>")==TRUE)
index_arrow0 <- c(index_linkAnchor_arrow,index_linkAnchor_arrow1,index_linkAnchor_arrow2)
yeast_map_new1 <- yeast_map[-index_arrow0]
writeLines(yeast_map_new1, file("result/model_test_protein-add rxn connection_using function_auto_layout_remove the linkAnchor.xml")) # save



#previous code
#yeast_map <- readLines(file("result/model_with_protein-add rxn connection.xml"))
#anchor <- readLines(file("data/anchor.xml"))
#index_linkAnchor <- which(str_detect(yeast_map,"celldesigner:linkAnchor")==FALSE)
#yeast_map_new <- yeast_map[index_linkAnchor]
#index_linkAnchor_arrow <- which(str_detect(yeast_map_new,"celldesigner:editPoints")==TRUE)
#index_linkAnchor_arrow1 <- which(str_detect(yeast_map_new,"<celldesigner:lineDirection index=\"1\" value=\"unknown\"/>")==TRUE)
#index_linkAnchor_arrow2 <- which(str_detect(yeast_map_new,"<celldesigner:lineDirection index=\"2\" value=\"unknown\"/>")==TRUE)
#index_arrow0 <- c(index_linkAnchor_arrow,index_linkAnchor_arrow1,index_linkAnchor_arrow2)
#yeast_map_new1 <- yeast_map_new[-index_arrow0]
#writeLines(yeast_map_new1, file("result/model_test_protein-add rxn connection_using function_auto_layout_remove the linkAnchor.xml")) # save