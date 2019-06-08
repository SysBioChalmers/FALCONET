# This code is used to prepare the data for the map from model in BIGG format
#source('model change.R')
#source('visualization_based_R.R')

library(FALCONET)
library(stringr)
library(tidyverse)
library(readxl)
library(igraph)
library(networkD3)
library(hongR)
library(DiagrammeR)

#part one
#prepare the reaction format
rxn <- read_excel("data/yeast_core_model.xlsx", sheet = "Reaction List")
metabolite <-  read_excel("data/yeast_core_model.xlsx", sheet = "Metabolite List")

# Unify the metabolite
metabolite <- select(metabolite, `Metabolite name`, `Metabolite formula`, Charge, KEGGID)
colnames(metabolite) <- c('Metabolite description', 'Metabolite formula', 'Charge', 'KEGGID')
#prepare the standard compartment
comparment <- unlist(str_extract_all(metabolite$`Metabolite description`, "_[:alpha:]$")) %>%
  str_replace_all(.,"_","[")
comparment <- paste(comparment, "]", sep = "")
for(i in seq_along(comparment)){
metabolite$`Metabolite description`[i] <- str_replace_all(metabolite$`Metabolite description`[i],"_[:alpha:]$", comparment[i])
}
metabolite$`Metabolite description` <- str_replace_all(metabolite$`Metabolite description`,"M_","")



rxn0 <- select(rxn, Abbreviation, Description_new)
colnames(rxn0) <- c('ID0', 'Equation')
rxn_split <- splitRxnToMetabolite(rxn0, sep0 = "<=>")
#first remove the exchange reaction
exchange_index <- which(is.na(rxn_split$MetID))
exchange_id <- rxn_split$ID[exchange_index]
others <- which(rxn_split$ID %in% exchange_id ==FALSE)
rxn_split_refine <- rxn_split[others,]
rxn_split_refine <- select(rxn_split_refine, ID, MetID, compostion)
colnames(rxn_split_refine) <- c('v2','v3','type')
rxn_split_refine$subsystem <- getSingleReactionFormula(rxn$Subsystem_new,rxn$Abbreviation,rxn_split_refine$v2)
rxn_split_refine$v2 <- str_replace_all(rxn_split_refine$v2,"R_","r_")

#prepare the standard compartment
comparment1 <- unlist(str_extract_all(rxn_split_refine$v3, "_[:alpha:]$")) %>%
  str_replace_all(.,"_","[")
comparment1 <- paste(comparment1, "]", sep = "")
for(i in seq_along(comparment1)){
  rxn_split_refine$v3[i] <- str_replace_all(rxn_split_refine$v3[i],"_[:alpha:]$", comparment1[i])
}
rxn_split_refine$v3 <- str_replace_all(rxn_split_refine$v3, "M_", "")



# choose the subsytem
subsystem1 <-  "sce00010  Glycolysis / Gluconeogenesis"
# for test we change the subsystem into 'one'
rxn_split_refine$subsystem <- "one"
subsystem1 <- "one"
# Define the currency metabolite in each subsystem
currency_metabolites <- DefineCurrencyMet(rxn_split_refine, 
                                          subsystem0=subsystem1,
                                          numberGEM=14,
                                          numberSubsystem=1)

# remove the reactions with only one metabolite
# if we do not remove the currency metabolite in the model then this step is mainly removed exchange reaction
rxn_split_refine <- removeRxnWithSingleMet(rxn_split=rxn_split_refine)


#--------------------------------------------------------------------------------------------
## define the base reactant and product for cellDesigner
#---------------------------------------------------------------------------------------------
rxn_split_refine <- addBaseTypeIntoRxn(rxn_split_refine, metabolite, currency_metabolites)

#---------------------------------------------------
# choose reaction based on the subsystem
#-----------------------------------------------------
rxn_core_carbon <-chooseRxnFromSubsystem(rxn_split_refine, subsystem0=subsystem1) 
rxnID_choose <- unique(rxn_core_carbon$v2)

#------------------------------------------------------------------
# produce the met, rxn and gpr used for the map production
#------------------------------------------------------------------
# prepare the metabolites formula
# this funcion is used to prepare the metabolite annotation for cell designer
met_annotation <- prepareMET(rxn_core_carbon, currency_metabolites,rxnID_choose)

# prepare the rxn formula
rxn_core_carbon_cellD0 <- prepareRXN(rxn_core_carbon,met_annotation,currency_metabolites)









#------------------------------------------------------------------
# visualization using R platform
#------------------------------------------------------------------
# creat nodes
currency_metabolite_R_map <- c("h", "h2o", "atp", "pi", "adp", "nadp", "nadph", "nad", "nadh",  "ppi",  "amp" )
node01 <- creatRnode(rxn_core_carbon_inf = rxn_core_carbon_cellD0, currency_metabolites_inf = currency_metabolite_R_map)
# creat edge
# prepare the reversible reaction list
rxn_rev <- filter(rxn, IsReversible==TRUE)
rxn_rev <- rxn_rev$Abbreviation
edge01 <- creatRedge(rxn_core_carbon_inf = rxn_core_carbon_cellD0, rxn_rev_list = rxn_rev, node_inf = node01)
# creat map
i_graph_3 <- creatRmap(nodes=node01, edges=edge01)
# visualization
render_graph(i_graph_3, layout = 'nicely', output = "visNetwork")
# visualization-general
render_graph(i_graph_3, layout = 'nicely')
# export the graph
i_graph_3 %>%
  export_graph(
    file_name = "result/graph.eps",
    title = "Simple Graph")









# note:
# find the path between two nodes
# for the big graph, this steps can be slow.
#i_graph_3 %>%
#  get_paths(
#    from = 54,
#    to = 49)