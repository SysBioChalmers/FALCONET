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
# function to define the node in a R graph
creatRnode <- function(rxn_core_carbon_inf, currency_metabolites_inf){
  #input
  #rxn_core_carbon_inf a dataframe contains the detailed annotation of rxn_core
  #currency_metabolites_inf a vector of currency metabolites
  #output
  #node_all a dataframe contains the node information for a R graph
  
  rxn_core_carbon_inf$rxnID <- str_replace_all(rxn_core_carbon_inf$rxnID, "r_", "R_")
  # analysis the metabolite and remove these currency metabolites
  # choose the rxn based met 
  met_analysis <- rxn_core_carbon_inf %>%
    count(name) %>%
    arrange(., desc(n)) 
  #currency_metabolite_R_map <- c("h", "h2o", "atp", "pi", "adp", "nadp", "nadph", "nad", "nadh",  "ppi",  "amp" )
  met_analysis$name_no_compartment <- str_replace_all(met_analysis$name, "\\[.*?\\]", "")
  met_analysis$currency_sign <- met_analysis$name_no_compartment %in% currency_metabolites_inf
  met_analysis <- filter(met_analysis, currency_sign==FALSE)
  
  
  # define the node
  label_all <- unique(c(rxn_core_carbon_inf$rxnID, met_analysis$name))
  id_all <- 1:length(label_all)
  
  node_all <- data.frame(id=id_all, label=label_all, stringsAsFactors = FALSE)
  
  from0 <- vector()
  to0 <- vector()
  for (i in seq_along(rxn_core_carbon_inf$rxnID)){
    print(i)
    if(rxn_core_carbon_inf$type[i] =='reactant'){
      from0[i] <- rxn_core_carbon_inf$name[i]
      to0[i] <- rxn_core_carbon_inf$rxnID[i]
    } else{
      from0[i] <- rxn_core_carbon_inf$rxnID[i]
      to0[i] <- rxn_core_carbon_inf$name[i]
    }
  }
  
  return(node_all)
  
}

# function to define the edge
creatRedge <- function(rxn_core_carbon_inf, rxn_rev_list, node_inf){
  #input
  #rxn_core_carbon_inf a dataframe contains the detailed annotation of rxn_core
  #rxn_rev_list, a vector contains the reversible reactions
  #node_inf, a dataframe contains the node definition for a R graph
  #output
  #edge4, a dataframe defines the edge between the nodes
  rxn_core_carbon_inf$rxnID <- str_replace_all(rxn_core_carbon_inf$rxnID, "r_", "R_")
  from0 <- vector()
  to0 <- vector()
  for (i in seq_along(rxn_core_carbon_inf$rxnID)){
    print(i)
    if(rxn_core_carbon_inf$type[i] =='reactant'){
      from0[i] <- rxn_core_carbon_inf$name[i]
      to0[i] <- rxn_core_carbon_inf$rxnID[i]
    } else{
      from0[i] <- rxn_core_carbon_inf$rxnID[i]
      to0[i] <- rxn_core_carbon_inf$name[i]
    }
  }
  
  edge0 <- data.frame(from1 = from0, to1 = to0, stringsAsFactors = FALSE)
  
  # add the reversiblity reaction
  rxn_rev <- rxn_rev_list
  edge0$rev_sign1 <- edge0$from1 %in% rxn_rev
  edge0$rev_sign2 <- edge0$to1 %in% rxn_rev
  edge2 <- filter(edge0, rev_sign1 ==TRUE | rev_sign2==TRUE)
  edge2_o <-  edge2[,1:2]
  edge2_n <- edge2_o
  edge2_n$from1 <- edge2_o$to1
  edge2_n$to1 <- edge2_o$from1
  
  edge3 <- rbind.data.frame(edge0[,1:2], edge2_n)
  edge3$from <- getSingleReactionFormula(node_inf$id,node_inf$label,edge3$from1)
  edge3$to <- getSingleReactionFormula(node_inf$id,node_inf$label,edge3$to1)
  
  edge_all <- edge3[,c('from','to')]
  edge_all$from <- as.numeric(edge_all$from)
  edge_all$to <- as.numeric(edge_all$to)
  return(edge_all)
}

# function to creat the R graph based on the R platform
creatRmap <- function (nodes, edges){
  # input
  # nodes, a dataframe contains the id and label information of each node
  # edges, a dataframe defining the connection of the metabolite
  # output
  # i_graph_x3, a graph format in R
  
  # Create the graph object
  i_graph_x1 <-
    create_graph()
  # Add the nodes to the graph
  i_graph_x2 <-
    i_graph_x1 %>%
    add_nodes_from_table(
      table = nodes,
      label_col = label)
  # Add the edges to the graph
  i_graph_x3 <-
    i_graph_x2 %>%
    add_edges_from_table(
      table = edges,
      from_col = from,
      to_col = to,
      from_to_map = id_external)
  return(i_graph_x3)
}


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

# find the path between two nodes
# for the big graph, this steps can be slow.
i_graph_3 %>%
  get_paths(
    from = 54,
    to = 49)








# Test code
# -------------------------------------------------------------------------------
rxn_core_carbon_cellD0$rxnID <- str_replace_all(rxn_core_carbon_cellD0$rxnID, "r_", "R_")
# analysis the metabolite and remove these currency metabolites
# choose the rxn based met 
met_analysis <- rxn_core_carbon_cellD0 %>%
  count(name) %>%
  arrange(., desc(n)) 
currency_metabolite_R_map <- c("h", "h2o", "atp", "pi", "adp", "nadp", "nadph", "nad", "nadh",  "ppi",  "amp" )
met_analysis$name_no_compartment <- str_replace_all(met_analysis$name, "\\[.*?\\]", "")
met_analysis$currency_sign <- met_analysis$name_no_compartment %in% currency_metabolite_R_map
met_analysis <- filter(met_analysis, currency_sign==FALSE)


# define the node
label_all <- unique(c(rxn_core_carbon_cellD0$rxnID, met_analysis$name))
id_all <- 1:length(label_all)

node0 <- data.frame(id=id_all, label=label_all, stringsAsFactors = FALSE)



# define the edge
from0 <- vector()
to0 <- vector()
for (i in seq_along(rxn_core_carbon_cellD0$rxnID)){
  print(i)
  if(rxn_core_carbon_cellD0$type[i] =='reactant'){
    from0[i] <- rxn_core_carbon_cellD0$name[i]
    to0[i] <- rxn_core_carbon_cellD0$rxnID[i]
  } else{
    from0[i] <- rxn_core_carbon_cellD0$rxnID[i]
    to0[i] <- rxn_core_carbon_cellD0$name[i]
  }
}

edge0 <- data.frame(from1 = from0, to1 = to0, stringsAsFactors = FALSE)

# add the reversiblity reaction
rxn_rev <- filter(rxn, IsReversible==TRUE)
rxn_rev <- rxn_rev$Abbreviation

edge0$rev_sign1 <- edge0$from1 %in% rxn_rev
edge0$rev_sign2 <- edge0$to1 %in% rxn_rev
edge2 <- filter(edge0, rev_sign1 ==TRUE | rev_sign2==TRUE)
edge2_o <-  edge2[,1:2]
edge2_n <- edge2_o
edge2_n$from1 <- edge2_o$to1
edge2_n$to1 <- edge2_o$from1

edge3 <- rbind.data.frame(edge0[,1:2], edge2_n)
edge3$from <- getSingleReactionFormula(node0$id,node0$label,edge3$from1)
edge3$to <- getSingleReactionFormula(node0$id,node0$label,edge3$to1)

edge4 <- edge3[,c('from','to')]
edge4$from <- as.numeric(edge4$from)
edge4$to <- as.numeric(edge4$to)

# creat the graph
# Create the graph object
i_graph_1 <-
  create_graph()
# Add the nodes to the graph
i_graph_2 <-
  i_graph_1 %>%
  add_nodes_from_table(
    table = node0,
    label_col = label)
# Add the edges to the graph
i_graph_3 <-
  i_graph_2 %>%
  add_edges_from_table(
    table = edge4,
    from_col = from,
    to_col = to,
    from_to_map = id_external)