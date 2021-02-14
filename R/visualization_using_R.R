#------------------------------------------------------------------
# visualization using R platform
#------------------------------------------------------------------
#' function to define the node in a R graph
#'
#'
#' @param rxn_core_carbon_inf a dataframe contains the detailed annotation of rxn_core
#' @param currency_metabolites_inf a vector of currency metabolites
#'
#' @return
#' @export node_all a dataframe contains the node information for a R graph
#'
#' @examples
creatRnode <- function(rxn_core_carbon_inf, currency_metabolites_inf){

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




#' function to define the edge
#'
#'
#' @param rxn_core_carbon_inf a dataframe contains the detailed annotation of rxn_core
#' @param rxn_rev_list a vector contains the reversible reactions
#' @param node_inf a dataframe contains the node definition for a R graph
#'
#' @return
#' @export edge4 a dataframe defines the edge between the nodes
#'
#' @examples
creatRedge <- function(rxn_core_carbon_inf, rxn_rev_list, node_inf){

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




#' function to creat the R graph based on the R platform
#'
#'
#' @param nodes a dataframe contains the id and label information of each node
#' @param edges a dataframe defining the connection of the metabolite
#'
#' @return
#' @export i_graph_x3 a graph format in R
#'
#' @examples
creatRmap <- function (nodes, edges){

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
