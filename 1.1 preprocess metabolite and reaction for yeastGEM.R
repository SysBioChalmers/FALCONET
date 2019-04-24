# This code is used to prepare the data for the map
# source('model change.R')
# source('transition for cellDesigner.R')

library(fastgraphGEM)
library(stringr)
library(tidyverse)
library(readxl)
library(igraph)
library(networkD3)
library(hongR)

# prepare the reaction format
rxn <- read_excel("data/yeastGEM_latest version1.xls", sheet = "Reaction List")
metabolite <-  read_excel("data/yeastGEM_latest version1.xls", sheet = "Metabolite List")
# Update the metabolite name in rxn sheet
rxn_split <- splitAndCombine0(rxn$Reaction,rxn$Abbreviation,sep=" ")
#using the detailed name to take place of the short name
rxn_split$v3 <- getSingleReactionFormula(metabolite$Description,metabolite$Abbreviation,rxn_split$v1)
for (i in 1:length(rxn_split$v2)){
  if(rxn_split$v3[i]=="NA"){
    rxn_split$v3[i] <- rxn_split$v1[i]
  } else{
    rxn_split$v3[i] <- rxn_split$v3[i]
  }
  
}
rxn_split$v3 <- str_replace_all(rxn_split$v3, " \\[.*?\\]", "")
rxn_split$v3 <- str_trim(rxn_split$v3, side = "both")
rxn_split$compartment <- str_extract(rxn_split$v1, "\\[.*?\\]")
for (i in 1:nrow(rxn_split)){
  if(!is.na(rxn_split$compartment[i])){
    rxn_split$v3[i] <- paste(rxn_split$v3[i],rxn_split$compartment[i], sep = "")
  } else{
    rxn_split$v3[i] <- rxn_split$v3[i]
  }
}
rxn$Description_new <- getMultipleReactionFormula(rxn_split$v3,rxn_split$v2,rxn$Abbreviation)
rxn$Description_new <- str_replace_all(rxn$Description_new,";"," ")
# Update the subsytem in rxn sheet to make sure each reaction belong to one subsystem
subsystem_V3 <- read_excel("data/subsystem for yeast8 map_V3.xlsx")
rxn$Subsystem_new <- getSingleReactionFormula(subsystem_V3$subsystem_map_v3,subsystem_V3$rxnID,rxn$Abbreviation)

# filter rxn based on the fluxes
# rxn <- filter(rxn, Flux != 0) # we can remove the reaction with zero fluxes
# which(rxn$Abbreviation =='r_2034')


# Update the metabolite formula
metabolite$Description <- str_replace_all(metabolite$Description, " \\[.*?\\]", "")
metabolite$Description <- str_trim(metabolite$Description, side = "both")
metabolite$compartment <- str_extract(metabolite$Abbreviation, "\\[.*?\\]")
for (i in 1:nrow(metabolite)){
  if(!is.na(metabolite$compartment[i])){
    metabolite$Description[i] <- paste(metabolite$Description[i],metabolite$compartment[i], sep = "")
  } else{
    metabolite$Description[i] <- metabolite$Description[i]
  }
}




# prepare the rxn format for the cellDesigner
colnames(metabolite) <-  c('Metabolite name','Metabolite description','Metabolite formula','Charge','Compartment','KEGGID','CHEBI')
rxn_split_refine <- splitRxnToMetabolite.Yeast(rxn, metabolite)

# analysis subsystem
analysis_subsystem <- rxn %>%
  count(Subsystem_new) %>%
  arrange(., desc(n)) 

# choose the subsytem
subsystem1 <- c("glycolysis / gluconeogenesis \\( sce00010 \\)")

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
rxn_core_carbon <- chooseRxnFromSubsystem_new(rxn_split_refine_inf = rxn_split_refine, subsystem0 = subsystem1)
rxnID_choose <- unique(rxn_core_carbon$v2)


#------------------------------------------------------------------
# produce the met, rxn and gpr used for the map production
#------------------------------------------------------------------
# prepare the metabolites formula
# this funcion is used to prepare the metabolite annotation for cell designer
met_annotation <- prepareMET(rxn_core_carbon, currency_metabolites,rxnID_choose)
# prepare the rxn formula
rxn_core_carbon_cellD0 <- prepareRXN(rxn_core_carbon,met_annotation,currency_metabolites)
# prepare the protein and gene
gpr <- prepareGPR(met_annotation)

#save the exampel data format for cell designer
#write.table(met_annotation,"result/met_annotation for example.txt", row.names = FALSE, sep = "\n")
#write.table(rxn_core_carbon_cellD0,"result/rxn_core_carbon_cellD0 for example.txt", row.names = FALSE, sep = "\n")
#write.table(gpr,"result/gpr for example.txt", row.names = FALSE, sep = "\n")

#------------------------------------------------------------------
# produce the file as the input for the cellDesigner
#------------------------------------------------------------------
produceInputForCellDesigner(met_annotation, 
                            gpr,
                            rxn_core_carbon_cellD0,
                            x_size=1200, 
                            y_size=2000)
