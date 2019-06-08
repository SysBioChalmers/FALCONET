# An example to show how we can map the fluxes fold changes on to the map
# the xml file obtained based on the automatic layout.
# otherwise the followed code can't be used
# input the flux data

flux_map <- read_excel("data/flux_map.xlsx")
fluxFoldMapping(input_map ="result/model_test_check.xml",
            flux_inf = flux_map,
            output_map = "result/map_with_flux_fold_changes.xml")
