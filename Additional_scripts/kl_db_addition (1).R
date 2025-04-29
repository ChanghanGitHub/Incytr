#LOAD ORIGINAL DB

setwd("/Users/csimpson/Desktop/Collaborations/Exfinder/Incytr_script_081924/Database")
# load the exFINDER database (can also download the updated database from the exFINDER github page)
load("DB_Layer1_mouse_filtered.rda")
load("DB_Layer2_mouse_filtered.rda")
load("DB_Layer3_mouse_filtered.rda")


#LOAD KL DATA
kl_call_siteQuant_19_export <- read_csv("~/Desktop/Collaborations/Exfinder/mc38/input_data/deconvolution_parsed/kl_call_siteQuant_19_export.csv")
kl_call_siteQuant_22_export2 <- read.csv("~/Desktop/Collaborations/Exfinder/mc38/input_data/deconvolution_parsed/kl_call_siteQuant_22_export2.csv")
names(kl_call_siteQuant_22_export2) <- names(kl_call_siteQuant_19_export)
kldata <- rbind(kl_call_siteQuant_19_export, kl_call_siteQuant_22_export2)

library(homologene)

mouse_kinases <- human2mouse(kldata$motif.geneName)

kldata <- merge(kldata, mouse_kinases, by.x = "motif.geneName", by.y = "humanGene")
kldata$motif.geneName <- kldata$mouseGene

interactions <- unique(kldata[,c(21, 23)])
interactions <- interactions[!apply(interactions, 1, function(x) any(x=="")),] 
interactions <- drop_na(interactions)

#START HERE
#kldata: kinase is mouseGene, substrate is gene

#add substrate -> kinase in layer 2 where substrate is already in "from" column
layer2_from <- unique(DB_Layer2_mouse_filtered$from)
paths <- unique(DB_Layer2_mouse_filtered$path)
interactions_to_add <- interactions[interactions$gene %in% layer2_from,]
new <- data.frame(from = interactions$gene, to = interactions$mouseGene, source = "KL")
new$path <- paste0(new$from, "*", new$to)
new <- new[is.na(match(new$path, paths)),]
DB_Layer2_mouse_filtered <- rbind(DB_Layer2_mouse_filtered, new)

#add kinase -> substrate in layer 2 where substrate is already in "to" column
layer2_to <- unique(DB_Layer2_mouse_filtered$to)
paths <- unique(DB_Layer2_mouse_filtered$path)
interactions_to_add <- interactions[interactions$gene %in% layer2_to,]
new <- data.frame(from = interactions$mouseGene, to = interactions$gene, source = "KL")
new$path <- paste0(new$from, "*", new$to)
new <- new[is.na(match(new$path, paths)),]
DB_Layer2_mouse_filtered <- rbind(DB_Layer2_mouse_filtered, new)



#add kinase -> substrate in layer 3 where substrate is already in "to" column
layer3_from <- unique(DB_Layer3_mouse_filtered$from)
paths <- unique(DB_Layer3_mouse_filtered$path)
interactions_to_add <- interactions[interactions$gene %in% layer3_from,]
new <- data.frame(from = interactions$gene, to = interactions$mouseGene, source = "KL")
new$path <- paste0(new$from, "*", new$to)
new <- new[is.na(match(new$path, paths)),]
DB_Layer3_mouse_filtered <- rbind(DB_Layer3_mouse_filtered, new)

#add kinase -> substrate in layer 3 where substrate is already in "to" column
layer3_to <- unique(DB_Layer3_mouse_filtered$to)
paths <- unique(DB_Layer3_mouse_filtered$path)
interactions_to_add <- interactions[interactions$gene %in% layer3_to,]
new <- data.frame(from = interactions$mouseGene, to = interactions$gene, source = "KL")
new$path <- paste0(new$from, "*", new$to)
new <- new[is.na(match(new$path, paths)),]
DB_Layer3_mouse_filtered <- rbind(DB_Layer3_mouse_filtered, new)


#New db

DB.M <- list()
DB.M[[1]] <- DB_Layer1_mouse_filtered
DB.M[[2]] <- DB_Layer2_mouse_filtered
DB.M[[3]] <- DB_Layer3_mouse_filtered



