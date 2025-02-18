# Load required packages
pacman::p_load(tidyverse, yaml)

# Read 'listAll' genes
listAll <- read_yaml("listAll.yaml")

# Prepare data for venn diagrams - per comparison across lines
list_F2tF2d_abdomen <- list( TRA=c(listAll[["listTRAa"]][["TRAa_F2_3N_F2_2N_up"]],listAll[["listTRAa"]][["TRAa_F2_3N_F2_2N_down"]]),
                             WGDL=c(listAll[["listWGDLa"]][["WGDLa_F2_3N_F2_2N_up"]],listAll[["listWGDLa"]][["WGDLa_F2_3N_F2_2N_down"]]),
                             WOM=c(listAll[["listWOMa"]][["WOMa_F2_3N_F2_2N_up"]],listAll[["listWOMa"]][["WOMa_F2_3N_F2_2N_down"]])
                             )
list_F2tF2d_head <- list( TRA=c(listAll[["listTRAh"]][["TRAh_F2_3N_F2_2N_up"]],listAll[["listTRAh"]][["TRAh_F2_3N_F2_2N_down"]]),
                          WGDL=c(listAll[["listWGDLh"]][["WGDLh_F2_3N_F2_2N_up"]],listAll[["listWGDLh"]][["WGDLh_F2_3N_F2_2N_down"]]),
                          WOM=c(listAll[["listWOMh"]][["WOMh_F2_3N_F2_2N_up"]],listAll[["listWOMh"]][["WOMh_F2_3N_F2_2N_down"]])
                          )
list_F14tF14d_abdomen <- list( TRA=c(listAll[["listTRAa"]][["TRAa_F14_3N_F14_2N_up"]],listAll[["listTRAa"]][["TRAa_F14_3N_F14_2N_down"]]),
                               WGDL=c(listAll[["listWGDLa"]][["WGDLa_F14_3N_F14_2N_up"]],listAll[["listWGDLa"]][["WGDLa_F14_3N_F14_2N_down"]]),
                               WOM=c(listAll[["listWOMa"]][["WOMa_F14_3N_F14_2N_up"]],listAll[["listWOMa"]][["WOMa_F14_3N_F14_2N_down"]])
                               #,WPL=c(listAll[["listWPLa"]][["WPLa_F_3N_F_2N_up"]],listAll[["listWPLa"]][["WPLa_F_3N_F_2N_down"]])
                               )
list_F14tF14d_head <- list( TRA=c(listAll[["listTRAh"]][["TRAh_F14_3N_F14_2N_up"]],listAll[["listTRAh"]][["TRAh_F14_3N_F14_2N_down"]]),
                            WGDL=c(listAll[["listWGDLh"]][["WGDLh_F14_3N_F14_2N_up"]],listAll[["listWGDLh"]][["WGDLh_F14_3N_F14_2N_down"]]),
                            WOM=c(listAll[["listWOMh"]][["WOMh_F14_3N_F14_2N_up"]],listAll[["listWOMh"]][["WOMh_F14_3N_F14_2N_down"]])
                            #,WPL=c(listAll[["listWPLh"]][["WPLh_F_3N_F_2N_up"]],listAll[["listWPLh"]][["WPLh_F_3N_F_2N_down"]])
                            )
list_F2dF14d_abdomen <- list( TRA=c(listAll[["listTRAa"]][["TRAa_F2_2N_F14_2N_up"]],listAll[["listTRAa"]][["TRAa_F2_2N_F14_2N_down"]]),
                             WGDL=c(listAll[["listWGDLa"]][["WGDLa_F2_2N_F14_2N_up"]],listAll[["listWGDLa"]][["WGDLa_F2_2N_F14_2N_down"]]),
                             WOM=c(listAll[["listWOMa"]][["WOMa_F2_2N_F14_2N_up"]],listAll[["listWOMa"]][["WOMa_F2_2N_F14_2N_down"]])
                             )
list_F2dF14d_head <- list( TRA=c(listAll[["listTRAh"]][["TRAh_F2_2N_F14_2N_up"]],listAll[["listTRAh"]][["TRAh_F2_2N_F14_2N_down"]]),
                          WGDL=c(listAll[["listWGDLh"]][["WGDLh_F2_2N_F14_2N_up"]],listAll[["listWGDLh"]][["WGDLh_F2_2N_F14_2N_down"]]),
                          WOM=c(listAll[["listWOMh"]][["WOMh_F2_2N_F14_2N_up"]],listAll[["listWOMh"]][["WOMh_F2_2N_F14_2N_down"]])
                          )
list_F2tF14t_abdomen <- list( TRA=c(listAll[["listTRAa"]][["TRAa_F2_3N_F14_3N_up"]],listAll[["listTRAa"]][["TRAa_F2_3N_F14_3N_down"]]),
                              WGDL=c(listAll[["listWGDLa"]][["WGDLa_F2_3N_F14_3N_up"]],listAll[["listWGDLa"]][["WGDLa_F2_3N_F14_3N_down"]]),
                              WOM=c(listAll[["listWOMa"]][["WOMa_F2_3N_F14_3N_up"]],listAll[["listWOMa"]][["WOMa_F2_3N_F14_3N_down"]])
                              )
list_F2tF14t_head <- list( TRA=c(listAll[["listTRAh"]][["TRAh_F2_3N_F14_3N_up"]],listAll[["listTRAh"]][["TRAh_F2_3N_F14_3N_down"]]),
                           WGDL=c(listAll[["listWGDLh"]][["WGDLh_F2_3N_F14_3N_up"]],listAll[["listWGDLh"]][["WGDLh_F2_3N_F14_3N_down"]]),
                           WOM=c(listAll[["listWOMh"]][["WOMh_F2_3N_F14_3N_up"]],listAll[["listWOMh"]][["WOMh_F2_3N_F14_3N_down"]])
                           )
list_F3dF3h_abdomen <- list( TRA=c(listAll[["listTRAa"]][["TRAa_F3_2N_F3_1N_up"]],listAll[["listTRAa"]][["TRAa_F3_2N_F3_1N_down"]]),
                              WGDL=c(listAll[["listWGDLa"]][["WGDLa_F5_2N_F5_1N_up"]],listAll[["listWGDLa"]][["WGDLa_F5_2N_F5_1N_down"]]),
                              WOM=c(listAll[["listWOMa"]][["WOMa_F3_2N_F3_1N_up"]],listAll[["listWOMa"]][["WOMa_F3_2N_F3_1N_down"]])
                              )
list_F3dF3h_head <- list( TRA=c(listAll[["listTRAh"]][["TRAh_F3_2N_F3_1N_up"]],listAll[["listTRAh"]][["TRAh_F3_2N_F3_1N_down"]]),
                           WGDL=c(listAll[["listWGDLh"]][["WGDLh_F5_2N_F5_1N_up"]],listAll[["listWGDLh"]][["WGDLh_F5_2N_F5_1N_down"]]),
                           WOM=c(listAll[["listWOMh"]][["WOMh_F3_2N_F3_1N_up"]],listAll[["listWOMh"]][["WOMh_F3_2N_F3_1N_down"]])
                           )
list_F15dF15h_abdomen <- list( TRA=c(listAll[["listTRAa"]][["TRAa_F15_2N_F15_1N_up"]],listAll[["listTRAa"]][["TRAa_F15_2N_F15_1N_down"]]),
                             WGDL=c(listAll[["listWGDLa"]][["WGDLa_F15_2N_F15_1N_up"]],listAll[["listWGDLa"]][["WGDLa_F15_2N_F15_1N_down"]]),
                             WOM=c(listAll[["listWOMa"]][["WOMa_F15_2N_F15_1N_up"]],listAll[["listWOMa"]][["WOMa_F15_2N_F15_1N_down"]])
                             #,WPL=c(listAll[["listWPLa"]][["WPLa_M_2N_M_1N_up"]],listAll[["listWPLa"]][["WPLa_M_2N_M_1N_down"]])
                             )
list_F15dF15h_head <- list( TRA=c(listAll[["listTRAh"]][["TRAh_F15_2N_F15_1N_up"]],listAll[["listTRAh"]][["TRAh_F15_2N_F15_1N_down"]]),
                          WGDL=c(listAll[["listWGDLh"]][["WGDLh_F15_2N_F15_1N_up"]],listAll[["listWGDLh"]][["WGDLh_F15_2N_F15_1N_down"]]),
                          WOM=c(listAll[["listWOMh"]][["WOMh_F15_2N_F15_1N_up"]],listAll[["listWOMh"]][["WOMh_F15_2N_F15_1N_down"]])
                          #,WPL=c(listAll[["listWPLh"]][["WPLh_M_2N_M_1N_up"]],listAll[["listWPLh"]][["WPLh_M_2N_M_1N_down"]])
                          )
list_F3hF15h_abdomen <- list( TRA=c(listAll[["listTRAa"]][["TRAa_F3_1N_F15_1N_up"]],listAll[["listTRAa"]][["TRAa_F3_1N_F15_1N_down"]]),
                             WGDL=c(listAll[["listWGDLa"]][["WGDLa_F5_1N_F15_1N_up"]],listAll[["listWGDLa"]][["WGDLa_F5_1N_F15_1N_down"]]),
                             WOM=c(listAll[["listWOMa"]][["WOMa_F3_1N_F15_1N_up"]],listAll[["listWOMa"]][["WOMa_F3_1N_F15_1N_down"]])
                             )
list_F3hF15h_head <- list( TRA=c(listAll[["listTRAh"]][["TRAh_F3_1N_F15_1N_up"]],listAll[["listTRAh"]][["TRAh_F3_1N_F15_1N_down"]]),
                          WGDL=c(listAll[["listWGDLh"]][["WGDLh_F5_1N_F15_1N_up"]],listAll[["listWGDLh"]][["WGDLh_F5_1N_F15_1N_down"]]),
                          WOM=c(listAll[["listWOMh"]][["WOMh_F3_1N_F15_1N_up"]],listAll[["listWOMh"]][["WOMh_F3_1N_F15_1N_down"]])
                          )
list_F3dF15d_abdomen <- list( TRA=c(listAll[["listTRAa"]][["TRAa_F3_2N_F15_2N_up"]],listAll[["listTRAa"]][["TRAa_F3_2N_F15_2N_down"]]),
                              WGDL=c(listAll[["listWGDLa"]][["WGDLa_F5_2N_F15_2N_up"]],listAll[["listWGDLa"]][["WGDLa_F5_2N_F15_2N_down"]]),
                              WOM=c(listAll[["listWOMa"]][["WOMa_F3_2N_F15_2N_up"]],listAll[["listWOMa"]][["WOMa_F3_2N_F15_2N_down"]])
                              )
list_F3dF15d_head <- list( TRA=c(listAll[["listTRAh"]][["TRAh_F3_2N_F15_2N_up"]],listAll[["listTRAh"]][["TRAh_F3_2N_F15_2N_down"]]),
                           WGDL=c(listAll[["listWGDLh"]][["WGDLh_F5_2N_F15_2N_up"]],listAll[["listWGDLh"]][["WGDLh_F5_2N_F15_2N_down"]]),
                           WOM=c(listAll[["listWOMh"]][["WOMh_F3_2N_F15_2N_up"]],listAll[["listWOMh"]][["WOMh_F3_2N_F15_2N_down"]])
                           )

forVennDiagram <- c("list_F2tF2d_abdomen", "list_F2tF2d_head", "list_F14tF14d_abdomen", "list_F14tF14d_head", "list_F2dF14d_abdomen", 
                    "list_F2dF14d_head", "list_F2tF14t_abdomen", "list_F2tF14t_head", "list_F3dF3h_abdomen", "list_F3dF3h_head", 
                    "list_F15dF15h_abdomen", "list_F15dF15h_head", "list_F3hF15h_abdomen", "list_F3hF15h_head", "list_F3dF15d_abdomen", "list_F3dF15d_head")

# Venn diagram for gene lists
library(ggvenn)
sigGenesList <- list_F2dF14d_abdomen
pdf("VennDiagram_F15dF15h_head2.pdf", width = 4, height = 4)
ggvenn(sigGenesList, fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
       stroke_size = 0.5, set_name_size = 4, text_size = 2.5)
dev.off()

## Get the actual sets from lists for the GO analysis
library(gplots)
itemsAbdomen <- venn(list_F15dF15h_abdomen, show.plot = FALSE)
vennAbdomen <- attributes(itemsAbdomen)$intersections
lengths(vennAbdomen)
itemsHead <- venn(list_F15dF15h_head, show.plot = FALSE)
vennHead <- attributes(itemsHead)$intersections
lengths(vennHead)
# Write to file
write(as.yaml(vennAbdomen), file = "VennDataAbdomen_F15dF15h.yaml")
write(as.yaml(vennHead), file = "VennDataHead_F15dF15h.yaml")

F2tF14t <- list(a=vennAbdomen$WOM, b=vennAbdomen$WGDL, c=vennHead$WOM, d=vennHead$WGDL, e=vennHead$TRA)
outList <- F2tF14t

# GO enrichment
term2gene <- read.delim("TERM2GENE.txt", sep = "\t", col.names = c("GO_ID","Gene_ID"))
term2name <- read.delim("TERM2NAME.txt", sep = "\t", col.names = c("GO_ID","GO_Name"))
library(clusterProfiler)
clustersComparedOwn <- compareCluster(geneCluster = outList, fun = "enricher", TERM2GENE=term2gene, TERM2NAME = term2name, pvalueCutoff = 0.05)
clustersComparedOwn2 <- as.data.frame(clustersComparedOwn)
dotplot(clustersComparedOwn, showCategory = 10, label_format = 50, font.size=10)
cnetplot(clustersComparedOwn, showCategory=20)

# Prepare data for venn diagrams - additional ones comparing lines to find polyploidy specific genes
list_female <- list( F2tF2d=unique(unlist(flatten(list_F2tF2d_abdomen), flatten(list_F2tF2d_head))),
                     wpl=unique(c(list_F14tF14d_abdomen$WPL, list_F14tF14d_head$WPL)),
                     F14tF14d=unique(unlist(flatten(list_F14tF14d_abdomen[-4]), flatten(list_F14tF14d_head[-4]))))
list_male <- list( F3dF3h=unique(unlist(flatten(list_F3dF3h_abdomen), flatten(list_F3dF3h_head))),
                     wpl=unique(c(list_F15dF15h_abdomen$WPL, list_F15dF15h_head$WPL)),
                     F15dF15h=unique(unlist(flatten(list_F15dF15h_abdomen[-4]), flatten(list_F15dF15h_head[-4]))))
pdf("VennDiagram_Male.pdf", width = 4, height = 4)
ggvenn(list_male, fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
       stroke_size = 0.5, set_name_size = 4, text_size = 2.5)
dev.off()

# GO for the female and male new venn diagram data from the previous step
# Female
itemsFemale <- venn(list_female, show.plot = FALSE)
vennFemale <- attributes(itemsFemale)$intersections
lengths(vennFemale)
# Male
itemsMale <- venn(list_male, show.plot = FALSE)
vennMale <- attributes(itemsMale)$intersections
lengths(vennMale)
# Give it to outlist and run the script for GO above
outList <- vennMale
write(as.yaml(vennFemale), file = "VennDiagram_Female.yaml")
write(as.yaml(vennMale), file = "VennDiagram_Male.yaml")



