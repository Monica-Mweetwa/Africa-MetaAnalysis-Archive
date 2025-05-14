

#setwd("~/Dropbox/Africa-MetaAnalysis")
setwd("C:/Users/Monica/Dropbox/Africa-MetaAnalysis")

pacman::p_load(
  taxonomizr,      # To convert taxid IDs to taxonomy
  cluster,         # For Hierarchical Clustering
  rio,             # import/export
  here,            # relative file paths
  tidyverse,       # general data management and visualization
  tibble,          # general data management
  ape,             # to import and export phylogenetic files
  ggtree,          # to visualize phylogenetic files
  treeio,          # to visualize phylogenetic files
  ggnewscale)

library(readxl)
TXID <- read_excel("CaCo_Taxids.xlsx", sheet = "taxid")
taxid <- c(TXID$taxid)
taxa_lineage <- getTaxonomy(taxid, sqlFile = "nameNode.sqlite")
#Export taxa_lineage for 2nd tree
taxa_lineage2 <- as.data.frame(taxa_lineage) %>%
  rownames_to_column()%>%
  mutate(Nodename = apply(., 1, function(x) tail(na.omit(x), 1))) ##To get the names from each taxonomic rank

write.csv(taxa_lineage2, "CaColineage.csv")

library(readxl)
review_tree_pre <- as_tibble(read_csv("CaColineage_edited.csv"))
review_tree_pre <- review_tree_pre %>%
  column_to_rownames(var = "rowname")
review_tree <- makeNewick(review_tree_pre, quote = "'", excludeTerminalNAs=TRUE)
myTree1 <- ape::read.tree(text = review_tree)
myTree10 <- myTree1
myTree10$tip.label <- str_sub(myTree10$tip.label, end=-2)
myTree10$tip.label <- str_sub(myTree10$tip.label, start =2)

myTree10$node.label <- str_sub(myTree10$node.label, end=-2)
myTree10$node.label <- str_sub(myTree10$node.label, start =2)

#Add Metadata
add <- as_tibble(myTree10)
write.csv(add, "Phytree_tabs_CaCo/Node_labels_CaCo.csv")

#Import metadata for tree plot
node_lab <- read_csv("Phytree_tabs_CaCo/Node_labels_CaCo_edited.csv")

node_lab <- node_lab %>%
  mutate(Name = paste0(label)) %>%
  mutate(Length = ifelse(Rank == "Species", 0.7,
                         ifelse(Rank == "Genus", 0.6,
                                ifelse(Rank == "Family", 0.5,
                                       ifelse(Rank == "Order", 0.4,
                                              ifelse(Rank == "Class", 0.3,
                                                     ifelse(Rank == "Phylum", 0.2,
                                                            ifelse(Rank == "Superkingdom", 0.1,0)
                                                            )
                                                     )
                                              )
                                       )
                                )
                         )) %>%
  mutate(Length = as.numeric(Length))

## Change the level of Rank
levels(as.factor(node_lab$Rank))
node_lab$Rank <- factor(node_lab$Rank, levels = c("Superkingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "NA"))
levels(node_lab$Rank)


#Join metadata with phylotree
y33<- left_join(add, node_lab, by = 'label')

y33$label <- ifelse(y33$node == 199 , "Root", y33$label)
y34 <- as.treedata(y33)
y34

library(viridis)
okabe <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
p1 <- ggtree(y34, layout='circular', branch.length = 'Length') +
  geom_tree(aes(color = as.numeric(Healthy)),  size=3, inherit.aes = TRUE) +
  scale_colour_gradient(low = "black", high = "orange", name = "Number of\nstudies associated") +
  geom_tiplab2(aes(label= ifelse(Healthy==0, NA,Name2)), color = "black", align = F, size=8, hjust= -0.1) +
  new_scale_color() + 
  geom_nodepoint(aes(color = Rank), size=4, show.legend = TRUE) +
  geom_tippoint(aes(color = Rank), size=4) + 
  scale_colour_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "red"), na.value = NA, name = "Rank", na.translate=FALSE) +
  theme(#legend.position = c(1.2, 0.5),
    legend.position = "bottom",
    legend.text = element_text(size = 18),
    legend.title = element_text(size=24, face="bold"),
    legend.key = element_blank(),
    legend.background = element_blank(),
    legend.box.background = element_rect(colour = "black", linewidth = 1),
    plot.title = element_text(size = 24, hjust = -0.2),
    panel.background = element_blank(),
    plot.background = element_blank(), 
    panel.grid = element_blank()) +
  geom_rootpoint(colour = "red", size=3, shape=8)
library(Cairo)
Cairo(file = "Figures/phytree_CaCo_health2.png", 
      type = "png", 
      units = "px",
      width = 9000, 
      height = 8500, dpi = 350, bg = "white")
p1 +  theme(legend.position = "none")
invisible(dev.off())


p2 <- ggtree(y34, layout='circular', branch.length = 'Length') +
  geom_tree(aes(color = as.numeric(Diseased)),  size=3, inherit.aes = TRUE) +
  scale_colour_gradient(low = "black", high = "orange", name = "Number of\nstudies associated") +
  geom_tiplab2(aes(label= ifelse(Diseased==0, NA,Name2)), color = "black", align = F, size=8, hjust= -0.1) +
  new_scale_color() + 
  geom_nodepoint(aes(color = Rank), size=4, show.legend = TRUE) +
  geom_tippoint(aes(color = Rank), size=4) + 
  scale_colour_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "red"), na.value = NA, name = "Rank", na.translate=FALSE) +
  theme(#legend.position = c(1.2, 0.5),
    legend.position = "bottom",
    legend.text = element_text(size = 18),
    legend.title = element_text(size=24, face="bold"),
    legend.key = element_blank(),
    legend.background = element_blank(),
    legend.box.background = element_rect(colour = "black", linewidth = 1),
    plot.title = element_text(size = 24, hjust = -0.2),
    panel.background = element_blank(),
    plot.background = element_blank(), 
    panel.grid = element_blank()) +
  geom_rootpoint(colour = "red", size=3, shape=8)

Cairo(file = "phytree_CaCo_dis2.png", 
      type = "png", 
      units = "px",
      width = 9000, 
      height = 8500, dpi = 350, bg = "white")
p2 +  theme(legend.position = "none")
invisible(dev.off())

#Get legend
Cairo(file = "phytree_CaCo_legend.png", 
      type = "png", 
      units = "px",
      width = 9000, 
      height = 8500, dpi = 350, bg = "white")
p2 +  theme(legend.position = "right")
invisible(dev.off())

p3 <- p2 +  theme(legend.position = "right")
grobs <- ggplotGrob(p3)$grobs
legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]

#Save legend
Cairo(file = "phytree_CaCo_legend.png", 
      type = "png", 
      units = "px",
      width = 7000, 
      height = 7000, dpi = 600, bg = "transparent")
b1 <- plot(legend)
b1 + theme(legend.background = element_blank(),
           legend.box.background = element_rect(colour = "black", linewidth = 1),
           plot.title = element_text(size = 24, hjust = -0.2),
           panel.background = element_blank(),
           plot.background = element_blank(), 
           panel.grid = element_blank())
invisible(dev.off())
