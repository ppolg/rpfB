# GO analysis 


########
# Func #
########

theme_mycopore <- function(base_size=10) {
  (theme_foundation(base_size=base_size, base_family="Arial")
   + theme(plot.title = element_text(face = "bold",
                                     size = rel(2.8), hjust = 0.5),
           text = element_text(color = "black"),
           panel.background = element_rect(colour = NA, fill = "white"),
           plot.background = element_rect(colour = NA, fill = "white"),
           panel.border = element_blank(),
           axis.title = element_text(face = "bold",size = rel(1.4)),
           axis.title.y = element_text(angle=90,vjust = 3, size = rel(1)),
           axis.title.x = element_text(vjust = -0.1, size = rel(1)),
           axis.text = element_text(face="bold",size = rel(1.3)), 
           axis.line = element_line(colour="black",size=1),
           axis.ticks = element_line(colour="black",size = 1),
           axis.ticks.length = unit(.35, "cm"),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           legend.key = element_rect(color = NA, fill = "white"),
           legend.position = "top",
           legend.background= element_rect(color = NA, fill = "white"),
           legend.direction = "horizontal",
           legend.key.size= unit(0.2, "cm"),
           legend.spacing.x = unit(0.3, "cm"),
           legend.text = element_text(color = "black", size = rel(1.2)),
           legend.title = element_text(face="italic", color = "black"),
           plot.margin=unit(c(12,6,6,6),"mm"),
           strip.background=element_rect(colour="grey90",fill="grey70"),
           strip.text = element_text(face="bold")
   ))
}


###########################################################################

########
# Init #
########

invisible(library(here))
invisible(library(DESeq2))
invisible(source(here("seq/R/mycopore_redux/mycopore_init.R")))


########
# Main #
########

# list for refseq/genbank numbers for smeg (any will do?)
msm_ncbi_list <- c("GCA_000767605.1",
                   "GCA_000015005.1",
                   "GCA_000283295.1",
                   "GCA_004307055.1",
                   "GCF_000767605.1",
                   "GCF_000015005.1",
                   "GCF_000283295.1",
                   "GCF_004307055.1")

# get RNAseq LFC
LFC_table <- read.csv(paste0(here(),"/Brindhaseqcounts/lfc_table_with_products.csv"))

# Get gffs (mycobrowser and NCBI, not sure which one works better)
gff <- get_gff_input("msmeg") %>%
  dplyr::mutate(locus_name = str_split_fixed(str_split_fixed(attributes, ";product=|;identified=|;this=", 2)[,1], "Name=", 2)[,2],) %>%
  dplyr::mutate(product = str_split_fixed(str_split_fixed(attributes, ";uniProt_AC", 2)[,1], "product=", 2)[,2],)

gff2 <- get_gff_input("msmeg3") %>%
  filter(type == "gene") %>%
  dplyr::mutate(locus_name = str_split_fixed(str_split_fixed(attributes, ";gbkey=", 2)[,1], "Name=", 2)[,2],) %>%
  dplyr::mutate(locus_name2 = sub(".*locus_tag=", "\\1",attributes)) %>%
  mutate(weird_locus = str_extract(attributes,"MSMEG_RS\\d{5}"))


gff2$realID <- apply(gff2[10:11],1, function(x) x[which.min(nchar(x))])

# Get full COG tables
COG_full <- read.csv(paste0(here(),"/COG/cog-20.cog.csv"), header = F) %>%
  dplyr::rename("Gene_ID" = 1,
                "NCBI_ID" = 2,
                "Protein_ID" = 3,
                "length" = 4,
                "footprint_coord" = 5,
                "footprint_length" = 6,
                "COG_ID" = 7,
                "reserved" = 8,
                "COG_class" = 9,
                "PSI_BLAST_bit" = 10,
                "PSI_BLAST_evalue" = 11,
                "COG_profile_length" = 12,
                "protein_footprint_coord" = 13)

COG_msm <- COG_full %>%
  filter(NCBI_ID %in% msm_ncbi_list)


# Definition and func cat per COG
COG_def <- read.csv(paste0(here(),"/COG/cog-20.def.tab.txt"), header = F, sep = "\t") %>%
  dplyr::rename("COG_ID" = 1,
                "func_cat" = 2,
                "COG_name" = 3,
                "gene_ass" = 4,
                "pathway_ass" = 5,
                "PubMed_ID" = 6,
                "PDB_ID" = 7)

# COG function group definitions
COG_functions <- read.csv(paste0(here(),"/COG/fun-20.tab.txt"), header = F, sep = "\t") %>%
  dplyr::rename("func_letter" = 1,
                "RGB_hex" = 2,
                "desc" = 3)


# append gff2
gff_msm <- gff2
gff_msm$COG <- COG_msm$COG_ID[match(x = gff_msm$locus_name2, table=COG_msm$Gene_ID)]
gff_msm$COG_str <- COG_def$func_cat[match(x = gff_msm$COG, table=COG_def$COG_ID)]
gff_msm$COG_main <- substring(gff_msm$COG_str,1,1)
gff_msm$COG_main_desc <- COG_functions$desc[match(x = gff_msm$COG_main, table=COG_functions$func_letter)]

gff_msm <- gff_msm %>%
  mutate(COG_main_desc = ifelse(is.na(COG_main_desc),"Not in COG database",COG_main_desc))

msm_sum <- gff_msm %>%
  group_by(COG_main_desc) %>%
  summarise(n = length(COG_main_desc))

gff_hits <- gff_msm %>%
  filter(realID %in% LFC_table$rowname)

hits_sum <- gff_hits %>%
  group_by(COG_main_desc) %>%
  summarise(n_LFC = length(COG_main_desc)) %>%
  add_row(COG_main_desc = "Cell motility",n_LFC= 0) %>%
  add_row(COG_main_desc = "Intracellular trafficking, secretion, and vesicular transport",n_LFC= 0) %>%
  add_row(COG_main_desc = "RNA processing and modification",n_LFC= 0)

msm_sum <- merge(msm_sum,hits_sum)

msm_sum <- msm_sum %>%
  mutate(k = n/sum(n),
         k_LFC = n_LFC/sum(n_LFC))

msm_melt <- msm_sum %>%
  select(1,4,5) %>%
  reshape2::melt()

# add COG to LFC table from DESeq2
LFC_table$category <- gff_msm$COG_main_desc[match(x = LFC_table$rowname, table=gff_msm$realID)]

# Separate up and down-reg in RNAseq
LFC_up <- LFC_table %>%
  filter(log2FoldChange > 0)

LFC_down <- LFC_table %>%
  filter(log2FoldChange < 0)

gff_hits_up <- gff_hits %>%
  filter(realID %in% LFC_up$rowname)

gff_hits_down <- gff_hits %>%
  filter(realID %in% LFC_down$rowname)

# sum per up/down
up_sum <- gff_hits_up %>%
  group_by(COG_main_desc) %>%
  summarise(n_LFC_up = length(COG_main_desc)) %>%
  add_row(COG_main_desc = "Cell motility",n_LFC_up= 0) %>%
  add_row(COG_main_desc = "Intracellular trafficking, secretion, and vesicular transport",n_LFC_up= 0) %>%
  add_row(COG_main_desc = "RNA processing and modification",n_LFC_up= 0)

down_sum <- gff_hits_down %>%
  group_by(COG_main_desc) %>%
  summarise(n_LFC_down = length(COG_main_desc)) %>%
  add_row(COG_main_desc = "Cell motility",n_LFC_down= 0) %>%
  add_row(COG_main_desc = "Intracellular trafficking, secretion, and vesicular transport",n_LFC_down= 0) %>%
  add_row(COG_main_desc = "RNA processing and modification",n_LFC_down= 0) %>%
  add_row(COG_main_desc = "Function unknown",n_LFC_down= 0) %>%
  add_row(COG_main_desc = "Mobilome: prophages, transposons",n_LFC_down= 0) 


updown_sum <- merge(merge(msm_sum, up_sum), down_sum) %>%
  select(-n_LFC,-k_LFC) %>%
  mutate(k_up = n_LFC_up/sum(n_LFC_up),
         k_down = n_LFC_down/sum(n_LFC_down))

updown_melt <- updown_sum %>%
  select(1,3,6,7) %>%
  reshape2::melt() %>%
  filter(COG_main_desc != "Not in COG database")

# for p-values, fold-change
updown_sum <- updown_sum %>%
  mutate(l2fc_up = log2(k_up/k),
         l2fc_down = log2(k_down/k),
         p_up = phyper(n_LFC_up,n,sum(n) - sum(n_LFC_up),sum(n_LFC_up), lower.tail=F),
         p_down = phyper(n_LFC_down,n,sum(n) - sum(n_LFC_down),sum(n_LFC_down), lower.tail=F)) %>%
  mutate(padj_up = p.adjust(p_up, method="BH"),
         padj_down = p.adjust(p_down, method="BH"))

l2fc_melt <- updown_sum %>%
  select(1,8,9) %>%
  reshape2::melt() %>%
  mutate(value = ifelse(value == -Inf,0,value))

# Figures
fig_COG_ratio <- ggplot(msm_melt, aes(x=COG_main_desc,y=value, colour=variable,fill=variable)) +
  geom_bar(stat= "identity", position = "dodge") +
  theme_mycopore() +
  xlab("") +
  ylab("% of all ORFs") +
  ggtitle("Percentage of ORFs in cat") +
  scale_fill_manual(values = c(colours$red,colours$purple)) +
  scale_colour_manual(values = c(colours$red,colours$purple)) +
  coord_flip()

print(fig_COG_ratio)

fig_COG_ratio_updown <- ggplot(updown_melt, aes(x=COG_main_desc,y=value, colour=variable,fill=variable)) +
  geom_bar(stat= "identity", position = "dodge") +
  theme_mycopore() +
  xlab("") +
  ylab("% of all ORFs") +
#  ggtitle("Percentage of ORFs in category") +
  scale_fill_manual(values = c(colours$red,colours$purple, colours$azure), labels = c("all","up","down")) +
  scale_colour_manual(values = c(colours$red,colours$purple, colours$azure), labels = c("all","up","down")) +
  scale_y_continuous(expand=c(0,0)) +
  coord_flip() +
  theme(legend.title = element_blank())

print(fig_COG_ratio_updown)


fig_COG_l2fc <- ggplot(l2fc_melt, aes(x=COG_main_desc,y=value, colour=variable,fill=variable)) +
  geom_bar(stat= "identity", position = "dodge") +
  theme_mycopore() +
  xlab("") +
  ylab("log2-fold enrichment") +
  ggtitle("COG enrichment log2-fold") +
  scale_fill_manual(values = c(colours$purple, colours$azure)) +
  scale_colour_manual(values = c(colours$purple, colours$azure)) +
  scale_y_continuous(expand=c(0,0)) +
  coord_flip() +
  theme(legend.title = element_blank())

print(fig_COG_l2fc)

fig_volcano_up <- ggplot(updown_sum, 
                         aes(x = l2fc_up, y = -log10(padj_up),
                             color=ifelse(padj_up<0.05,"grey60",colours$alex_R1),
                             fill=ifelse(padj_up<0.05,"grey60",colours$alex_R1))) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour=colours$alex_R1) +
  theme_mycopore() +
  scale_color_manual(values = c("grey60",colours$alex_R3)) +
  theme(legend.position = "none") +
  xlab("log2-fold enrichment") +
  ylab("-log10 p(adjusted)") +
  ggtitle(paste0("volcano - enrichment for up"))
  

print(fig_volcano_up)

fig_volcano_down <- ggplot(updown_sum, 
                         aes(x = l2fc_down, y = -log10(padj_down),
                             color=ifelse(padj_down<0.05,"grey60",colours$alex_R1),
                             fill=ifelse(padj_down<0.05,"grey60",colours$alex_R1))) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour=colours$alex_R1) +
  theme_mycopore() +
  scale_color_manual(values = c("grey60",colours$alex_R3)) +
  theme(legend.position = "none") +
  xlab("log2-fold enrichment") +
  ylab("-log10 p(adjusted)") +
  ggtitle(paste0("volcano - enrichment for down"))


print(fig_volcano_down)



# Write
fwrite(LFC_table, paste0(here::here(),"/LFC_categories.csv"))
fwrite(updown_sum, paste0(here::here(),"/updown.csv"))
