library(readr)
library(dplyr)
library(SVbyEye)
library(pafr)
library(stringr)
library(GenomicRanges)

##### ------------------------------------
##### Set working directory
##### ------------------------------------

setwd("/Users/gallegosda/Library/CloudStorage/OneDrive-NationalInstitutesofHealth/2-14-26_AncientKZFPclusterFigure_clean")

##### ------------------------------------
##### Load all six PAF files
##### ------------------------------------

# Target: Human | Query: Dog
dog_human_paf <- readPaf(
  paf.file = "./paf/QUERY_Dog_canFam3_chr16_TARGET_Human_hg38_chr7_ancientKZFPcluster_lastz_trimmed.paf",
  include.paf.tags = TRUE, restrict.paf.tags = "cg"
)
dog_human_paf$q.name <- "Dog_chr16"
dog_human_paf$t.name <- "Human_chr7"

# Target: Dog | Query: Mouse
mouse_dog_paf <- readPaf(
  paf.file = "./paf/QUERY_Mouse_mm10_chr6_TARGET_Dog_canFam3_chr16_ancientKZFPcluster_lastz_trimmed.paf",
  include.paf.tags = TRUE, restrict.paf.tags = "cg"
)
mouse_dog_paf$q.name <- "Mouse_chr6"
mouse_dog_paf$t.name <- "Dog_chr16"

# Target: Mouse | Query: Opossum
opossum_mouse_paf <- readPaf(
  paf.file = "./paf/QUERY_Opossum_monDom5_chr8_TARGET_Mouse_mm10_chr6_ancientKZFPcluster_lastz_trimmed.paf",
  include.paf.tags = TRUE, restrict.paf.tags = "cg"
)
opossum_mouse_paf$q.name <- "Opossum_chr8"
opossum_mouse_paf$t.name <- "Mouse_chr6"

# Target: Opossum | Query: Platypus
platypus_opossum_paf <- readPaf(
  paf.file = "./paf/QUERY_Platypus_ornAna1_Contig1664_TARGET_Opossum_monDom5_chr8_ancientKZFPcluster_lastz_trimmed.paf",
  include.paf.tags = TRUE, restrict.paf.tags = "cg"
)
platypus_opossum_paf$q.name <- "Platypus_contig1664"
platypus_opossum_paf$t.name <- "Opossum_chr8"

# Target: Platypus | Query: Chicken
chicken_platypus_paf <- readPaf(
  paf.file = "./paf/QUERY_Chicken_galGal4_chr2_TARGET_Platypus_ornAna1_Contig1664_ancientKZFPcluster_lastz_trimmed.paf", # Chicken chr 2, Platypus Contig1664
  include.paf.tags = TRUE, restrict.paf.tags = "cg"
)
chicken_platypus_paf$q.name <- "Chicken_chr2"
chicken_platypus_paf$t.name <- "Platypus_contig1664"

# Target: Chicken | Query: Chinese Softshell Turtle [Pelodiscus sinensis]
turtle_chicken_paf <- readPaf(
  paf.file = "./paf/QUERY_Turtle_PelSin1_JH212076.1_TARGET_Chicken_galGal4_chr2_ancientKZFPcluster_lastz_trimmed.paf", # Chicken chr 2, Turtle JH212076.1
  include.paf.tags = TRUE, restrict.paf.tags = "cg"
)
turtle_chicken_paf$q.name <- "Turtle_JH212076_1"
turtle_chicken_paf$t.name <- "Chicken_chr2"

##### ------------------------------------
##### Load all seven custom GTF files to specify where to place rectangles for each KZFP
##### ------------------------------------

# Human hg38 annotation
humanHg38_anno <- read_tsv("./gtf/hg38_kzfpClusterFigureAnnotation_labeledGenes.gtf",
                           comment = "#", col_names = FALSE)

# Dog canFam3 annotation
dogCanFam3_anno <- read_tsv("./gtf/canFam3_kzfpCluster.gtf",
                            comment = "#", col_names = FALSE)

# Mouse mm10 annotation
mouseMm10_anno <- read_tsv("./gtf/mm10_kzfpCluster.gtf",
                           comment = "#", col_names = FALSE)

# Opossum monDom5 annotation
opossumMonDom5_anno <- read_tsv("./gtf/monDom5_kzfpCluster_Liu_Stubbs_2014.gtf",
                                comment = "#", col_names = FALSE)

# Platypus ornAna1 annotation
platypusOrnAna1_anno <- read_tsv("./gtf/ornAna1_kzfpCluster_Imbeault_2017.gtf",
                                 comment = "#", col_names = FALSE)

# Chicken galGal4 annotation
chickenGalGal4_anno <- read_tsv("./gtf/galGal4_kzfpCluster_Imbeault_2017.gtf",
                                comment = "#", col_names = FALSE)

# Turtle PelSin_1 annotation
turtlePelSin_1_anno <- read_tsv("./gtf/Pelodiscus_sinensis_kzfpCluster_Imbeault_2017.gtf",
                                comment = "#", col_names = FALSE)

##### ------------------------------------
##### Format annotations for plot
##### ------------------------------------
 
# Human Annotation
humanHg38_chr7_anno <- humanHg38_anno %>%
  dplyr::mutate(gene_name = str_extract(X9, '(?<=gene_name ")[^"]+')) %>%
  dplyr::filter(., X1=='chr7'
                & X4>149059641
                & X5<149507809 ) %>%
  mutate(., chr='Human_chr7')%>%
  mutate(., start=X4-149059641) %>%
  mutate(., end=X5-149059641) %>%
  select(., chr, start, end, X7, gene_name) %>%
  arrange(., start) %>%
  dplyr::rename(., strand = X7)

# Dog Annotation
dogCanFam3_chr16_anno <- dogCanFam3_anno %>%
  dplyr::mutate(gene_name = str_extract(X9, '(?<=gene_name ")[^"]+')) %>%
  dplyr::filter(., X1=='chr16'
                & X4>14134497
                & X5<14436586 ) %>%
  mutate(., chr='Dog_chr16')%>%
  mutate(., start=X4-14134497) %>%
  mutate(., end=X5-14134497) %>%
  select(., chr, start, end, X7, gene_name) %>%
  arrange(., start) %>%
  dplyr::rename(., strand = X7)

# Mouse Annotation
mouseMm10_chr6_anno <- mouseMm10_anno %>%
  dplyr::mutate(gene_name = str_extract(X9, '(?<=gene_name ")[^"]+')) %>%
  dplyr::filter(., X1=='chr6'
                & X4>47809266
                & X5<48096593 ) %>%
  mutate(., chr='Mouse_chr6')%>%
  mutate(., start=X4-47809266) %>%
  mutate(., end=X5-47809266) %>%
  select(., chr, start, end, X7, gene_name) %>%
  arrange(., start) %>%
  dplyr::rename(., strand = X7)

# Opossum Annotation
opossumMonDom5_chr8_anno <- opossumMonDom5_anno %>%
  dplyr::mutate(gene_name = str_extract(X9, '(?<=gene_name ")[^"]+')) %>%
  dplyr::filter(., X1=='chr8'
                & X4>213448780
                & X5<213953073 ) %>%
  mutate(., chr='Opossum_chr8')%>%
  mutate(., start=X4-213448780) %>%
  mutate(., end=X5-213448780) %>%
  select(., chr, start, end, X7, gene_name) %>%
  arrange(., start) %>%
  dplyr::rename(., strand = X7)

# Platypus
platypusOrnAna1_contig1664_anno <- platypusOrnAna1_anno %>%
  dplyr::mutate(gene_name = str_extract(X9, '(?<=gene_name ")[^"]+')) %>%
  dplyr::filter(., X1=='Contig1664'
                & X4>4339
                & X5<68062 ) %>%
  mutate(., chr='Platypus_contig1664')%>%
  mutate(., start=X4-4339) %>%
  mutate(., end=X5-4339) %>%
  select(., chr, start, end, X7, gene_name) %>%
  arrange(., start) %>%
  dplyr::rename(., strand = X7)

# Chicken
chickenGalGal4_chr2_anno <- chickenGalGal4_anno %>%
  dplyr::mutate(gene_name = str_extract(X9, '(?<=gene_name ")[^"]+')) %>%
  dplyr::filter(., X1=='chr2'
                & X4>479515
                & X5<565923 ) %>%
  mutate(., chr='Chicken_chr2')%>%
  mutate(., start=X4-479515) %>%
  mutate(., end=X5-479515) %>%
  select(., chr, start, end, X7, gene_name) %>%
  arrange(., start) %>%
  dplyr::rename(., strand = X7)

# Turtle
turtlePelSin_1_JH212076_1_anno <- turtlePelSin_1_anno %>%
  dplyr::mutate(gene_name = str_extract(X9, '(?<=gene_id ")[^"]+')) %>%
  dplyr::filter(., X1=='JH212076.1'
                & X4>257571
                & X5<434685 ) %>%
  mutate(., chr='Turtle_JH212076_1')%>%
  mutate(., start=X4-257571) %>%
  mutate(., end=X5-257571) %>%
  select(., chr, start, end, X7, gene_name) %>%
  arrange(., start) %>%
  dplyr::rename(., strand = X7)

##### ------------------------------------
##### Bind all dataframes in order
##### ------------------------------------

ancientKRABcluster_list <- list(dog_human_paf = dog_human_paf, mouse_dog_paf = mouse_dog_paf, opossum_mouse_paf = opossum_mouse_paf, platypus_opossum_paf = platypus_opossum_paf, chicken_platypus_paf = chicken_platypus_paf, turtle_chicken_paf = turtle_chicken_paf)
allAncientKRABclusterPafs <- bind_rows(ancientKRABcluster_list)

all_ancientKRABcluster_anno <- dplyr::bind_rows(humanHg38_chr7_anno, dogCanFam3_chr16_anno, mouseMm10_chr6_anno, opossumMonDom5_chr8_anno, platypusOrnAna1_contig1664_anno, chickenGalGal4_chr2_anno, turtlePelSin_1_JH212076_1_anno)

# Retain gene_name metadata
all_ancientKRABcluster_anno_gr <- GenomicRanges::makeGRangesFromDataFrame(
  all_ancientKRABcluster_anno, 
  seqnames.field = "chr",
  start.field    = "start",
  end.field      = "end",
  strand.field   = "strand",
  keep.extra.columns = TRUE   # <- key
)

# Assign the assembly order in the plot
seqnames.order <- c("Human_chr7", "Dog_chr16", "Mouse_chr6", "Opossum_chr8", "Platypus_contig1664", "Chicken_chr2", "Turtle_JH212076_1")

# Check that worked
head(allAncientKRABclusterPafs)

# Assemble base plot
all_ancientKRABcluster_dir <- plotAVA(
  paf.table = allAncientKRABclusterPafs,
  color.by = "direction",
  # color.by = "identity",
  color.palette = c("+" = "#F5DEB3", "-" = "#C4C4C2"),
  seqnames.order = seqnames.order
)

# Add rectangle KZFP gene annotations, and optionally gene_name labels
all_ancientKRABcluster_dir_anno <- addAnnotation(
  ggplot.obj = all_ancientKRABcluster_dir, 
  annot.gr = all_ancientKRABcluster_anno_gr, 
  shape = "rectangle",
  coordinate.space = 'self', 
  y.label.id = 'seqnames', 
  # label.by = 'gene_name',
  annotation.level = 0
)

# Check plot
all_ancientKRABcluster_dir_anno

# Inspect layers to adjust gene_name label size as needed
all_ancientKRABcluster_dir_anno$layers
# 
# all_ancientKRABcluster_dir_anno$layers[[4]]$aes_params$size <- 1.5
# all_ancientKRABcluster_dir_anno
# 
# all_ancientKRABcluster_dir_anno$layers[[4]]$mapping$y <-
#   rlang::expr(y.offset + 0)
# 
# all_ancientKRABcluster_dir_anno


##### ------------------------------------
##### Save plot
##### ------------------------------------

# Save pdf of plot
pdf("./plots/2-14-26_Ancient_KZFP_Cluster_Conservation_in_Vertebrates_trimmed_000.pdf", width = 10, height = 5)
print(all_ancientKRABcluster_dir_anno)
dev.off()

# Add title and other adjustments here
plotWithTitle <- all_ancientKRABcluster_dir_anno + labs(
  title = "Ancient KRAB ZFP Cluster Alignments in Vertebrates",
  x = "Ancient KRAB ZFP Cluster Locus Size (bp)",
  # y = "Sequence Blocks",
  color = "Strand Direction"
)

# Save with adjustments
pdf("./plots/2-14-26_Ancient_KZFP_Cluster_Conservation_in_Vertebrates_withTitle_000.pdf", width = 10, height = 5)
print(plotWithTitle)
dev.off()

##### ------------------------------------
##### View final plot
##### ------------------------------------

plotWithTitle