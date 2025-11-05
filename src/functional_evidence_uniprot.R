

#!/usr/bin/env Rscript
# Title: UniProt Evidence Plot for Single Variant Example
# Description: Standalone script to illustrate functional evidence from UniProt for a single hard-coded variant.
# Author: Switzerland Omics

filename <- "functional_evidence_uniprot"

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(rtracklayer)
  library(stringr)
  library(wesanderson)
  library(data.table)
  library(ggpubr)
  library(gridExtra)
  library(tidyr)
  library(ggrepel)
})

#-------------------------------------------------------------
# Hard-coded single variant example
#-------------------------------------------------------------

df_report <- data.frame(
  VARIANT_SCORE = 8,
  sample.id = "SAMPLE001",
  SYMBOL = "TNFAIP3",
  Protein_position = 307,
  CDS_position = 920,
  variant_type = "SNV",
  HGVS_p = "p.Leu307Ter",
  HGVS_c = "c.920T>A",
  chr = 6,
  pos = 137877190,
  ref = "T",
  alt = "A",
  consequence = "stop_gained",
  transcripts = "ENST00000612899.5 (MANE Select), ENST00000237289.8, ENST00000485192.1",
  consequence_summary = "pLoF High-confidence",
  gnomAD_version = "v4.1.0",
  gnomAD_AC = 1,
  gnomAD_AN = 1613922,
  gnomAD_AF = 6.2e-7,
  gnomAD_hom = 0,
  gnomAD_filter = "PASS",
  Uniprot_acc = "P21580",
  Uniprot_entry = "TNAP3_HUMAN",
  stringsAsFactors = FALSE
)

#-------------------------------------------------------------
# Automatically generate and wrap variant_info
#-------------------------------------------------------------
variant_info <- df_report %>%
  mutate(
    variant_str = paste0(chr, "-", pos, " ", ref, ">", alt, " (GRCh38); ", variant_type, "."),
    gene_str = paste0("Gene: ", SYMBOL, "; Transcripts affected: ", transcripts, "."),
    consequence_str = paste0("Consequence: ", consequence, " (", HGVS_p, ", ", HGVS_c, "); ", consequence_summary, "."),
    gnomad_str = paste0(
      "Reference population gnomAD ", gnomAD_version, ": AC=", gnomAD_AC, ", AN=", gnomAD_AN,
      ", AF=", format(gnomAD_AF, scientific = TRUE),
      "; Homozygotes=", gnomAD_hom, "; Filter=", gnomAD_filter, "."
    ),
    subtitle = paste(variant_str, gene_str, consequence_str, gnomad_str, sep = "\n")
  ) %>%
  pull(subtitle) %>%
  str_wrap(width = 80)

#-------------------------------------------------------------
# Import UniProt data
#-------------------------------------------------------------
df_uniprot <- readGFF("../data/uniprot/uniprot-filtered-organism%3A%22Homo+sapiens+%28Human%29+%5B9606%5D%22.gff")
df_uniprot_meta <- read.csv("../data/uniprot/uniprot-filtered-organism%3A%22Homo+sapiens+%28Human%29+%5B9606%5D%22.tab", sep = "\t")

# rename 'Entry' column to match expected field name
colnames(df_uniprot_meta)[colnames(df_uniprot_meta) == 'Entry'] <- 'seqid'

# handle possible variation in column name for Gene.names
if (!"Gene.names" %in% names(df_uniprot_meta)) {
  alt_col <- grep("^Gene\\.names", names(df_uniprot_meta), value = TRUE)[1]
  names(df_uniprot_meta)[names(df_uniprot_meta) == alt_col] <- "Gene.names"
}

# tidy metadata table safely
df_uniprot_meta_tidy <- df_uniprot_meta %>%
  tidyr::separate_rows(all_of("Gene.names"), sep = " ") %>%
  dplyr::rename(SYMBOL = all_of("Gene.names"))

#-------------------------------------------------------------
# Match UniProt data to target gene
#-------------------------------------------------------------
target_gene <- df_report$SYMBOL[1]

seqids_for_target <- df_uniprot_meta_tidy %>%
  filter(SYMBOL == target_gene) %>%
  pull(seqid)

df_report$seqid <- seqids_for_target[1]

dt_uniprot_filtered <- as.data.table(df_uniprot) %>%
  filter(seqid %in% seqids_for_target)

df_report_position <- df_report %>%
  select(SYMBOL, seqid, Protein_position, CDS_position, VARIANT_SCORE)

filtered_dt_uniprot <- merge(
  df_report %>% select(SYMBOL, seqid) %>% distinct(),
  dt_uniprot_filtered
)

#-------------------------------------------------------------
# Colours
#-------------------------------------------------------------
wes_pal <- rep(c(
  "#00A08A", "#F2AD00", "#F98400", "#5BBCD6",
  "#ECCBAE", "#046C9A", "#D69C4E", "#ABDDDE",
  "#DD8D29", "#E2D200", "#46ACC8", "#E58601", "#B40F20",
  "#F1BB7B", "#FD6467", "#5B1A18", "#D67236",
  "#E6A0C4", "#C6CDF7", "#D8A499", "#7294D4",
  "#446455", "#FDD262", "#D3DDDC", "#C7B19C"
), 7)

#-------------------------------------------------------------
# Plot function with non-overlapping annotations
#-------------------------------------------------------------


create_vlines <- function(data) {
  geom_vline(
    data = data,
    aes(xintercept = Protein_position),
    colour = "black", linewidth = 1
  )
}

create_plot <- function() {
  x <- filtered_dt_uniprot
  pos <- df_report_position
  symbol <- unique(x$SYMBOL)
  
  x$position_label <- rowMeans(x[, c("start", "end")], na.rm = TRUE)
  x$end[x$end == x$start] <- x$end + 1
  
  Domain <- c("Chain", "Domain", "Region", "Motif")
  Structure <- c("Helix", "Turn", "Beta strand")
  
  # classify before filtering for labels
  x <- x %>%
    mutate(label = case_when(
      type %in% Domain ~ "Family & Domain",
      type %in% Structure ~ "Structure",
      TRUE ~ "Features"
    )) %>%
    group_by(label) %>%
    mutate(type = if_else(label == "Family & Domain",
                          make.unique(as.character(type), sep = "_"), type)) %>%
    mutate(across(c(Note, Dbxref), ~ str_replace_all(., "character\\(0\\)", ""))) %>%
    mutate(across(c(Note, Dbxref, evidence), ~ str_wrap(., width = 30)))
  
  # pre-wrap for better spacing
  x <- x %>%
    mutate(Note = str_wrap(Note, width = 25))
  
  ggplot(x, aes(y = type, x = start)) +
    geom_segment(aes(xend = end, yend = type, colour = type),
                 linewidth = 4, show.legend = FALSE) +
    facet_grid(vars(label), scales = "free", space = "free") +
    ggrepel::geom_text_repel(
      data = x %>% filter(label == "Family & Domain"),
      aes(label = Note, x = position_label, y = type),
      size = 2.3,
      direction = "both",
      lineheight = 0.9,
      box.padding = unit(1.0, "lines"),
      point.padding = unit(0.5, "lines"),
      segment.size = 0.3,
      segment.colour = "grey60",
      min.segment.length = 0,
      force = 10,
      force_pull = 0.5,
      nudge_y = 0.8,
      max.overlaps = Inf
    ) +
    geom_vline(data = pos, aes(xintercept = Protein_position),
               colour = "black", linewidth = 1) +
    ylab("") + xlab("Protein position") +
    labs(
      title = paste0(symbol, " evidence tracks"),
      subtitle = variant_info
    ) +
    theme_bw() +
    theme(
      panel.background = element_blank(),
      plot.subtitle = element_text(size = 10, colour = "grey20", face = "italic")
    ) +
    scale_color_manual(values = wes_pal)
}

#-------------------------------------------------------------
# Run and save
#-------------------------------------------------------------
p <- create_plot()
print(p)

ggsave(paste0("../images/", filename, ".pdf"), plot = p,
       height = 7, width = 8, limitsize = FALSE)

message("Saved plot")
