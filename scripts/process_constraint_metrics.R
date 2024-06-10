## Process Cell line and Population Sequencing Constraint metrics
protein.coding.genes <- read.fst("./data/raw/protein.coding.genes.fst")

# Depmap
models <- read.fst("./data/raw/models.fst")
cns_brain <- models %>% filter(OncotreeLineage == "CNS/Brain")

# cell0
# Transform and turn into matrix of 1 and 0 where 1 is essential and 0 is not (-0.5 gene effect threshold)
cell0 <- read.fst("./data/raw/cell0.fst")
colnames(cell0) <-  str_split(names(cell0), "\\ ", simplify=T)[,1]
colnames(cell0)[colnames(cell0) == "...1"] <- "model_id"
cell1  = cell0[,-1]
cell2 <- t(cell1)
cell3 <- cell2
colnames(cell3) = cell0$V1
cell4 = cell3
cell4[cell4 >= -0.5] <- 0
cell4[cell4 != 0] <- 1
cell4 <- data.frame(cell4)
# cell 5 = with n_essential, percentage_essential columns
cell5 <- cell4 %>%
  mutate(n_essential = rowSums(.)) %>%
  mutate(percentage_essential = (n_essential/ncol(cell4))*100)
# cell 7 = mean scores
cell7 <- data.frame(cell3) %>%
  mutate(mean_score_all = rowMeans(.)) %>%
  mutate(gene_symbol = row.names(cell3)) %>%
  dplyr::select(gene_symbol, mean_score_all) %>%
  `rownames<-`( NULL )

depmap.data1 <- cell5 %>%
  mutate(gene_symbol = row.names(.)) %>%
  dplyr::select(gene_symbol, percentage_essential) %>%
  mutate(gene_symbol = gsub("\\.\\..*?\\.", "", gene_symbol))

depmap.data2 <- cell7 %>%
  dplyr::select(gene_symbol, mean_score_all) %>%
  mutate(gene_symbol = gsub("\\.\\..*?\\.", "", gene_symbol))

depmap.data3 <- depmap.data1 %>%
  full_join(depmap.data2, by = "gene_symbol") %>%
  mutate(percentage_essential = round(percentage_essential, 3),
         mean_score_all = round(mean_score_all, 3))

# Mef laminin 
mef <- read.fst("./data/raw/mef.fst")
mef <- mef[,c(1, 2, 7)]
names(mef) <- c("gene_symbol", "bf_mef", "fdr_mef")
mef <- mef %>%
  mutate(fdr_mef = round(fdr_mef, 3)) 

laminin <- read.fst("./data/raw/laminin.fst")
laminin <- laminin[,c(1, 2, 7)]
names(laminin) <- c("gene_symbol", "bf_lam", "fdr_lam")
laminin <- laminin %>%
  mutate(fdr_lam = round(fdr_lam, 3)) 

# gnomAD ---- CONSIDER THE METRICS WE KEEP --- What to do with dups?
# Get protein coding genes
# Get gnomad v4 data
gnomad_v4 <- read.fst("./data/raw/gnomad_v4.fst")
gnomad_v4 <- gnomad_v4[, c(1, 2, 3, 21, 30, 46)]
names(gnomad_v4) <- c("gene", "gnomad_transcript", "gnomad_mane_select", 
                      "gnomad_lof_upper_90_ci", "gnomad_mis_upper_90_ci", 'gnomad_constraint_flags')

# Get genes without a mane transcript (715 genes)
gnomad_v4_genes_with_no_mane_transcript <- gnomad_v4 %>%
  group_by(gene) %>% 
  summarise(gnomad_mane_select = paste(gnomad_mane_select, collapse=", ")) %>%
  filter(!grepl("true", gnomad_mane_select))

# Subset original data for only genes with no mane transcript
gnomad_v4_no_mane_df <- gnomad_v4 %>%
  filter(gene %in% gnomad_v4_genes_with_no_mane_transcript$gene)
length(unique(gnomad_v4_no_mane_df$gene))

# Get all Ensembl canonical transcripts
# mart <- useEnsembl("ensembl", dataset="hsapiens_gene_ensembl")
# BM.info <- getBM(attributes=c("ensembl_gene_id","ensembl_transcript_id","hgnc_symbol", "hgnc_id", "transcript_is_canonical"), mart=mart)
canon_transcripts <- read.fst("./data/raw/canon_transcripts.fst")
canonical_transcripts <- subset(canon_transcripts, transcript_is_canonical == 1) %>%
  dplyr::rename(transcript = ensembl_transcript_id) %>%
  dplyr::rename(gene = hgnc_symbol) %>%
  dplyr::select(gene, transcript)

# Subset the genes with no mane transcript for only canonical transcripts (550 genes -> 166 genes with no canonical transcript -> 86 of which are protein coding)
gnomad_no_mane_canon_transcript_only <- subset(gnomad_v4_no_mane_df, gnomad_transcript %in% canonical_transcripts$transcript)

# Genes with no canonical transcript (166)
lost_genes <- dplyr::setdiff(gnomad_v4_genes_with_no_mane_transcript$gene, gnomad_no_mane_canon_transcript_only$gene)
# Find how many are protein coding (86)
lost_genes_df <- data.frame(lost_genes) %>%
  mutate(is_protein_coding = ifelse(lost_genes %in% protein.coding.genes$symbol, 'yes', 'no')) 
lost_genes_count <- lost_genes_df %>%
  group_by(is_protein_coding) %>%
  tally()
protein_coding_lost_genes <- lost_genes_df %>%
  filter(is_protein_coding == 'yes')

# Check that the 86 protein coding genes with supposedly no canonical transcript are not in the list of all canonical transcripts - we should expect 0 genes and this is what we get
transcripts_to_check <- gnomad_v4_no_mane_df %>%
  subset(gene %in% protein_coding_lost_genes$lost_genes)
intersecting_genes <- base::intersect(transcripts_to_check$gnomad_transcript, canonical_transcripts$transcript)

# How to select for these 86 gene with no canonical?

# Create final dataset with 2 cols flags (1. missing data reason & 2. transcript type)
# MANE stands for Matched Annotation from NCBI and EMBL-EBI and MANE Select identifies 
#   one high-quality representative transcript per protein-coding gene that is well-supported 
#   by experimental data and represents the biology of the gene. More background can be found on the MANE project's website.
gnomad_v4_mane_only <- gnomad_v4[gnomad_v4$gnomad_mane_select == 'true',]
gnomad_v4_mane_and_ens_id_only <- gnomad_v4_mane_only[grepl('^ENST', gnomad_v4_mane_only$gnomad_transcript), ]

gnomad_v4_mane <- gnomad_v4_mane_and_ens_id_only %>%
  mutate(flag = 'gnomad_mane_transcript')

gnomad_v4_canon <- gnomad_v4 %>%
  filter(gnomad_transcript %in% gnomad_no_mane_canon_transcript_only$gnomad_transcript) %>%
  mutate(flag = 'ensembl_canonical_transcript')

gnomad_with_flags <- bind_rows(gnomad_v4_mane, gnomad_v4_canon)

gnomad.v4.data <- gnomad_with_flags
gnomad.v4.data$gnomad_constraint_flags <- gsub("\\[|\\]", "", gnomad.v4.data$gnomad_constraint_flags)
gnomad.v4.data <- gnomad.v4.data %>%
  dplyr::rename(gene_symbol = gene) %>%
  filter(!is.na(gene_symbol))
c.gnomad <- gnomad.v4.data %>%
  count(gene_symbol) %>%
  filter(n > 1)
gnomad.v4.data2 <- gnomad.v4.data %>%
  filter(!gene_symbol %in% c.gnomad$gene_symbol)

# Alpha missense
alpham2 <- read.fst("./data/raw/alpham2.fst") %>%
  mutate(ensembl_transcript_id = gsub("\\.\\d+", "", transcript_id)) %>%
  dplyr::select(-transcript_id)

# # get hgnc ens transcript id mappings from biomart 
# ensembl <- useMart("ENSEMBL_MART_ENSEMBL",
#                    dataset="hsapiens_gene_ensembl")
# values1 <- alpham2$ensembl_transcript_id
# att1 <- c("ensembl_transcript_id", "hgnc_id")
# filter1 <- c("ensembl_transcript_id")
# searchResults <-getBM(att1,
#                       filters= filter1,
#                       values= values1, mart=ensembl,
#                       checkFilters = FALSE, # dont think this is needed
#                       verbose = TRUE)
ens_transcript_hgnc_id_mappings <- read.fst("./data/raw/ens_transcript_hgnc_id_mappings.fst")

# Add hgnc_id
alpha_missense <- alpham2 %>%
  left_join(ens_transcript_hgnc_id_mappings)

alpha_missense_2 <- alpha_missense %>%
  filter(!is.na(hgnc_id) & hgnc_id != "")

# Filter duplicated genes for canonical
c.am <- alpha_missense_2 %>%
  count(hgnc_id) %>%
  filter(n > 1)
# canon transcrpts 
canonical_transcripts <- subset(canon_transcripts, transcript_is_canonical == 1) %>%
  dplyr::rename(transcript = ensembl_transcript_id) %>%
  dplyr::rename(gene = hgnc_symbol) %>%
  dplyr::select(hgnc_id, transcript)
# get canon transcripts with am score
am_dups <- alpha_missense_2 %>%
  filter(hgnc_id %in% c.am$hgnc_id)
am_canon <- am_dups %>%
  filter(ensembl_transcript_id %in% canonical_transcripts$transcript)

alpha_missense_3 <- alpha_missense_2 %>%
  filter(!hgnc_id %in% am_dups$hgnc_id) %>%
  bind_rows(am_canon) %>%
  mutate(mean_am_pathogenicity = round(mean_am_pathogenicity, 3))

# shet rgcme ----
shet_rgcme <- read.fst("./data/raw/shet_rgcme.fst")
shet_rgcme <- shet_rgcme[, c(1, 2, 4, 6, 7)]
names(shet_rgcme) <- c("gene_symbol", "ens_gene_id", "shet_rgcme_mean", "shet_rgcme_lower", "shet_rgcme_upper")
shet_rgcme <- shet_rgcme %>%
  mutate(
    shet_rgcme_mean = round(shet_rgcme_mean, 3),
    shet_rgcme_lower = round(shet_rgcme_lower, 3),
    shet_rgcme_upper = round(shet_rgcme_lower, 3)
  )

# shet post ----
shet_post <- read.fst("./data/raw/shet_post.fst")
shet_post <- shet_post[, c(1, 2, 7, 8, 9)]
names(shet_post) <- c("ens_gene_id", "hgnc_id", "shet_post_mean", "shet_post_lower", "shet_post_upper")
shet_post <- shet_post %>%
  mutate(
    shet_post_mean = round(shet_post_mean, 3),
    shet_post_lower = round(shet_post_lower, 3),
    shet_post_upper = round(shet_post_upper, 3)
  ) %>%
  dplyr::select(-ens_gene_id)

# scones and domino ----
scones_domino <- read.fst("./data/raw/scones_domino.fst")
scones_domino <- scones_domino[, c(1, 17, 19)] %>%
  dplyr::rename(gene_symbol = Gene)
scones_domino.2 <- subset(scones_domino, gene_symbol %in% protein.coding.genes$symbol)
scones_domino.2$DOMINO <- as.numeric(scones_domino.2$DOMINO)
scones_domino <- scones_domino.2 %>%
  dplyr::rename(domino = DOMINO) %>%
  dplyr::rename(scones = SCoNeS) %>%
  mutate(domino = round(domino, 3))

# combine
pcg <- protein.coding.genes %>%
  dplyr::select(symbol, hgnc_id, entrez_id, name, gene_group) %>%
  dplyr::rename(gene_symbol = symbol)

constraint <- pcg %>%
  left_join(depmap.data3) %>%
  left_join(gnomad.v4.data2) %>%
  left_join(alpha_missense_3) %>%
  left_join(shet_rgcme) %>%
  left_join(shet_post) %>%
  left_join(scones_domino) %>%
  left_join(mef) %>%
  left_join(laminin) %>%
  dplyr::select(-gnomad_transcript, -gnomad_mane_select, -gnomad_constraint_flags, -flag, -ensembl_transcript_id, 
                -ens_gene_id, -shet_rgcme_lower, -shet_rgcme_upper, -shet_post_lower, -shet_post_upper, -hgnc_id, -fdr_lam, -fdr_mef)
constraint$entrez_id <- as.character(constraint$entrez_id)
# final = constraint

