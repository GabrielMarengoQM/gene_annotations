## Mouse viability and phenotypes from knockouts

# Human mouse gene mappings from EBI ----
pcg <- read.fst("./data/raw/protein.coding.genes.fst")

human_mouse_genes <- pcg %>%
  separate_rows(mgd_id, sep = "\\|") %>%
  dplyr::select(symbol, mgd_id) %>%
  dplyr::rename(gene_symbol = symbol, mgi_id = mgd_id)

# IMPC viability ----
mouse.viability.impc <- read.fst("./data/raw/mouse.viability.impc.fst") 

mouse.viability.impc2 <- mouse.viability.impc %>%
  dplyr::filter(Comment == "") %>% # Remove conflicting evidence
  dplyr::select('Gene Accession Id', 'Viability Phenotype HOMs/HEMIs') %>%
  dplyr::rename(mgi_id = 'Gene Accession Id',
                impc_viability = 'Viability Phenotype HOMs/HEMIs') %>%
  distinct() %>%
  left_join(human_mouse_genes) %>%  # Join gene symbols to mouse viability
  dplyr::filter(!is.na(gene_symbol)) %>%
  dplyr::select(-mgi_id) %>%
  distinct()
  
conflicts <- mouse.viability.impc2 %>% # remove conflicts that arise from one2many mappings
  distinct() %>%
  count(gene_symbol) %>%
  filter(n == 1)

mouse.viability.impc3 <- mouse.viability.impc2 %>%
  filter(gene_symbol %in% conflicts$gene_symbol) %>%
  distinct()
# final = mouse.viability.impc3 

# IMPC phenotypes ----
mouse.phenotypes.impc <- read.fst("./data/raw/mouse.phenotypes.impc.fst")

mouse.phenotypes.impc2 <- mouse.phenotypes.impc %>%
  dplyr::select(marker_accession_id, zygosity, mp_term_name, mp_term_id) %>%
  dplyr::rename(mgi_id = marker_accession_id,
                impc_zygosity = zygosity,
                impc_phenotypes = mp_term_name,
                mp_id = mp_term_id) %>%
  dplyr::distinct() %>%
  left_join(human_mouse_genes) %>%
  dplyr::filter(!is.na(gene_symbol))
# final = mouse.phenotypes.impc2

# IMPC WoL ----
wol.categories <- read.fst("./data/raw/wol.categories.fst")
wol.categories2 <- wol.categories %>%
  dplyr::select(mgi_id, wol) %>%
  left_join(human_mouse_genes) %>%
  dplyr::filter(!is.na(gene_symbol))
# final = wol.categories2

# MGI Viability ----
mouse.viability.mgi <- read.fst("./data/raw/mouse.viability.mgi.fst")
mp.lethal.terms <- read.fst("./data/raw/mouse.lethal.terms.fst")

mp.lethal.terms2 <- mp.lethal.terms %>%
  dplyr::select(2) %>%
  dplyr::rename(mp_id = 1)

mouse.viability.mgi2 <- mouse.viability.mgi %>%
  dplyr::select(7,5) %>%
  dplyr::rename(mgi_id = "V7", mp_id = "V5") %>%
  mutate(mp_term_lethal = ifelse(mp_id %in% mp.lethal.terms2$mp_id, "lethal", "viable")) %>%
  dplyr::select(mgi_id, mp_term_lethal) %>%
  distinct() %>%
  left_join(human_mouse_genes)

# remove genes that have both lethal and viable annotations
mouse.viability.mgi3 <- mouse.viability.mgi2 %>%
  dplyr::select(gene_symbol, mp_term_lethal) %>%
  distinct() %>%
  count(gene_symbol) %>%
  filter(n == 1)

mouse.viability.mgi4 <- mouse.viability.mgi2 %>%
  dplyr::filter(gene_symbol %in% mouse.viability.mgi3$gene_symbol) %>%
  dplyr::select(mp_term_lethal, gene_symbol) %>%
  dplyr::rename(mgi_viability = mp_term_lethal) %>%
  distinct()
# final = mouse.viability.mgi4
# MGI Viability - bug fix 22/7/24
mouse.viability.mgi <- read.fst("./data/raw/mouse.viability.mgi.fst")
mp.lethal.terms <- read.fst("./data/raw/mouse.lethal.terms.fst")

mp.lethal.terms2 <- mp.lethal.terms %>%
  dplyr::select(2) %>%
  dplyr::rename(mp_id = 1)

mouse.viability.mgi2 <- mouse.viability.mgi %>%
  dplyr::select(7,5) %>%
  dplyr::rename(mgi_id = "V7", mp_id = "V5") %>%
  mutate(mp_term_lethal = ifelse(mp_id %in% mp.lethal.terms2$mp_id, "lethal", "viable")) %>%
  dplyr::select(mgi_id, mp_term_lethal) %>%
  distinct() %>%
  left_join(human_mouse_genes)

# genes with no conflicts - viable and lethal
mouse.viability.mgi3 <- mouse.viability.mgi2 %>%
  dplyr::select(gene_symbol, mp_term_lethal) %>%
  distinct() %>%
  count(gene_symbol) %>%
  filter(n == 1) %>%
  pull(gene_symbol)

# genes with evidence for lethality (conflicting)
mouse.viability.mgi4 <- mouse.viability.mgi2 %>%
  dplyr::select(gene_symbol, mp_term_lethal) %>%
  distinct() %>%
  count(gene_symbol) %>%
  filter(n > 1) %>%
  mutate(mp_term_lethal = 'lethal') %>%
  pull(gene_symbol)

mouse.viability.mgi5 <- mouse.viability.mgi2 %>%
  dplyr::select(gene_symbol, mp_term_lethal) %>%
  mutate(mp_term_lethal = ifelse(gene_symbol %in% mouse.viability.mgi4, 'lethal', mp_term_lethal)) %>%
  distinct()

mouse.viability.mgi4 <- mouse.viability.mgi5

c <- mouse.viability.mgi4 %>%
  dplyr::count(mp_term_lethal)

# final = mouse.viability.mgi4
# 6/8/24 - Using mouse genes
lethal_terms <- mouse.lethal.terms$x
mouse_proteincoding_genes <- mouse.protein.coding.genes.mgi$`MGI Accession ID`
mouse.viability.mgi2 <- mouse.viability.mgi %>%
  select(7,5) %>%
  distinct() %>%
  rename(mgi_id = V7, mp_term = V5) %>%
  filter(mgi_id %in% mouse_proteincoding_genes) %>%
  mutate(mp_term_lethal = ifelse(mp_term %in% lethal_terms, "y","n")) %>%
  select(mgi_id, mp_term_lethal) %>%
  distinct() %>%
  arrange(mgi_id, mp_term_lethal) %>%
  group_by(mgi_id) %>%
  summarise(mgi_lethal_term = paste0(unique(mp_term_lethal), collapse = "|")) %>%
  mutate(viability_mgi = ifelse(mgi_lethal_term == "n","viable","lethal")) %>%
  select(mgi_id, viability_mgi) %>%
  distinct() 
mouse.viability.mgi3 <- mouse.viability.mgi2 %>% 
  left_join(human_mouse_genes) %>%
  filter(!is.na(gene_symbol))
mouse.viability.mgi4 <- mouse.viability.mgi5
# final = mouse.viability.mgi4
# c <- mouse.viability.mgi3 %>%
#   dplyr::count(viability_mgi)

# Human-Mouse Orthologs ----
orth.mouse.high.conf <- read.fst("./data/raw/human.mouse.orth.fst")

pcg2 <- pcg %>%
  dplyr::select(ensembl_gene_id, symbol) %>%
  dplyr::rename(gene_symbol = symbol)

orthologs <- orth.mouse.high.conf %>%
  left_join(pcg2) %>%
  dplyr::filter(!is.na(gene_symbol)) %>%
  dplyr::rename(ortholog_mapping = mmusculus_homolog_orthology_type) %>%
  dplyr::select(gene_symbol, ortholog_mapping) %>%
  mutate(ortholog_mapping = str_remove(ortholog_mapping, "^ortholog_"))
# final = orthologs

# Lethal gene
mouse.lethal.terms <- read.fst("./data/raw/mouse.lethal.terms.fst")
