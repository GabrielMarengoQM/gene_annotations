# Gene ontology
# Retrieve GO Ids for each Ontology (BP, MF, CC) and Entrez IDs
go.data.raw <- read.fst("./data/raw/go.fst")

go.data <- go.data.raw %>%
  group_by(gene_id, Ontology) %>%
  summarise(go_ids = paste0(go_id, collapse = "|"),
            go_terms = paste0(go_term, collapse = "|")) %>%
  pivot_wider(names_from = Ontology,
              values_from = c("go_ids","go_terms")) %>%
  dplyr::select(gene_id, go_ids_BP, go_terms_BP,
                go_ids_MF, go_terms_MF,
                go_ids_CC, go_terms_CC) %>%
  dplyr::rename(entrez_id = gene_id)

# BP table
go.bp <- go.data %>%
  dplyr::select(entrez_id, go_ids_BP, go_terms_BP) %>%
  separate_rows(go_ids_BP, go_terms_BP, sep = "\\|") %>%
  mutate(go_bp = paste(go_terms_BP, " (", go_ids_BP, ")", sep = "")) %>%
  dplyr::select(entrez_id, go_bp) %>%
  group_by(entrez_id) %>%
  summarise(go_bp = paste0(go_bp, collapse = " | ")) %>%
  filter(go_bp != 'NA (NA)')

# MF table
go.mf <- go.data %>%
  dplyr::select(entrez_id, go_ids_MF, go_terms_MF) %>%
  separate_rows(go_ids_MF, go_terms_MF, sep = "\\|") %>%
  mutate(go_mf = paste(go_terms_MF, " (", go_ids_MF, ")", sep = "")) %>%
  dplyr::select(entrez_id, go_mf) %>%
  summarise(go_mf = paste0(go_mf, collapse = " | ")) %>%
  filter(go_mf != 'NA (NA)')

# CC table
go.cc <- go.data %>%
  dplyr::select(entrez_id, go_ids_CC, go_terms_CC) %>%
  separate_rows(go_ids_CC, go_terms_CC, sep = "\\|") %>%
  mutate(go_cc = paste(go_terms_CC, " (", go_ids_CC, ")", sep = "")) %>%
  dplyr::select(entrez_id, go_cc) %>%
  summarise(go_cc = paste0(go_cc, collapse = " | ")) %>%
  filter(go_cc != 'NA (NA)')

# Join
pcg <- read.fst('./data/raw/protein.coding.genes.fst') %>%
  dplyr::select(symbol, entrez_id) %>%
  dplyr::rename(gene_symbol = symbol)
pcg$entrez_id <- as.character(pcg$entrez_id)

go.data2 <- go.bp %>%
  left_join(go.mf) %>%
  left_join(go.cc)  %>%
  left_join(pcg) %>%
  dplyr::select(-entrez_id) %>%
  filter(!is.na(gene_symbol))


# final = go.data2





