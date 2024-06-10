## Pantherdb annotations 
panther_data_shortcols <- read.fst('./data/raw/panther.fst')
panther.data <- panther_data_shortcols %>%
  distinct(UNIPROT, .keep_all = TRUE) %>%
  dplyr::rename(uniprot_id = UNIPROT)

pcg <- read.fst('./data/raw/protein.coding.genes.fst') %>%
  dplyr::select(symbol, uniprot_ids) %>%
  dplyr::rename(gene_symbol = symbol, uniprot_id = uniprot_ids)

panther.data2 <- panther.data %>%
  mutate(CLASS_TERM = paste(CLASS_TERM, " (", CLASS_ID, ")", sep = "")) %>%
  dplyr::select(-FAMILY_ID, -FAMILY_TERM, -CLASS_ID) %>%
  left_join(pcg) %>%
  dplyr::select(-uniprot_id) %>%
  filter(!is.na(gene_symbol)) %>%
  distinct()
# final = panther.data2