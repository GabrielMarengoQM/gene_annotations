reactome_annotations <- read.fst('./data/raw/reactome_annotations.fst')

pcg <- read.fst('./data/raw/protein.coding.genes.fst') %>%
  dplyr::select(symbol, entrez_id) %>%
  dplyr::rename(gene_symbol = symbol)
pcg$entrez_id <- as.character(pcg$entrez_id)

reactome.annotations <- reactome_annotations %>%
  filter(grepl("Homo sapiens:", path_name)) %>%
  mutate(path_name = sub("Homo sapiens: ", "", path_name)) %>%
  dplyr::rename(entrez_id = 'gene_id')  %>%
  mutate(pathway = paste(path_name, " (", DB_ID, ")", sep = "")) %>%
  dplyr::select(entrez_id, pathway) %>%
  group_by(entrez_id) %>%
  summarise(pathway = paste0(pathway, collapse = " | ")) %>%
  left_join(pcg) %>%
  dplyr::select(-entrez_id) %>%
  filter(!is.na(gene_symbol))
# final = reactome.annotations