## Disease phenotypes and lethality from OMIM & DDG2P

# OMIM genemap ----
genemap <- read.fst("./data/raw/genemap.fst")

# Split phenotypes across new rows for one2many gene-pheno mappings
genemap_2 <- genemap %>%
  separate_rows(Phenotypes, sep = "; ")

# Grep mois
moi_keywords <- c("Autosomal recessive", "Autosomal dominant", "Digenic recessive", "Digenic dominant",
                  "Somatic mutation", "Multifactorial", "Isolated cases", "Mitochondrial",
                  "Pseudoautosomal dominant", "Pseudoautosomal recessive", "X-linked recessive",
                  "X-linked dominant", "X-linked", "Y-linked")

genemap_3 <- genemap_2 %>%
  mutate(moi = sapply(str_extract_all(Phenotypes, paste(moi_keywords, collapse = "|")), 
                      paste, collapse = "; ")) %>%
  separate_rows(moi, sep = "; ")

# Grep symbols
genemap_4 <- genemap_3 %>%
  mutate(symbol_key = case_when(
    grepl("^\\{", Phenotypes) ~ "Mutations that contribute to susceptibility to multifactorial disorders",
    grepl("^\\[", Phenotypes) ~ "Nondiseases",
    grepl("^\\?", Phenotypes) ~ "The relationship between the phenotype and gene is provisional",
    TRUE ~ NA_character_  # Handle other cases if needed
  ))

# Grep numbers
genemap_5 <- genemap_4 %>%
  mutate(number_key = case_when(
    grepl("\\(1)", Phenotypes) ~ "The disorder is placed on the map based on its association with a gene, but the underlying defect is not known",
    grepl("\\(2)", Phenotypes) ~ "The disorder has been placed on the map by linkage or other statistical method; no mutation has been found",
    grepl("\\(3)", Phenotypes) ~ "The molecular basis for the disorder is known; a mutation has been found in the gene",
    grepl("\\(4)", Phenotypes) ~ "A contiguous gene deletion or duplication syndrome, multiple genes are deleted or duplicated causing the phenotype",
    TRUE ~ NA_character_  # Handle other cases if needed
  ))

# Get phenotype id
genemap_6 <- genemap_5 %>%
  mutate(phenotype_id = str_extract(Phenotypes, "\\d{6}"))

# Get Phenotype
genemap_7 <- genemap_6 %>%
  mutate(phenotypes = str_extract(Phenotypes, "^(.*?)(?=\\d{6})")) %>%
  mutate(phenotypes = str_replace_all(phenotypes, "[\\{\\[\\?\\]\\}]", "")) %>%
  mutate(phenotypes = str_replace(phenotypes, ", $", "")) 

genemap8 <- genemap_7 %>%
  dplyr::select(`Approved Gene Symbol`, moi, symbol_key, number_key, phenotype_id, phenotypes) %>%
  filter(!is.na(`Approved Gene Symbol`))
# final = genemap8

# Lethal phenotypes ----
lethal.phenotypes <- read.fst("./data/raw/lethal.phenotypes.fst")

# lethal.phenotypes2 <- lethal.phenotypes %>%
#   dplyr::select(
#     omim_phenotype,
#     disease_gene_lethal,
#     earliest_lethality_category
#   ) %>%
#   mutate(
#          earliest_lethality_category = case_when(
#            is.na(earliest_lethality_category) ~ NA_character_,
#            earliest_lethality_category == "L1" ~ "Prenatal death (L1)",
#            earliest_lethality_category == "L2" ~ "Neonatal death (L2)",
#            earliest_lethality_category == "L3" ~ "Death in infancy (L3)",
#            earliest_lethality_category == "L4" ~ "Death in childhood (L4)",
#            earliest_lethality_category == "L5" ~ "Death in adolescence (L5)",
#            earliest_lethality_category == "L6" ~ "Death in adulthood (L6)",
#            earliest_lethality_category == "LU" ~ "Not determined (LU)",
#            earliest_lethality_category == "NL" ~ "Non lethal (NL)",
#            TRUE ~ earliest_lethality_category
#          )) %>%
#   mutate_all(., ~ ifelse(. == "-", NA, .))
lethal.phenotypes2 <- lethal.phenotypes %>%
  dplyr::select(
    omim_phenotype,
    disease_gene_lethal,
    earliest_lethality_category
  ) %>%
  mutate(
    earliest_lethality_category = case_when(
      is.na(earliest_lethality_category) ~ NA_character_,
      earliest_lethality_category == "L1" ~ "Prenatal death; L1",
      earliest_lethality_category == "L2" ~ "Neonatal death; L2",
      earliest_lethality_category == "L3" ~ "Death in infancy; L3",
      earliest_lethality_category == "L4" ~ "Death in childhood; L4",
      earliest_lethality_category == "L5" ~ "Death in adolescence; L5",
      earliest_lethality_category == "L6" ~ "Death in adulthood; L6",
      earliest_lethality_category == "LU" ~ "Not determined; LU",
      earliest_lethality_category == "NL" ~ "Non lethal; NL",
      TRUE ~ earliest_lethality_category
    )) %>%
  mutate_all(., ~ ifelse(. == "-", NA, .))
# final = lethal_phenotypes2

# DDG2P ----
ddg2p <- read.fst("./data/raw/dd.g2p.fst")

ddg2p_2 <- ddg2p %>%
  tidyr::separate_rows(`organ specificity list`, sep = ";")

ddg2p_3 <- ddg2p_2 %>%
  dplyr::select(`gene symbol`, `disease name`, `allelic requirement`, `confidence category`, `organ specificity list`) %>%
  dplyr::rename(gene_symbol = `gene symbol`)

ddg2p_4 <- ddg2p_3 %>%
  dplyr::rename(
    disease_name = `disease name`,
    allelic_requirement = `allelic requirement`,
    confidence_category = `confidence category`,
    organ_specificity = `organ specificity list`
  ) %>%
  distinct()
# final = ddg2p_4