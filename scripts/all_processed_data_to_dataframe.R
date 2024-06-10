## Convert processed data to unified dataframe

# Compress rows to gene level ----
# Mouse knockouts ----
mouse.viability.impc3

mouse.phenotypes.impc3 <- mouse.phenotypes.impc2 %>%
  pivot_wider(names_from = impc_zygosity, values_from = impc_phenotypes,
              names_prefix = "phenotype_", values_fill = list(impc_phenotypes = NA)) %>%
  dplyr::rename(homozygote_phenotypes = phenotype_homozygote,
                heterozygote_phenotypes = phenotype_heterozygote) %>%
  dplyr::select(gene_symbol, homozygote_phenotypes, heterozygote_phenotypes) %>%
  group_by(gene_symbol) %>%
  summarise(
    homozygote_phenotypes = paste(na.omit(homozygote_phenotypes), collapse = " | "),
    heterozygote_phenotypes = paste(na.omit(heterozygote_phenotypes), collapse = " | ")
  ) %>%
  mutate_all(., ~ ifelse(. == "", NA, .)) %>%
  mutate_all(., ~ ifelse(. == " ", NA, .))

wol.categories3 <- wol.categories2 %>%
  dplyr::select(-mgi_id)

mouse.viability.mgi4

orthologs

# Disease data ----
genemap9 <- genemap8 %>%
  mutate_all(., ~ ifelse(. == "NA", NA, .)) %>%
  filter(!is.na(phenotypes)) %>%
  dplyr::rename(omim_phenotype = phenotypes) %>%
  left_join(lethal.phenotypes2) %>%
  mutate(omim_phenotype = paste(omim_phenotype, " (", phenotype_id, ")", " (", moi, ")", " (", na.omit(disease_gene_lethal), ")", " (", na.omit(earliest_lethality_category), ")", sep = "")) %>%
  dplyr::rename(gene_symbol = `Approved Gene Symbol`) %>%
  dplyr::select(gene_symbol, number_key, omim_phenotype) %>% 
  mutate(omim_phenotype = gsub("\\(\\)", "", omim_phenotype)) %>%
  group_by(gene_symbol) %>%
  summarise(
    omim_phenotype = paste(na.omit(omim_phenotype), collapse = " | "),
    number_key = paste(na.omit(number_key), collapse = " | ")
    )

ddg2p_5 <- ddg2p_4 %>%
  group_by(gene_symbol) %>%
  summarise(
    organ_specificity = paste(na.omit(organ_specificity), collapse = " | "),
    disease_name = paste(na.omit(disease_name), collapse = " | "),
    allelic_requirement = paste(na.omit(allelic_requirement), collapse = " | "),
    confidence_category = paste(na.omit(confidence_category), collapse = " | ")
    )

# Constraint metrics ----
constraint

# Panther protein data ----
panther.data2

# Reactome pathway data ----
reactome.annotations

# GO data ----
go.data2

# Join all ----
pcg <- read.fst('./data/raw/protein.coding.genes.fst') %>%
  dplyr::select(symbol, name, gene_group, hgnc_id) %>%
  dplyr::rename(gene_symbol = symbol) %>%
  mutate_all(., ~ ifelse(. == "", NA, .))

all_data <- pcg %>%
  left_join(mouse.viability.impc3) %>%
  left_join(mouse.phenotypes.impc3) %>%
  left_join(wol.categories3) %>%
  left_join(mouse.viability.mgi4) %>%
  left_join(orthologs) %>%
  left_join(genemap9) %>%
  left_join(ddg2p_5) %>%
  left_join(constraint) %>%
  left_join(panther.data2) %>%
  left_join(reactome.annotations) %>%
  left_join(go.data2)

write.fst(all_data, './data/processed/annotated.protein.coding.genes.fst')
