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

orthologs2

# Disease data ----
genemap9 <- genemap8 %>%
  mutate_all(~ ifelse(. == "NA", NA, .)) %>%
  filter(!is.na(phenotypes)) %>%
  dplyr::rename(omim_phenotype = phenotypes) %>%
  left_join(lethal.phenotypes2) %>%
  mutate(omim_phenotype = pmap_chr(list(omim_phenotype, phenotype_id, moi, disease_gene_lethal, earliest_lethality_category), 
                                   function(omim_phenotype, phenotype_id, moi, disease_gene_lethal, earliest_lethality_category) {
                                     paste(omim_phenotype, 
                                           if (!is.na(phenotype_id)) paste0(" (", phenotype_id, ")") else "", 
                                           if (!is.na(moi)) paste0(" (", moi, ")") else "", 
                                           if (!is.na(disease_gene_lethal)) paste0(" (", disease_gene_lethal, ")") else "", 
                                           if (!is.na(earliest_lethality_category)) paste0(" (", earliest_lethality_category, ")") else "", 
                                           sep = "")
                                   })) %>%
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
constraint2 <- constraint %>%
  dplyr::select(-entrez_id, -name, -gene_group)

# Panther protein data ----
panther.data2

# Reactome pathway data ----
reactome.annotations

# GO data ----
go.data2

# HCA & GTEx gene expression ----
HCA_expression_data2 <-  HCA_expression_data %>%
  group_by(gene_symbol) %>%
  mutate(HCA_high_expression_in = paste(tissue, " (", cell_type, ")", sep = "")) %>%
  summarise(
    HCA_high_expression_in = paste(na.omit(HCA_high_expression_in), collapse = " | ")
    )

GTEx_expression_data3 <- GTEx_expression_data2 %>%
  mutate(GTEx_average_TPM_across_tissues = round(average_TPM_across_tissues, 3)) %>%
  dplyr::select(-average_TPM_across_tissues)

expression_data <- HCA_expression_data2 %>%
  full_join(GTEx_expression_data3)

# STRING PPI data ----
string_ppi_df2 <- string_ppi_df %>%
  mutate(string_ppi = paste(protein1_string_id, "&", protein1_string_id, ":", combined_score, sep = "")) %>%
  group_by(gene_symbol) %>%
  summarise(
    string_ppi = paste(na.omit(string_ppi), collapse = " | ")
  ) %>%
  dplyr::select(gene_symbol, string_ppi)

# Join all ----
pcg <- read.fst('./data/raw/protein.coding.genes.fst') %>%
  dplyr::select(symbol, name, gene_group, hgnc_id) %>%
  dplyr::rename(gene_symbol = symbol) %>%
  mutate_all(., ~ ifelse(. == "", NA, .))

all_data <- pcg %>%
  left_join(mouse.viability.impc3) %>%
  left_join(mouse.phenotypes.impc3) %>%
  left_join(wol.categories3) %>%
  # left_join(mouse.viability.mgi4) %>%
  left_join(orthologs2) %>%
  left_join(genemap9) %>%
  left_join(ddg2p_5) %>%
  left_join(constraint2) %>%
  left_join(panther.data2) %>%
  left_join(reactome.annotations) %>%
  left_join(go.data2) %>%
  left_join(expression_data) %>%
  left_join(string_ppi_df2) %>%
  distinct()

all_data <- all_data %>%
  relocate(bf_mef, bf_lam, .after = mean_score_all)

nice_col_names <- c(
  'Gene symbol', 'Gene name', 'Gene group', 'HGNC ID',
  'IMPC Viability', 'IMPC Homozygote phenotypes', 'IMPC Heterozygote phenotypes', 'Window of Lethality',
  # 'MGI viability',
  'Human-Mouse ortholog',
  'OMIM Phenotype', 'Molecular basis',
  'Affected organ', 'Disease', 'Allelic requiremnet', 'Gene-Disease association confidence', 
  'Percentage essential lines', 'Mean CERES score', 'Bayes factor (MEF)', 'Bayes factor (Laminin)',
  'gnomAD LOEUF', 'gnomAD Mis', 'Alpha Missense', 'Shet (GeneBayes)', 'DOMINO', 'SCoNeS', 'GISMO median', 'GISMO decile','GISMO-mis median', 'GISMO-mis decile',
  'PANTHER Class', 'PANTHER Subfamily', 
  'Reactome pathway',
  'GO Biological process', 'GO Molecular function', 'GO Cellulcar component',
  'HCA Highly expressed tissues', 'GTEx average expression across tissues',
  'STRING PPIs'
)

names(all_data) <- nice_col_names


file_path <- paste0('./data/processed/annotated.protein.coding.genes_', Sys.Date(), '.fst')
write.fst(all_data, file_path)
