# Create EGE df
# Look at all_data and seperate rows and might need to add a unique() to the data in the df
mouse.viability.impc33 <- mouse.viability.impc3 %>%
  select(gene_symbol, impc_viability) %>%
  mutate(
    impc_viability = as.factor(impc_viability)
  )
mouse.phenotypes.impc22 <- mouse.phenotypes.impc2 %>%
  select(gene_symbol, impc_zygosity, impc_phenotypes)%>%
  mutate(
    impc_zygosity = as.factor(impc_zygosity)
  )
wol.categories22 <- wol.categories2 %>%
  select(gene_symbol, wol) %>%
  mutate(
    wol = as.factor(wol)
  )
mouse.viability.mgi44 <- mouse.viability.mgi4 %>%
  mutate(
    mgi_viability = as.factor(mgi_viability)
  )
# orthologs2 <- orthologs %>%
#   mutate(
#     ortholog_mapping = as.factor(ortholog_mapping)
#   )
orthologs2 <- orthologs %>%
  rename(gene_symbol = V1, ortholog_mapping = `Human Support Count Threshold`) %>%
  select(gene_symbol, ortholog_mapping) %>%
  distinct() %>%
  mutate(
    ortholog_mapping = as.factor(ortholog_mapping)
  )

genemap88 <- genemap8 %>%
  mutate(across(everything(), ~ ifelse(. == "NA", NA, .))) %>%
  mutate(
    moi = as.factor(moi),
    symbol_key = as.factor(symbol_key),
    number_key = as.factor(number_key)
  ) %>%
  rename(gene_symbol = `Approved Gene Symbol`)
ddg2p_44 <- ddg2p_4 %>%
  mutate(
    allelic_requirement = as.factor(allelic_requirement),
    confidence_category = as.factor(confidence_category),
    organ_specificity = as.factor(organ_specificity)
  )
# add genes, remove genes with ambiguous disease_gene_lethal value (check how many we retain/lose)
lethal.phenotypes22 <- genemap88 %>%
  rename(omim_phenotype = phenotypes) %>%
  left_join(lethal.phenotypes2, by = 'omim_phenotype')

# only keep genes with lethal phenotypes
lethal.phenotypes222 <- lethal.phenotypes22 %>%
  select(gene_symbol, disease_gene_lethal) %>%
  distinct() %>%
  count(gene_symbol) %>%
  filter(n == 1)

lethal.phenotypes2222 <- lethal.phenotypes22 %>%
  filter(gene_symbol %in% lethal.phenotypes222$gene_symbol) %>%
  select(gene_symbol, disease_gene_lethal, earliest_lethality_category) %>%
  mutate(
    disease_gene_lethal = as.factor(disease_gene_lethal),
    earliest_lethality_category = as.factor(earliest_lethality_category)
  )
  
constraint

all_data_ege <- mouse.viability.impc33 %>%
  full_join(mouse.phenotypes.impc22) %>%
  full_join(wol.categories22) %>%
  full_join(mouse.viability.mgi44) %>%
  full_join(orthologs2) %>%
  full_join(genemap88) %>%
  full_join(lethal.phenotypes2222) %>%
  full_join(ddg2p_44) %>%
  full_join(constraint) 

# Reorder and rename for app
names(all_data_ege)
reordered_list <- c(
  'gene_symbol', 'name', 'gene_group', 'entrez_id', 
  'impc_viability', 'impc_zygosity', 'impc_phenotypes', 'wol', 'mgi_viability', 'ortholog_mapping', 
  'moi', 'symbol_key', 'number_key', 'phenotype_id', 'phenotypes', 'disease_gene_lethal', 'earliest_lethality_category', 'disease_name', 'allelic_requirement', 'confidence_category', 'organ_specificity',
  'percentage_essential', 'mean_score_all', 'bf_mef', 'bf_lam', 'gnomad_lof_upper_90_ci', 'gnomad_mis_upper_90_ci', 'mean_am_pathogenicity', 'shet_rgcme_mean', 'shet_post_mean', 'domino', 'scones'
  )
reordered_ege_df <- all_data_ege %>% select(all_of(reordered_list))
names(reordered_ege_df)
  
new_col_names <- c(
  'Gene symbol', 'Gene name', 'Gene group', 'Entrez ID', 
  'IMPC viability', 'IMPC zygosity', 'IMPC phenotypes', 'Windows of Lethality', 'MGI viability', 'Ortholog mapping',
  'Mode of Inheritance', 'Disease type', 'Molecular basis', 'OMIM Phenotype ID', 'Phenotype', 'Lethal gene', 'Earliest age of death category', 'Disease', 'Allelic requirement', 'Confidence category', 'Organ specificity',
  'Percentage of samples essential in', 'Mean DepMap Gene Effect score', 'Bayes factor (MEF)', 'Bayes factor (Laminin)', 
  'gnomAD LOEUF', 'gnomAD missense score', 'AlphaMissense mean pathogenicity', 'Shet (RGCME)', 'Shet (GeneBayes)', 'DOMINO', 'SCoNeS'
  )

names(reordered_ege_df) <- new_col_names

write.fst(reordered_ege_df, './data/processed/all_data_ege.fst')
