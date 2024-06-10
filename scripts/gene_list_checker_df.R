# Gene list checker df
protein.coding.genes2 <- protein.coding.genes %>%
  select(symbol, alias_symbol, prev_symbol) %>%
  rename(gene_symbol = symbol) %>%
  separate_rows(alias_symbol, sep = "\\|") %>%
  separate_rows(prev_symbol, sep = "\\|")

write.fst(protein.coding.genes2, "./data/processed/pcg.checker.df.fst")
