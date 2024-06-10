## Symbol/ID checker
# Check for input genes which are aliases or previous symbols and Return a list of up-to-date gene symbols & HGNC IDs

# Library
library(readxl)
library(tools)

# Read file
readFile <- function(file_path, pcg) {
  # Get file extension
  file_extension <- tools::file_ext(file_path)
  all_genes <- c(pcg$symbol, pcg$alias_symbol, pcg$prev_symbol)
  
  # Read file based on file extension
  if (file_extension == "fst") {
    # Read fst or RData file
    data <- read.fst(file_path)
  } else if (file_extension == "csv") {
    # Read csv file
    data <- read.csv(file_path)
  } else if (file_extension == "xlsx") {
    # Read Excel file
    data <- readxl::read_excel(file_path)
  } else if (file_extension == "tsv") {
    # Read tab-separated values (TSV) file
    data <- read.delim(file_path, sep = "\t")
  } else if (file_extension == "txt") {
    # Read text (TXT) file
    data <- read.table(file_path) # Adjust arguments based on your TXT format
  } else if (file_extension %in% c("rda", "rds", "RData")) {
    # Read RDS (R Data Serialization) file
    data <- readRDS(file_path)
  } else {
    # Unsupported file type
    stop("Unsupported file type.")
  }
  
  data <- data %>%
    pull(1)  %>%
    trimws()
  
  return(data)
}

checkInputs2 <- function(gene_list, pcg) {
  
  # Read list
  # gene_list <- readFile(input_gene_list_path, pcg)
  
  official_genes <- pcg$symbol
  alias_genes <- pcg$alias_symbol
  prev_genes <- pcg$prev_symbol
  all_genes <- c(pcg$symbol, pcg$alias_symbol, pcg$prev_symbol)
  
  official_genes2 <- intersect(official_genes, gene_list)
  alias_genes2 <- intersect(alias_genes, gene_list)
  # Only retain alias gene to convert (genes can be both alias & official/current)
  alias_genes2 <- setdiff(alias_genes2, official_genes2)
  prev_genes2 <- intersect(prev_genes, gene_list)
  # Only retain alias gene to convert (genes can be both prev & official/current) Another way to do this is only intersect by genes not captured by official_genes2
  prev_genes2 <- setdiff(prev_genes2, official_genes2)
  no_match <- setdiff(gene_list, all_genes)
  
  prev2official_tbl <- pcg %>%
    filter(prev_symbol %in% prev_genes2) %>%
    dplyr::select(symbol, prev_symbol) %>%
    distinct()
    
  prev2official2 <- pcg %>%
    filter(prev_symbol %in% prev_genes2) %>%
    pull(symbol) %>%
    unique()
  
  alias2official_tbl <- pcg %>%
    filter(alias_symbol %in% alias_genes2) %>%
    dplyr::select(symbol, alias_symbol) %>%
    distinct()
  
  alias2official2 <- pcg %>%
    filter(alias_symbol %in% alias_genes2) %>%
    pull(symbol) %>%
    unique()
  
  returned_list <- unique(c(official_genes2, prev2official2, alias2official2))
  
  print(returned_list)
  print(prev2official_tbl)
  print(alias2official_tbl)
  print(no_match)
  
  return(
    list(
      "gene_list" = returned_list,
      "prev2official" = prev2official_tbl,
      "alias2official" = alias2official_tbl,
      "no_match" = no_match
    )
  )
}


