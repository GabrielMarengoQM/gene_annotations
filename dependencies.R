## Dependencies for the project ##
packages <- c(
  # For raw data retrieval 
  'fst', 'data.table', 'readr', 'readxl', 'biomaRt', 'PANTHER.db', 'org.Hs.eg.db', 'GO.db', 'reactome.db',
  # For wrangling
  'dplyr', 'tidyr', 'stringr')
lapply(packages, library, character.only = TRUE)