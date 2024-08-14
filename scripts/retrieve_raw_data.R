## Call APIs and read local static files then store as FST files in ./data/raw/

## -----------------------------------------------------------------------------
# Data retrieved:
# 1 - Human protein coding genes (EBI)
# 2 - Mouse genes, viability and phenotypes (IMPC, MGI)
# 3 - Disease phenotypes, moi, genetic basis, lethality categories and tissue (OMIM, DDG2P)
# 4 - Gene constraint metrics from cell lines (hPSCs, DepMap)
# 5 - Gene constraint metrics from population sequencing (gnomAD, Shet, AlphaMissense, DOMINO, SCoNeS)
# 6 - Protein families (PantherDB)
# 7 - Reacome pathways (Reactome)
# 8 - Gene Ontology (GO)

# !!! Before running !!! 
# - Check sites and repos for newest release version
# - Ensure OMIM api key is provided
api_key <- "INSERT_YOUR_OWN_API_KEY_HERE"
## -----------------------------------------------------------------------------

## -----------------------------------------------------------------------------
# Data retrieval ----
# 1 - Human protein coding genes (EBI) ----
protein.coding.genes <- read.delim(
  "ftp://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/locus_types/gene_with_protein_product.txt"
)

# 2 - Mouse genes, viability and phenotypes (IMPC, MGI) ----
# IMPC viability, phenotypes, WoL & completeness
mouse.viability.impc <- fread(
  "http://ftp.ebi.ac.uk/pub/databases/impc/all-data-releases/release-20.1/results/viability.csv.gz"
)
mouse.phenotypes.impc <- fread(
  "http://ftp.ebi.ac.uk/pub/databases/impc/all-data-releases/release-20.1/results/genotype-phenotype-assertions-ALL.csv.gz"
)
wol.categories <- read.table("../../Downloads/sup_file_1.txt", header = TRUE)
mouse.procedure.completeness <- fread(
  "http://ftp.ebi.ac.uk/pub/databases/impc/all-data-releases/release-20.1/results/procedureCompletenessAndPhenotypeHits.csv.gz"
)
# MGI viability, phenotypes, protein coding genes & lethal terms
mouse.viability.mgi <- fread(
  "https://www.informatics.jax.org/downloads/reports/MGI_GenePheno.rpt"
)
mouse.phenotypes.mgi <- fread(
  header = FALSE,
  "https://www.informatics.jax.org/downloads/reports/VOC_MammalianPhenotype.rpt"
)
mouse.protein.coding.genes.mgi <- fread(
  "https://www.informatics.jax.org/downloads/reports/MRK_List2.rpt"
)
mouse.lethal.terms <- read.csv("/Users/gabrielm/Desktop/pipeline_raw_data/lethal_terms.csv")

# Orthologs from BioMart
# orthologs
human<- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

attributes<-  c("ensembl_gene_id", "external_gene_name",
                "mmusculus_homolog_ensembl_gene", 
                "mmusculus_homolog_associated_gene_name",
                "mmusculus_homolog_orthology_type",
                "mmusculus_homolog_orthology_confidence",
                "mmusculus_homolog_perc_id_r1")

listAttributes(human) %>% 
  filter(stringr::str_detect(name, "mmusculus_homolog_"))

orth.mouse<-  getBM(attributes, filters="with_mmusculus_homolog",
                    values=TRUE, mart = human, uniqueRows=TRUE)

listFilters(human)%>% head()
listFilters(human)%>% 
  filter(stringr::str_detect(name, "mmusculus"))

orth.mouse.high.conf <- orth.mouse %>%
  filter(mmusculus_homolog_orthology_confidence == 1) %>%
  dplyr::select(ensembl_gene_id, mmusculus_homolog_orthology_type) %>%
  distinct()

# orthologs from gentar
orthologs <- fread("https://www.gentar.org/orthology-api/api/ortholog/write_to_tsv_file")

# 3 - Disease phenotypes, moi, genetic basis, lethality categories and tissue (OMIM, DDG2P) ----
# OMIM data
url <- "https://data.omim.org/downloads/"
mim.Titles <- read_delim(
  paste0(
    url,
    api_key,
    "/mimTitles.txt"
  ),
  delim = "\t",
  col_names = TRUE, 
  skip = 2
)
morbidmap <- read_delim(
  paste0(
    url,
    api_key,
    "/morbidmap.txt"
  ),
  delim = "\t",
  col_names = TRUE, 
  skip = 3
)
genemap <- read_delim(
  paste0(
    url,
    api_key,
    "/genemap2.txt"
  ),
  delim = "\t",
  col_names = TRUE, 
  skip = 3
)
# Requires new token - go to url and copy paste
lethal.phenotypes <- fread("https://raw.github.qmul.ac.uk/whri-phenogenomics/catalogue_lethal_genes/master/catalogue_lethal_genes_app/data/omim_curation.tsv?token=INSERT_YOUR_TOKEN_HERE")
# DDG2P
dd.g2p <- fread(
  "http://ftp.ebi.ac.uk/pub/databases/gene2phenotype/28_07_2023/DDG2P_28_7_2023.csv.gz"
)

# 4 - Gene constraint metrics from cell lines (hPSCs, DepMap) ----
# 5 - Gene constraint metrics from population sequencing (gnomAD, Shet, AlphaMissense, DOMINO, SCoNeS) ----
raw_data_files <- c("mef_metrics.xlsx",
                    "laminin_metrics.xlsx",
                    "shet_rgcme.txt",
                    "shet_posterior.tsv",
                    "scones_domino.xlsx",
                    "alphamissense_metrics.csv",
                    "am_path_scores_and_classes.csv",
                    "gnomad.v4.0.constraint_metrics.tsv")
file_path <- "/Users/gabrielm/Desktop/december_package/inst/extdata/raw/constraint_metrics/"
# Depmap
models <- read.csv("/Users/gabrielm/Desktop/december_package/inst/extdata/raw/constraint_metrics/depmap/Model.csv")
cell0 <- read.csv("/Users/gabrielm/Desktop/december_package/inst/extdata/raw/constraint_metrics/depmap/CRISPRGeneEffect.csv")
# mef & laminin
mef <- read_excel(paste0(file_path, paste0("mef_laminin/", raw_data_files[1])))
laminin <- read_excel(paste0(file_path, "mef_laminin/", raw_data_files[2]))
# gnomad
gnomad_v4 <- read.delim(paste0(file_path, "gnomad_v4/", 'gnomad.v4.0.constraint_metrics.tsv'))
# human canonical transcripts
mart <- useEnsembl("ensembl", dataset="hsapiens_gene_ensembl")
BM.info <- getBM(attributes=c("ensembl_gene_id","ensembl_transcript_id","hgnc_symbol", "hgnc_id", "transcript_is_canonical"), mart=mart)
# shet rgcme
shet_rgcme <- read.delim(paste0(file_path, "shet_rgcme/", raw_data_files[3]))
# shet post
shet_post <- read.delim(paste0(file_path, "shet_post/", raw_data_files[4]))
# scones & domino
scones_domino <- read_xlsx(paste0(file_path, "scones_domino/", raw_data_files[5]))
# alphamissense
alpham2 <- read_tsv("../../Downloads/AlphaMissense_gene_hg38.tsv.gz", skip = 3)
# get hgnc ens transcript id mappings from biomart 
ensembl <- useMart("ENSEMBL_MART_ENSEMBL",
                   dataset="hsapiens_gene_ensembl")
values1 <- alpham2$ensembl_transcript_id
att1 <- c("ensembl_transcript_id", "hgnc_id")
filter1 <- c("ensembl_transcript_id")
searchResults <-getBM(att1,
                      filters= filter1,
                      values= values1, mart=ensembl,
                      checkFilters = FALSE, # dont think this is needed
                      verbose = TRUE)

# 6 - PANTHERdb
pthOrganisms(PANTHER.db) <- "HUMAN"

# define columns to get from PANTHERdb
cols_short <- c("UNIPROT", "CLASS_ID", "CLASS_TERM", "FAMILY_ID", "FAMILY_TERM", "SUBFAMILY_TERM")
uniprot_ids <- keys(PANTHER.db, keytype="UNIPROT")

# Get data from PANTHER.db and merge columns to contain only unique unpirot ids without losing rows data
panther_data_shortcols <- PANTHER.db::select(PANTHER.db, 
                                             keys=uniprot_ids, 
                                             columns=cols_short, 
                                             keytype="UNIPROT")

# 7 - Reactome ----
genes_universe <- protein.coding.genes %>%
  dplyr::select(symbol, entrez_id)

universe <- as.character(genes_universe$entrez_id)

# create object with mappings for entrez ids to path ids to path names
entrez_pathid <- as.data.frame(reactomeEXTID2PATHID)
pathid_name <- as.data.frame(reactomePATHID2NAME)

reactome_annotations <- entrez_pathid %>%
  full_join(pathid_name)

# 8 - GO ----
go_anot <- toTable(org.Hs.egGO) %>%
  dplyr::select(gene_id, go_id, Ontology) %>%
  distinct()

keys <- unique(go_anot$go_id)

go_terms <- AnnotationDbi::select(GO.db, keys = keys, keytype="GOID", columns=c("TERM")) %>%
  dplyr::rename(go_id = GOID,
                go_term = TERM)

go_anot_term_raw <- go_anot %>%
  left_join(go_terms)
## -----------------------------------------------------------------------------

## -----------------------------------------------------------------------------
# Save data to ./data/raw
write.fst(protein.coding.genes, "./data/raw/protein.coding.genes.fst")
write.fst(mouse.viability.impc, "./data/raw/mouse.viability.impc.fst")
write.fst(mouse.phenotypes.impc, "./data/raw/mouse.phenotypes.impc.fst")
write.fst(wol.categories, "./data/raw/wol.categories.fst")
write.fst(mouse.procedure.completeness, "./data/raw/mouse.phenotypes.completeness.fst")
write.fst(mouse.viability.mgi, "./data/raw/mouse.viability.mgi.fst")
write.fst(mouse.phenotypes.mgi, "./data/raw/mouse.phenotypes.mgi.fst")
write.fst(mouse.protein.coding.genes.mgi, "./data/raw/mouse.protein.coding.genes.mgi.fst")
write.fst(mouse.lethal.terms, "./data/raw/mouse.lethal.terms.fst")
write.fst(orth.mouse.high.conf, "./data/raw/human.mouse.orth.fst")
write.fst(orthologs, "./data/raw/orthologs.fst")
write.fst(mim.Titles, "./data/raw/mim.Titles.fst")
write.fst(morbidmap, "./data/raw/morbidmap.fst")
write.fst(genemap, "./data/raw/genemap.fst")
write.fst(lethal.phenotypes, "./data/raw/lethal.phenotypes.fst")
write.fst(dd.g2p, "./data/raw/dd.g2p.fst")
write.fst(models, "./data/raw/models.fst")
write.fst(cell0, "./data/raw/cell0.fst")
write.fst(mef, "./data/raw/mef.fst")
write.fst(laminin, "./data/raw/laminin.fst")
write.fst(gnomad_v4, "./data/raw/gnomad_v4.fst")
write.fst(BM.info, "./data/raw/canon_transcripts.fst")
write.fst(shet_rgcme, "./data/raw/shet_rgcme.fst")
write.fst(shet_post, "./data/raw/shet_post.fst")
write.fst(scones_domino, "./data/raw/scones_domino.fst")
write.fst(alpham2, "./data/raw/alpham2.fst")
write.fst(searchResults, "./data/raw/ens_transcript_hgnc_id_mappings.fst")
write.fst(panther_data_shortcols, './data/raw/panther.fst')
write.fst(reactome_annotations, './data/raw/reactome_annotations.fst')
write.fst(go_anot_term_raw, "./data/raw/go.fst")
