# Setup

if (!require(here, quietly=TRUE)) {
  if (!require(devtools, quietly=TRUE)) {
    install.packages("devtools")
    library(devtools)
  }
  devtools::install_github("krlmlr/here")
  library(here)
}

here::here()

source(here::here("phylify/ontology-explorer.r"))

# load the ontology
# TODO: add option to generate COdat.Rdata
load(here::here("phylify/COdat.RData"))

# Load the fantom functions
source(here::here("fantom_import/fantom_main.R"))

# Search for FANTOM.5.ONTOLOGY.ID
# say we want to search for
sourcecell <- "eye"
target <- "fibroblast"

# to get the IDs
fantomSearch(sourcecell)
fantomSearch(target)

# From the search, we can choose
sourceFF <- "FF:10272-104E2"
targetFF <- "FF:11268-116G8"

# Create a request vector
requestVector <- c(sourceFF, targetFF)

# Download raw expression for the FFs (result stored in fantomResults)
fantomOntology(requestVector)

# Filter out with raw expressions less than (results stored in fantomCounts)
# Normalize gene names as well
fantomSummarize(5)

# load up contrast
source(here::here("contrast/contrast.R"))

# Exclude first col
counts <- fantomCounts[ , -1]
rownames(counts) <- fantomCounts$short_description
colnames(counts) <- c("eye", "fibroblast")

diffExp <- contrast(counts)

## Section 3.4.1

# Load up normalized WEAVE

WEAVE_FILE <- here::here("WEAVE/curatedOutput.RData")

if (!file.exists(WEAVE_FILE)) {
  source(here::here("normalize/normalizeWeave.R"))
} else {
  load(WEAVE_FILE)
}

# Get high confidence interactions only
high_conf_interactions = src[src$combined_score > 700,]

STRGRAPH <- graph_from_data_frame(high_conf_interactions, directed = FALSE)

# Load up TRRUST
source(here::here("TRRUST_network/TRRUST_network.R"))

trrust <- loadTRRUST()
trrust <- fixColumns(trrust)

setMode(1)
genes <- trrust[,1]

trrust_nodes <- getNodes(trrust)
trrust_edges <- getEdges (trrust)
trrust_nodes <- getWeights(trrust_nodes, trrust_edges)

TRRUST_GRNGRAPH <- graph_from_data_frame(trrust_edges, directed = TRUE)

# Load up REGNET. Graph is stored as REGNETGRAPH
source(here::here("REGNET/REGNET.R"))

# WEAVE (combine) all the subgraphs together

# WEAVE-all contains this function
getTFSubgraph <- function(TF, order=1, GRAPH=STRGRAPH) {
  
  return(make_ego_graph(GRAPH, order, TF))
}


WEAVE_Path = here::here("WEAVE/WEAVE.RDATA")
save(STRGRAPH, REGNETGRAPH, TRRUST_GRNGRAPH, getTFSubgraph, file=WEAVE_Path)

SubgraphList <- getTFSubgraph("MYC", 2, REGNETGRAPH)	

## 3.5
# Need to find out current_gene
Lrn <- shortest.paths(SubgraphList, ,"MYC")