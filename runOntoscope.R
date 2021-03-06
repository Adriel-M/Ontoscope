# runOntoscope.R
#
# Purpose:  Run an "Ontoscope" query, and analysis based on
#                and extending the mogrify workflow of
#                Rackham et al. (2016).
#
# V 0.00001
# Date:     March 29. 2016
# Author:   Boris Steipe and BCB420 class
#
# V 0.00001 Very first crude code to integrate modules
#
# ==========================================================


if (!require(here, quietly=TRUE)) {
  if (!require(devtools, quietly=TRUE)) {
    install.packages("devtools")
    library(devtools)
  }
  devtools::install_github("krlmlr/here")
  library(here)
}

here::here()

# ==== INITS =============================
source(here::here("phylify/ontology-explorer.r"))



# ========================================




# load the ontology
load(here::here("phylify/COdat.RData"))

mogList <- getMogrifyIDs()

targetParent <- "FF:0000592"


# use gather module to get a list of background Fantom IDs
x1 <- neighborhood(COdat, order=1, nodes= targetParent, mode="out")
x2 <- neighborhood(COdat, order=2, nodes= targetParent, mode="out")
x3 <- neighborhood(COdat, order=3, nodes= targetParent, mode="out")

parent <- "FF:0101581"
bck <- neighborhood(COdat, order=2, nodes=parent, mode="in")

# added
source(here::here("phylify/createCOdat.R"))
# source Fantom functions
source(here::here("fantom_import/fantom_main.R"))
fantom

x <- fantomSearch("cardiac")
cell <- levels(x$FANTOM.5.Ontology.ID)[400]
bck <- levels(x$FANTOM.5.Ontology.ID)[1:2]

request <- c(cell, bck)

fantomOntology(request)

fantomSummarize(5)  #parameter is min. expression count threshold

# load up CONTRAST

source(here::here("contrast/contrast.R"))

# NOTE: change format!
counts <- fantomCounts[ , -1]
rownames(counts) <- fantomCounts$short_description
colnames(counts) <- c("fib", "card1", "card2")

diffExp <- contrast(counts)


# Build a network

# Normalize the STRING file by executing functions in ./WEAVE/normalizeWeave.R
# If this has been run before, you only need to reload the file.
# load(file = STRINGnew)

# use only high-confidence edges (score > 900)
# TODO: find out which src
# in normalize/normalizeWEAVE.R
# Finish save to file
tmp <- src
src <- src[src$combined_score > 900, ]
nrow(src)

# build the WEAVE network using WEAVE-STRING.R fucntions
# STRGRAPH <- graph_from_data_frame(src, directed = FALSE)



load(here::here("MARA.RData"))
# needs code fpath that builds a network for the cell-types we are working with.
# ....


# using tools from TRRUST_network to load TRRUST GRN
#works if ontoscope is the working directory
#This is the entire TRRUST network with no subsetting
#Directed from Transcription Factor to Gene

source(here::here("TRRUST_network/TRRUST_network.R"))
trrust <- loadTRRUST()
setMode(1)
genes <- trrust[,1]
trrust_nodes <- getNodes(trrust)
trrust_edges <- getEdges (trrust)
trrust_nodes <- getWeights(trrust_nodes,trrust_edges)

TRRUST_GRNGRAPH <- graph_from_data_frame(trrust_edges, directed = TRUE)



# [END]