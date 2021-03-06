# phylify

## Design

The overall design of phylify is to create a module which satisfies the following:

- easy OBO parsing without having to look up an OBO parser
- quick ontology summary
- ontology as an igraph object
- easy to customize interative visualization (implemented with visNetwork atm)
- helper functions to ease filtering ontology data
- example workflow to illustrate API use
- concise `COdat` creation script for VCS, Ontoscope purposes
- implements helper methods for filtering by various heuristics, such as
  - comparing to Mogrify
  - using human and mouse samples tables with cell lines, primary cells, time courses, etc
  - `filterByBad` and `filterByGood`

Summarily:
- produce COdat while maintaining an API that can be picked up by someone else
to experiment with and create and alternate COdat

## COdat

See [workflow](https://github.com/hyginn/Ontoscope/blob/master/phylify/workflow.r)
for how to make your own `COdat.RData`. The current version takes all `FF:X` IDs
from the Mogrify website, plus all ~~~"cell lines"+"primary cells" categories
in humans from the set of~~~ `FF:A-B` IDs for which we can pull samples to run
through DESeq. (Human cell lines + primary cells was too restrictive, and no harm
having a little more)

## OBO Parsing Heuristics

See [here](https://github.com/hyginn/Ontoscope/blob/master/phylify/obo-heuristic.md).

## ontology-explorer

#### `getOBO(<string>)`

Returns an `OBOCollection`.

#### `summarizeOBO(<OBOCollection>, head=<logical>, n=6L)`

Returns a list with `numTerms`, `numEdges`, `termIDs`, `termTypes`, and `propertyTypes`.
Use `head=TRUE` to take a glance.

#### `makeVisNetwork(*params)`

Params:
- `smooth` <logical> whether or not to use curved edges. warning: slower!
- `useLabel` <logical> whether or not to take label attribute. use false if your nodes don't have a label
- `cluster` <logical> whether or not to cluster (with clusterAlg)
- `clusterAlg` <FUN> the igraph function to cluster by
- `clusterAsUndirected` <logical> whether to cast graph to cluster to be undirected. some algs. only work on undirected
- `customGroups` <logical> whether or not to use your own custom $group vertex attributes
- `hierarchicalLayout` <logical> whether or not to let vis.js compute a hierarchical layout
- `levelSeparation` <number> for vis.js hierarchical layout
- `direction` <string> direction of tree for hierarchical layout
- `igraphLayout` <logical> whether or not to compure layout on igraph side vs vis.js side
- `layout` <string> the name of the igraph layout function to use

## visNetwork

To make visNetwork full screen, create a bookmark on your bookmarks bar and put
this in the URL:

```js
javascript:(() => { const styles=document.createElement('style');  styles.innerHTML = 'body { padding: 0 !important; } html, body, #htmlwidget_container, .visNetwork { height: 100% !important; }'; document.body.appendChild(styles) })()
```

Then you can click it whenever you are viewing a visNetwork.

The code is a *self-executing anonymous function*. Here it is formatted properly:
```js
(() => {
  const styles = document.createElement('style');
  styles.innerHTML = 'body { padding: 0 !important; } html, body, #htmlwidget_container, .visNetwork { height: 100% !important; }';
  document.body.appendChild(styles)
})()
```

Quite simple and clean, just a little CSS, though it seems visNetwork authors want you to specifiy specific height..

## Plots

Some visNetwork plots as images [here](https://github.com/hyginn/Ontoscope/blob/master/phylify/plots).
TODO: export html plots.

## ontoCAT

ontoCAT provides methods that can be useful for filtering the ontology.

**Note!** ontoCAT converts `:` into `_` for term IDs. For example,
`FF:11436-118H5` becomes `FF_11436-118H5`. There already apparently unique
terms where one is `_` and the other `:` (only a few), so I wonder how this
conversion might mess things up then..

Quick API walkthrough:

```r
> fantomCAT <- getOntology(normalizePath("ff-phase2-140729.obo"))

> head(getAllTerms(fantomCAT)) # no guaranteed order, also getAllTermIds
[[1]]
FF_0000045: human plasmacytoid dendritic cell sample

[[2]]
CL_0002597: smooth muscle cell of bladder

[[3]]
FF_11436-118H5: Smooth Muscle Cells - Internal Thoracic Artery, donor3

[[4]]
FF_13540-145H3: cerebellar granule cells, embryo E18, biol_rep3 (E13R3)

[[5]]
FF_11793-124C2: CD4+CD25+CD45RA+ naive regulatory T cells expanded, donor1

[[6]]
DOID_0060058: lymphoma

> getTermById(fantomCAT, "FF:11436-118H5") # can still give it : version
FF_11436-118H5: Smooth Muscle Cells - Internal Thoracic Artery, donor3

> getTermNameById(fantomCAT, "FF:11436-118H5")
[1] "Smooth Muscle Cells - Internal Thoracic Artery, donor3"

> getTermParentsById(fantomCAT, "FF_11436-118H5")
[[1]]
EFO_0002091: biological replicate

[[2]]
FF_0000173: human smooth muscle cell of the internal thoracic artery sample

> getTermChildrenById(fantomCAT, "FF_0000173")
[[1]]
FF_11364-117I5: Smooth Muscle Cells - Internal Thoracic Artery, donor2

[[2]]
FF_11287-116I9: Smooth Muscle Cells - Internal Thoracic Artery, donor1

[[3]]
FF_11436-118H5: Smooth Muscle Cells - Internal Thoracic Artery, donor3

> getAllTermParentsById(fantomCAT, "FF_11436-118H5")
[[1]]
FF_0000102: sample by type

[[2]]
owl:Role: role

[[3]]
FF_0000001: sample

[[4]]
FF_0000210: human sample

[[5]]
FF_0000167: smooth muscle cell sample

[[6]]
FF_0000002: in vivo cell sample

[[7]]
owl:SpecificallyDependentContinuant: material property

[[8]]
FF_0000173: human smooth muscle cell of the internal thoracic artery sample

[[9]]
FF_0000101: sample by species

[[10]]
EFO_0000001: experimental factor

[[11]]
FF_0000175: human smooth muscle cell of subclavian artery sample

[[12]]
EFO_0002091: biological replicate

[[13]]
EFO_0000683: replicate

> length(getAllTermChildrenById(fantomCAT, "EFO_0002091"))
[1] 2255

> showHierarchyDownToTermById(fantomCAT, "EFO_0002091")
BFO_0000040 material entity
DOID_4 disease
Disposition Disposition
FF_0000001 sample
GO_0005575 cellular_component
MaterialEntity MaterialEntity
NCBITaxon_1 root
ProcessualEntity ProcessualEntity
Role Role
SpecificallyDependentContinuant SpecificallyDependentContinuant
UBERON_0001062 anatomical entity
EFO_0000001 experimental factor
-owl:MaterialEntity material entity
-IAO_0000030 information entity
-owl:ProcessualEntity process
-owl:SpecificallyDependentContinuant material property
--owl:Disposition disposition
--owl:Role role
---EFO_0001461 control
---EFO_0000683 replicate
----EFO_0002090 technical replicate
----EFO_0002091 biological replicate

> getRootIds(fantomCAT)
 [1] "BFO_0000040"                     "DOID_4"                          "Disposition"                     "EFO_0000001"                    
 [5] "FF_0000001"                      "GO_0005575"                      "MaterialEntity"                  "NCBITaxon_1"                    
 [9] "ProcessualEntity"                "Role"                            "SpecificallyDependentContinuant" "UBERON_0001062"
```
