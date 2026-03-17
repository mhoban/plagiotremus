# Origins and relationships of endemic Hawaiian saber-toothed blenny genus *Plagiotremus* (Blennioidei: Blenniidae)
This repository contains the raw manuscript and associated code and data for the paper:

Origins and relationships of endemic Hawaiian saber-toothed blennies in the genus *Plagiotremus* (Blennioidei: Blenniidae)

# Contents

## General

  * `manuscript.Rmd` &mdash; primary manuscript Rmarkdown file
  * `Makefile` &mdash; makefile to build docx output
  * `plagiotremus.Rproj` &mdash; RStudio project
  * `renv.lock` &mdash; `renv` lockfile, tracks `R` package dependencies

## Dependencies

  * R package dependencies can be installed by opening an R console in the repository root directory and running:
    ```r
    renv::restore()
    ```
  * Other external dependencies
    * [pandoc](https://pandoc.org/)
    * [pandoc-crossref](https://lierdakil.github.io/pandoc-crossref/)

## Code

  * `spp.R` &mdash; `R6` class to build taxonomic citations with authority via WoRMS
  * `tree_alignments.R` &mdash; `R` script to generate alignments for phylogenetic analyses

## Data ([data/](data/))

  * `samples.csv` &mdash; sample metadata
  * `meristics.csv` &mdash; morphological data for *P. goslinei/tapeinosoma/azaelus* assembled from literature
  * `plagiotremus-coi.fasta` &mdash; COI sequences for *Plagiotremus* spp.
  * `plagiotremus-cytb.fasta` &mdash; cytb sequences for *P. goslinei/tapeinosoma*

## Trees ([data/trees](data/trees))

  * `coi_c*.nex.con.tre` &mdash; COI MrBayes consensus trees (c00 = all sites, c99 = 99%, c50 = 50%)
  * `*_clade.tre` &mdash; clades exported for convenience from `coi_c50.nex.con.tre`, used for plotting
  * `tapeinosoma_biogeography.tre` &mdash; BEAST tree of *P. goslinei/tapeinosoma/azaleus*

### Phylogenetic analysis configurations ([data/trees/setup](data/trees/setup/))
  * `coi_c**.nex` &mdash; NEXUS files to configure MrBayes runs (c00 = all sites, c99 = 99%, c50 = 50%)
  * `tapeinosoma_biogeography.xml` &mdash; BEAST configuration file for ancestral reconstruction analysis
