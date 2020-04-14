
### Please note that the CovidGenotyper backend and UI are still under active development.

## Overview

[CovidGenotyper](https://hsmaan.shinyapps.io/CovidGenotyper/) is an
R-Shiny based web application that allows researchers to upload SARS-CoV-2 consensus genome sequences (fasta)
and compare with publically available genome sequence
data on [GISAID](https://www.gisaid.org/). Genomic distance is
visualized using manifold projection (UMAP) and minimum spanning tree (MST) network analysis, and genotypes based on highly prevalent non-synonymous single nucleotide polymorphisms (SNPs) within structural genes (E, N, S, M).

## Methodology

#### Sequence and metadata retrieval

Processed consensus fasta files of SARS-CoV-2 viral genome sequences are retrieved
from the [GISAID](https://www.gisaid.org/) EpiCoV database, a
public database for sharing of viral genome consensus sequence data. Metadata for
GISAID viral genomes are obtained from [nextstrain’s ncov
build](https://github.com/nextstrain/ncov/blob/master/data/metadata.tsv).
Viral genome data and metadata are updated on a daily basis.

#### Genome sequence alignment

GISAID sequences are filtered to only include sequence data accompanied by metadata from
nextstrain. Consensus sequences obtained from GISAID are pre-aligned before being uploaded
to the server. GISAID consensus fasta sequences are read and written using the
`Biostrings` package. Gap removal and multiple-sequence alignment is
performed using `DECIPHER`. Post alignment processing is done using
`ape`. User uploaded fasta sequences are processed similarly, with the
exception of performing a complete alignment - the user sequences are aligned to the
pre-aligned GISAID sequence alignment profile using `AlignProfiles` from `DECIPHER`.

#### DNA distance

For both pre-aligned and profile aligned sequence data, DNA distance is
determined using `ape` and the Kimura-80 model of nucleotide
substitutions. Currently only Kimura-80 is supported, but integrating
other evolutionary distance metrics will be part of a future release.

#### Uniform manifold projection and approximation (UMAP)

UMAP is performed to visualize DNA distances between public data and
user uploaded data. In this context, UMAP aims to cluster groups of
closely-linked viral genomes together based on DNA distance. The `uwot`
implementation of UMAP is utilized using the pre-computed DNA distance
matrix. Defaults for the `uwot` implementation are employed, with the
exception of the following parameters:

`init = "spectral"`<br/> `metric = "cosine"`<br/> `n_neighbors
= 50`<br/> `min_dist = 0.001`<br/> `spread = 30`<br/>
`local_connectivity = 10`<br/>

#### Network analysis and minimum spanning tree

Network representations of DNA distance are visualized using the
`igraph` package. DNA distance matrices are utilized to create graph
representations. Connectivity is determined by employing the minimum
spanning tree, which aims to determine a path that connects all vertices
while minimizing distance. This method, along with the graphopt
force-directed layout aims to visualize connectivity of viral genomes
from affected individuals, potential paths of transmission based on
connectivity, and large network hubs which may be associated with
outbreak epicenters.

#### Single-nucleotide polymorphisms

Genotype profiles of viral genomes are determined using high prevalence
non-synonymous SNPs (minor allele frequence \> 0.05) within structural genes (E, N, S, M) present in the GISAID consensus sequence data.
SNPs were called using *snp-sites* and linkage analysis was performed
using *vcftools*. Both coding and non-coding variants are considered in
the genotyping analysis, with the exception of the 3’ and 5’ UTRs. All
alleles of SNPs with a significant minor allele are visualized, and
compared to the alleles of the user uploaded viral genome.

## References

### To be completed

## License

[GNU General Public
License 3.0](https://www.gnu.org/licenses/gpl-3.0.en.html)
