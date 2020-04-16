
# COVID-19 Genotyping Tool (CGT)

### Please note that the Covid-19 Genotyping Tool backend and UI are under continuous development.

## Overview

[The Covid-19 Genotyping
Tool](https://hsmaan.shinyapps.io/CovidGenotyper/) (CGT) is an R-Shiny
based web application that allows researchers to upload fasta sequences
of Covid-19 viral genomes and compare with public sequence data
available on [GISAID](https://www.gisaid.org/). Genomic distance is
visualized using manifold projection and network analysis, and genotype
information with respective to high-prevalence SNPs is determined.

## Details and methodology

#### UI and visualizations

The CGT application was developed using the `shiny` R package and
framework. Visualizations are generated using the `ggplot2`,
`ggnetwork`, and `plotly` R packages. User uploaded fasta files are
considered input, and the visualizations reactively adjust to user data.
To expedite the loading of plots using public data, the alignment, DNA
distance, and plot data fetching steps are all pre-processed for GISAID
public sequence data. Once users upload in-house sequencing data, these
steps are re-performed with the concatenation of user and public data.

#### Sequence and metadata retrieval

Processed fasta files of Covid-19 viral genome sequence are retrieved
from the [GISAID](https://www.gisaid.org/) EpiCoV database, which is a
public database for sharing of viral genome sequence data. Metadata for
GISAID viral genomes are obtained from [nextstrain’s ncov
build](https://github.com/nextstrain/ncov/blob/master/data/metadata.tsv).
Viral genome data and metadata are updated on a weekly basis.

#### Genome sequence alignment

GISAID sequences are subset for those that have metadata from
nextstrain. Public sequencing data is pre-aligned before being uploaded
to the server. Fasta sequences are read and written using the
`Biostrings` package. Gap removal and multiple-sequence alignment is
performed using `DECIPHER`. Post alignment processing is done using
`ape`. User uploaded fasta sequences are processed similarly, with the
exception of complete alignment - the user sequence is aligned to the
pre-aligned public data profile using `AlignProfiles` from `DECIPHER`.

#### DNA distance

For both pre-aligned and profile aligned data, DNA distance is
determined using `ape` and the Kimura-80 model of nucleotide
substitution. Currently only Kimura-80 is supported, but integrating
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

#### Network analysis and minimum spanning tree (MST)

Network representations of DNA distance are visualized using the
`igraph` package. DNA distance matrices are utilized to create graph
representations. Connectivity is determined by employing the minimum
spanning tree, which aims to determine a path that connects all vertices
while minimizing distance. This method, along with the graphopt
force-directed layout aims to visualize connectivity of viral genomes
from affected individuals, potential paths of transmission based on
connectivity, and large network hubs which may be associated with
outbreak epicenters.

#### Single-nucleotide polymorphisms (SNP)

Genotype profiles of viral genomes are determined using high prevalence
non-synonymous SNPs (based on minor allele frequency) within structural
protein (E, M, N, S) genome regions in the public sequencing data. SNPs
were called using *snp-sites*, annotated using *snpEff*, and frequency
analysis was performed using *vcftools*. The top 9 most frequent
non-synonymous (missense, nonsense) SNPs within structural protein
genome regions are presented.

## References

  - Hadfield J, Megill C, Bell SM, Huddleston J, Potter B, Callender C,
    et al. NextStrain: Real-time tracking of pathogen evolution.
    Bioinformatics. 2018;34(23):4121–3.

  - Elbe S, Buckland-Merrett G. Data, disease and diplomacy: GISAID’s
    innovative contribution to global health. Glob Challenges.
    2017;1(1):33–46.

  - Chang W, Cheng J, Allaire JJ, Xie Y, McPherson J. Shiny: Web
    Application Framework for R. R package version 1.4.0.2. 2020.
    Available from: <https://CRAN.R-project.org/package=shiny>

  - Wickham H. ggplot2: Elegant Graphics for Data Analysis.
    Springer-Verlag New York, 2016.

  - Briatte F. ggnetwork: Geometries to Plot Networks with ‘ggplot2’. R
    package version 0.5.8. 2020.

  - C Sievert. Interactive Web-Based Data Visualization with R, plotly,
    and shiny. Chapman and Hall/CRC Florida, 2020.

  - McInnes L, Healy J, Melville J. UMAP: Uniform Manifold Approximation
    and Projection for Dimension Reduction. arXiv \[Internet\]. 2018;
    Available from: <http://arxiv.org/abs/1802.03426>

  - Mamun A, Rajasekaran S. An efficient Minimum Spanning Tree
    algorithm. In: 2016 IEEE Symposium on Computers and Communication
    (ISCC). 2016. p. 1047–52.

  - Benson DA, Karsch-Mizrachi I, Lipman DJ, Ostell J, Wheeler DL.
    GenBank. Nucleic Acids Res. 2007;35(SUPPL. 1):21–5.

  - Pagès H, Aboyoun P, Gentleman R, DebRoy S. Biostrings: Efficient
    manipulation of biological strings. R package version 2.54.0. 2019.

  - Wright ES. Using DECIPHER v2.0 to analyze big biological sequence
    data in R. R J. 2016;8(1):352–9.

  - Paradis E, Claude J, Strimmer K. APE: Analyses of phylogenetics and
    evolution in R language. Bioinformatics. 2004;20(2):289–90.

  - Melville J. uwot: The Uniform Manifold Approximation and Projection
    (UMAP). Method for Dimensionality Reduction. R package version
    0.1.8. 2020.

  - Csardi G, Nepusz T. The igraph software package for complex network
    research. InterJournal, Complex Systems. 1695. 2006.
    <http://igraph.org>

  - Page AJ, Taylor B, Delaney AJ, Soares J, Seemann T, Keane JA, et
    al. SNP-sites: rapid efficient extraction of SNPs from multi-FASTA
    alignments. Microb genomics. 2016;2(4):e000056.

  - Cingolani P, Platts A, Wang LL, Coon M, Nguyen T, Wang L, et al. A
    program for annotating and predicting the effects of single
    nucleotide polymorphisms, SnpEff: SNPs in the genome of Drosophila
    melanogaster strain w1118; iso-2; iso-3. Fly. 2012.

  - Danecek P, Auton A, Abecasis G, Albers CA, Banks E, DePristo MA, et
    al. The variant call format and VCFtools. Bioinformatics.
    2011;27(15):2156–8.

## License

[GNU General Public
License 3.0](https://www.gnu.org/licenses/gpl-3.0.en.html)
