# MANGO: Multi-case Active oNtology-based GO organizer
MANGO is an R package for Gene Ontology (GO) Biological Process enrichment analysis that reduces redundancy in top-ranked results by restructuring enriched terms into ontology-guided term trees based on the GO DAG. It defines and filters active trees using coverage/consistency criteria to suppress structurally driven false positives arising from hierarchical dependencies. MANGO supports single- and multiple- case study designs by integrating enrichment outputs across conditions, providing scoring options and visualization utilities to summarize common and condition-specific biological processes. Optional features include ORA-style soft filtering and fold-change-aware weighting for tree and term prioritization.

![Alt text](./FIG/MANGO_F1.jpg "MANGO")
>- **Problem:** GO Biological Process enrichment results are often dominated by biologically similar top terms, and GO’s hierarchical DAG + gene sharing can create **structural false positives**, reducing interpretability and confidence.
>- **Solution (MANGO):** An analytical framework that reorganizes GO outputs into **ontology-informed term trees** and applies **structural filtering** based on **within-tree consistency**.
>- **Key idea 1 — Tree organizing:** MANGO groups similar terms into trees to reduce redundancy and compress results into interpretable units.
>- **Key idea 2 — Active-tree filtering:** MANGO suppresses isolated, single-term–driven signals using an **active-tree criterion**.
>- **Multiple-case support:** MANGO aligns results from multiple comparisons into a **shared tree structure** and classifies **common vs condition-specific** trees using **HWES** and relative signal strength across conditions.
>- **Parameter recommendations:** A parameter search provides **input-size–specific recommended settings**.
>- **Benchmarking:** Evaluated across **31 single-case analyses** against **ORA** and **GSEA** using **adjusted p-values**, **fraction of terms assigned**, and **rich factors** to assess interpretability.
>- **Validation:** Demonstrated reproduction of reported signals across public datasets spanning **knockout**, **cancer**, **differentiation**, **dose–response**, and **cohort** designs.
>- **Usability:** Provides **R-based visualization functions** for both single- and multiple-case settings to aid interpretation.
>- **Take-home:** MANGO reduces redundancy and structural false positives driven by the GO DAG, enabling **reproducible term-pattern summarization** for multi-comparison studies.

## Installation

### Linux / macOS (terminal)

```bash
## Using conda OR micromamba

micromamba create -n MANGO
micromamba activate MANGO

micromamba install -c conda-forge -c r -c bioconda -y \
jupyter_core jupyter_client jupyterlab_pygments jupyter_server r-irkernel jupyterlab r=4.3.1

micromamba install -c conda-forge -c r -c bioconda -y \
r-png r-data.table r-systemfonts r-gdtools r-ggforce r-ggiraph bioconductor-xvector bioconductor-sparsearray \
bioconductor-biostrings bioconductor-delayedarray bioconductor-summarizedexperiment bioconductor-annotationdbi \
bioconductor-go.db bioconductor-keggrest bioconductor-fgsea bioconductor-deseq2 bioconductor-gosemsim pigz \
bioconductor-dose bioconductor-enrichplot bioconductor-clusterprofiler bioconductor-ggtree \
bioconductor-org.mm.eg.db 

mkdir path_to_MANGOanalysis
wget -P path_to_MANGOanalysis https://github.com/user-attachments/files/25766902/MANGO_FORMAT.tar.gz
### ex) wget -P /home/RNA/gitMANGO/ https://github.com/user-attachments/files/25766902/MANGO_FORMAT.tar.gz

cd path_to_MANGOanalysis
pigz -dc -p 4 MANGO_FORMAT.tar.gz | tar -xf -

mv MANGO_FORMAT MANGO_project_name
### ex) mv MANGO_FORMAT MANGO_DIFF

cd MANGO_project_name/MANGO
pwd
### result of pwd is filepath we use in R

cd count
cp path_to_input_countMATRIX RAWcount.txt
### ex) cp RAWcount_DIFF.txt RAWcount.txt
```

### R


```r
# install.packages("remotes")
remotes::install_github("ERASMUSlab/MANGO")
```

## Documentation

Full documentation and tutorials:
https://erasmuslab.github.io/MANGO


## Citation

If you use MANGO in your research, please cite:
https://erasmuslab.github.io/MANGO/authors.html#citation
