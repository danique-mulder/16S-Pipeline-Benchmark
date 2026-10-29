# When Pipelines Matter (and When They Do Not)

Analysis code for:

> Mulder D, Natasha EE, Bastiaanssen T, Gosepa G, Szopinska-Tokov J, Poelmans G, Bloemendaal M, Homberg JR, Arias Vasquez A. *When Pipelines Matter (and When They Do Not): Benchmarking 16S rRNA Gene Bioinformatic Pipelines for Human Gut Microbiome Analysis.* (in preparation) · 

The scripts are shared for transparency, to show how each step in the benchmark was carried out. They were written for the Snellius cluster with hard-coded paths and settings.

Code archive: 10.5281/zenodo.23767633

Data archive: 10.5281/zenodo.23597756

## Workflow

```mermaid
flowchart LR
    prep["Data preparation<br/>mock · clinical · simulated"]
    qc["Quality control<br/>FastQC / MultiQC"]
    pip["15 pipelines<br/>default settings"]
    cls["Harmonized classification<br/>Q2-NBC · IDTAXA<br/>SILVA · Greengenes2"]
    pic["Functional prediction<br/>PICRUSt2"]
    tse["Tree summarized<br/>experiments (mia)"]
    bench["Benchmark metrics<br/>detection · agreement"]
    rep["Reports and figures<br/>R Markdown"]

    prep --> qc --> pip
    pip --> cls --> tse
    pip --> pic --> tse
    tse --> bench --> rep
```

## Structure

| Folder | Contents |
|---|---|
| `data/` | Processing mock communities (Mock 5, Mock 12): download, demultiplexing, primers, re-multiplexing and SILVA v138 reference taxonomy. Clinical cohorts (GEMSIM, MINDSet, Parkinson's) |
| `simulation/` | In silico communities, implanted group differences, chimeras and sequencing errors |
| `qc/` | FastQC/MultiQC |
| `pipelines/` | One script per pipeline; `runallpipelines.sh` submitted them all with the same settings |
| `classify/` | Harmonized classification (Q2-NBC, IDTAXA) |
| `picrust/` | PICRUSt2 functional prediction |
| `helper/` | Primers, re-multiplexing, creating combined `TreeSummarizedExperiment` objects |
| `benchmark/` | Benchmark metrics; statistical analyses and figures (`pipelines.Rmd`, `pipeline_plots.Rmd`) |

## Script dependencies

Solid arrows: one script calls another. Dashed arrows: a script reads files written by another.

### Running the pipelines

```mermaid
flowchart LR
    run[runallpipelines.sh]
    run --> AMPtk.sh & APSCALE.sh & DADA2.sh & DADAist2.sh & dadasnake.sh
    run --> eDNAFlow.sh & LotuS2.sh & MetaWorks.sh & Natrix.sh & nfcore_ampliseq.sh
    run --> NGTax2.sh & QIIME2.sh & QIIME2_deblur.sh & Tourmaline.sh & VSEARCH.sh
    DADA2.sh --> DADA2.R
    Natrix.sh --> mih[helper/make_illumina_headers.py]
    VSEARCH.sh --> silva2sintax.sh
```

### Taxonomic classification

```mermaid
flowchart LR
    classify.sh --> q2_naivebayes.sh
    classify.sh --> IDTAXA.sh
    IDTAXA.sh --> IDTAXA.R
    learnTaxa.sh --> learnTaxa.R
    learnTaxa.sh -.-> IDTAXA.sh
```

### Data preparation and QC

```mermaid
flowchart LR
    g[data/processgemsim.sh] --> qc[qc/qc.sh]
    p[data/processparkinson.sh] --> qc
    m[data/processmindset.sh] --> qc
```

### Mock communities

```mermaid
flowchart LR
    dl[data/downloadmock.sh] -.-> m5[data/process_mock5.sh]
    dl -.-> m12[data/process_mock12.sh]
    m5 --> ap[helper/addprimers.sh]
    m5 --> remux[helper/remux.sh]
    m12 --> ap
```

### Simulation

```mermaid
flowchart TB
    sim[simulate.sh]
    sim --> filter_meta.py & filter_features.py & filter_fasta.py & checkdup.py
    sim --> filter_blastn.py & subset_db.py & alternativetax.py & filter_uniquegenus.py
    sim --> merge_counts.py & mergetax.py & remap_asvs.py & simba.R
    sim --> generatereads.sh & picrust2_exp.sh & classifysims.sh
    sim --> qc[qc/qc.sh]
    sim --> ap[helper/addprimers.sh]
    classifysims.sh --> idt[classify/IDTAXA.R]
```

### Benchmarking

```mermaid
flowchart LR
    asv_precisionrecall.sh --> precision_NGTax_asv.py
    asv_precisionrecall.sh --> recall_NGTax_asv.py
    cor_pip.sh --> cor_pip.R
    permanova.sh --> permanova.R
    permanova_group.sh --> permanova_group.R
    tse_filter.sh --> tse_filter.R
    runrmd.sh --> pipelines.Rmd
```

### Building the summarized experiments

```mermaid
flowchart LR
    a[helper/asv_tse.sh] --> u[helper/unified_asv.sh]
    a --> r[helper/asv_tse.R]
    t[helper/trim_asv_tse.sh] --> tu[helper/trimmed_unified_asv.sh]
    t --> r
    e[helper/tse_exp.sh] --> er[helper/tse_exp.R]
```

### Inputs to the R Markdown reports

```mermaid
flowchart LR
    pc[picrust/picrust2.sh] -.-> pt[helper/picrust_tse.R]
    pt -.-> rmd[pipelines.Rmd]
    at[helper/asv_tse.sh] -.-> rmd
    cor_pip.sh -.-> rmd
    permanova_group.sh -.-> rmd
    tse_filter.sh -.-> rmd
    getseqlength.sh -.-> rmd
    readsin.sh -.-> rmd
    rc[data/reclassify_mock.sh] -.->|reference taxonomy| rmd
    permanova.sh -.-> rmd
    permanova.sh -.-> plots[pipeline_plots.Rmd]
    rmd -.-> plots
```

## Data availability

- Mock communities: [Mockrobiota](https://github.com/caporaso-lab/mockrobiota)
- Parkinson's cohort: NCBI SRA [PRJNA601994](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA601994) (dataset 2)
- Simulated datasets, pipeline outputs for all datasets (including the clinical cohorts) and benchmark metrics: [10.5281/zenodo.23597756](https://doi.org/10.5281/zenodo.23597756)
- MINDSet and GEMSIM raw sequences and clinical metadata: closed access in the Radboud Data Repository [link], available on request through [contact]
