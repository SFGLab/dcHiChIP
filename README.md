# dcHiChIP

**dcHiChIP** is a modular and scalable pipeline for the **multi-scale analysis of chromatin architecture from HiChIP data**, developed using Nextflow DSL2. It supports different experimental scenarios, handles a wide variety of input formats, and integrates state-of-the-art tools for loop calling, stripe detection, compartment analysis, and 3D genome modelling.

---

## Key Features

The pipeline includes the following modules:

- **Read Alignment & Filtering** (BWA, SAMtools, pairtools)
- **Loop Calling** (MAPS)
- **Stripe Calling** (gStripe)
- **Compartment & TAD Analysis** (cooltools)
- **Motif & Peak Annotation** (HOMER, Bioframe)
- **Coverage & Correlation Plots** (DeepTools)
- **3D Genome Modelling** (MultiMM)
- **Subcompartment Calling** (CALDER)

---

## User Guide

A complete and interactive documentation is available at:  
📘 **https://sfglab.github.io/dcHiChIP/**

---

## Example Use Cases

dcHiChIP supports three types of input design:

- **Case 1:** HiChIP-only input → pseudo-ChIP signal
- **Case 2:** HiChIP + peak files (.narrowPeak) from ChIP-seq
- **Case 3:** HiChIP + raw ChIP-seq (with inputs) → calls peaks internally

Example datasets and sample designs are available in the documentation.

---

## Input Requirements

To run dcHiChIP, you need:

- HiChIP FASTQ files (paired-end)
- Optionally:
  - ChIP-seq FASTQ files (with or without inputs), or
  - Pre-processed `.narrowPeak` files
- A design CSV describing samples
- BWA index for the reference genome
- Chromosome sizes file
- Genomic features file (e.g., gene annotation BED)

---

## Installation

The pipeline uses [Nextflow](https://www.nextflow.io/) and supports [nf-core](https://nf-co.re/) compatible modularity.

To run locally or on a cluster:

```bash
git clone https://github.com/SFGLab/dcHiChIP.git
cd dcHiChIP
nextflow run main.nf -profile standard --input design.csv --fasta /path/to/genome.fa ...
