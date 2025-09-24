---
title: Documentation and User Guide for dcHiChIP
contributors: [Abhishek Agarwal, Ziad Al-Bkhetan]
description:
toc: false
---
## Introduction

<!-- pipeline figure -->
<p align="center">
  <img src="{{ site.baseurl }}/assets/dcHiChIP_workflow.png" alt="dcHiChIP workflow" width="100%">
</p>

**dcHiChIP** is a modular and reproducible Nextflow-based pipeline designed for high-throughput analysis of chromatin architecture using HiChIP data. The pipeline supports a wide range of tools for mapping, loop calling, peak annotation, chromatin feature extraction, and 3D modeling.

The pipeline includes integrated modules for:

- Mapping & filtering of HiChIP/ChIP-seq reads using `BWA`, `SAMtools`, and `MACS3`
- Loop calling using `MAPS` and conversion to BEDPE format
- Feature extraction with `cooltools`, `HOMER`, `gStripe`, and `CADLER`
- Multi-resolution compartment and TAD analysis
- 3D modeling of genome architecture using `MultiMM`

<hr/>

## Working Test Cases

The pipeline supports three fully documented test cases with representative input formats and figure overviews.

👉 Visit [Working Test Cases](./pages/working_cases.md) for more.

<hr/>

## Design & Sample Sheet

Explore the experimental design principles and detailed sample sheet formats compatible with the pipeline.

👉 Visit [Design and Sample Sheet Guide](./pages/design_samplesheet.md) for more.

<hr/>

## Module Parameters

Each pipeline module supports customizable parameters such as resolution, CPUs, input formats, and modelling levels.

👉 Visit [Module Parameters](./pages/parameters.md) for full details.

<hr/>

## Contact Us

For help, questions, or feedback, feel free to reach out.

👉 Visit [Contact Page](./pages/contact_us.md)

<hr/>

## Citation

If you use dcHiChIP in your research or publication, please cite the following manuscript (in preparation):

> Agarwal, A., Al-Bkhetan, Z., Plewczyński, D. *dcHiChIP: A comprehensive Nextflow-based pipeline for Multi-Scale Analysis of Chromatin Architecture from HiChIP Data.*

