---
title: Design of the Experiment
contributors: [Abhishek Agarwal,Ziad Al-Bkhetan]
description: Guide to prepare design.csv for all three working cases in dcHiChIP.
toc: true
type: guides
---

This section provides a detailed explanation of the `design.csv` input format for the three working modes of **dcHiChIP**.

---

## 🧪 Case 1: Only HiChIP raw FASTQ files (no ChIP-seq input)

> **Use this format when you do not have ChIP-seq data (neither raw nor peak files).**

|id  | sample | hichip_r1                    | hichip_r2                    | chipseq_r1 | chipseq_r2 | narrowpeak |
|----|--------|------------------------------|------------------------------|------------|------------|------------|
|S-1 | S1     | ./hichip/SAMPLE1_R1.fastq.gz | ./hichip/SAMPLE1_R2.fastq.gz |            |            |            |
|S-2 | S1     | ./hichip/SAMPLE2_R1.fastq.gz | ./hichip/SAMPLE2_R2.fastq.gz |            |            |            |


**📝 Notes:**
- The last column (`chipseq`) & (`narrowpeak`) must be left empty.
- Pseudo-ChIP-seq peaks will be automatically generated from the HiChIP data.

---

## 📊 Case 2: HiChIP FASTQ + processed ChIP-seq peak files (narrowPeak)

> **Use this format when you have pre-processed ChIP-seq peak files in narrowPeak format.**

|id  | sample | hichip_r1                    | hichip_r2                    | chipseq_r1 | chipseq_r2 | narrowpeak                |
|----|--------|------------------------------|------------------------------|------------|------------|---------------------------|
|S-1 | S1     | ./hichip/SAMPLE1_R1.fastq.gz | ./hichip/SAMPLE1_R2.fastq.gz |            |            | ./chip/SAMPLE1.narrowpeak |
|S-2 | S1     | ./hichip/SAMPLE2_R1.fastq.gz | ./hichip/SAMPLE2_R2.fastq.gz |            |            | ./chip/SAMPLE2.narrowpeak |    

**📝 Notes:**
- Ensure peak files are in **BED6+4** format.
- Chromosome names **must** follow the `chrX` naming convention (e.g., `chr1`, `chr21`).

---

## 🔬 Case 3: HiChIP FASTQ + raw ChIP-seq FASTQ files (needs peak calling)

> **Use this format when you have raw ChIP-seq data and wish to process it within the pipeline.**

|id  | sample | hichip_r1                    | hichip_r2                    | chipseq_r1                    | chipseq_r2                    | narrowpeak |
|----|--------|------------------------------|------------------------------|-------------------------------|-------------------------------|------------|
|S-1 | S1     | ./hichip/SAMPLE1_R1.fastq.gz | ./hichip/SAMPLE1_R2.fastq.gz | ./chipseq/SAMPLE1_R1.fastq.gz | ./chipseq/SAMPLE1_R2.fastq.gz |            |
|S-2 | S1     | ./hichip/SAMPLE2_R1.fastq.gz | ./hichip/SAMPLE2_R2.fastq.gz | ./chipseq/SAMPLE2_R1.fastq.gz | ./chipseq/SAMPLE2_R1.fastq.gz |            |

**📝 Notes:**
- ChIP-seq input files must correspond to the same sample and replicate scheme.
- Peak calling will be automatically performed using MACS3.

---

Each mode triggers a distinct branch of the pipeline, ensuring flexible integration based on the type and availability of ChIP-seq data.

