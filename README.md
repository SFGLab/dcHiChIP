# dcHiChIP

**dcHiChIP** is a modular and scalable pipeline for the **multi-scale analysis of chromatin architecture from HiChIP data**, developed using Nextflow DSL2.
It supports different experimental scenarios, handles a wide variety of input formats, and integrates state-of-the-art tools for loop calling, stripe detection, compartment analysis, and 3D genome modelling.

---

## Key Features

The pipeline includes the following modules:

- **Read Alignment & Filtering** (BWA, SAMtools, pairtools)
- **Peak Calling** (MACS3)
- **Loop Calling** (MAPS)
- **Stripe Calling** (gStripe)
- **Compartment & TAD Analysis** (cooltools)
- **Motif & Peak Annotation** (HOMER, Bioframe)
- **Coverage & Correlation Plots** (DeepTools)
- **3D Genome Modelling** (MultiMM)
- **Subcompartment Calling** (CALDER)

---

## Documentation

A complete and interactive documentation is available at:
[dcHiChIP Documentation](https://sfglab.github.io/dcHiChIP/)

## Example Use Cases

dcHiChIP supports three types of input design:

- **Case 1:** HiChIP-only input → pseudo-ChIP signal
- **Case 2:** HiChIP + peak files (.narrowPeak) from ChIP-seq
- **Case 3:** HiChIP + raw ChIP-seq (with inputs) → calls peaks internally

Example datasets and sample designs are provided in the "Working Cases" section of [documentation](https://sfglab.github.io/dcHiChIP/working_cases).

---

## Input Requirements

To run dcHiChIP, you will need:

- Paired-end **HiChIP FASTQ files**
- Optionally:
  - ChIP-seq FASTQ files (with/without inputs), or
  - Pre-processed `.narrowPeak` files
- A **design CSV** describing samples
- **BWA index** for the reference genome
- **Chromosome sizes** file
- **Genomic features** file (e.g., BED gene annotation)

---

## Quick Start

1. **Install Nextflow (>=22.10.1)**

   ```
   bash
   curl -s https://get.nextflow.io | bash
   ```

2. **Set up a software environment** using one of:

- Docker
- Singularity (recommended for HPC)
- Podman
- Shifter
- Charliecloud
- Conda (only as a last resort)

Note - You can chain multiple config profiles, e.g., `profile test, docker`.

3. **Run your own analysis**:

```
nextflow run SFGLab/dcHiChIP \
             -profile <docker/singularity/podman/shifter/charliecloud/conda/institute> \
             --input samplesheet.csv \
             --outdir results
```

Note -

- Use **nf-core/configs** to check if your institute already has a config file.
- For **Singularity**, pre-download images using `nf-core download` and set a cache directory (`NXF_SINGULARITY_CACHEDIR`).
- For **Conda**, set a cache directory (`NXF_CONDA_CACHEDIR`) to avoid re-installing environments.
