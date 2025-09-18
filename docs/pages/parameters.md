---
title: Parameters for the Modules
contributors: [Abhishek Agarwal, Ziad Al-Bkhetan]
description: Detailed overview of configurable parameters for each module in the dcHiChIP pipeline.
toc: true
type: guides
---

---

The `dcHiChIP` pipeline provides flexible module-level control through a comprehensive set of parameters. Below is an overview of key parameters categorized by module, with default values and usage notes.---

### Mandatory input files

- **`max_cpus`**: Maximum Allocated CPU
- **`max_memory`**: Maximum Allocated Memory
- **`max_time`**: Maximum Execution Time
- **`input`**: Location of the samplesheet
- **`fasta`**: Location of the bwa index file
- **`genomics_features`**: Location of the genomic fetures
- **`outdir`**: Location of the output directory
- **`genome_size`**: reference genome (eg, hs or mm)

### Other Parameters

## bwa mem
- **`mapq`**: MAPping Quality (Phred-scaled)
- **`-t`**: Number of threads
---

## samtools
- **`--threads`**: Number of threads
---

## MACS3
- **`--qvalue`**: q-value (minimum FDR) cutoff to call significant regions [Default = 0.05]
---

## MAPS
- **`bin_size`**: resolution. Usually 5000 or 10000. [Default = 5000]
- **`binning_range`**: How far 3D interactions can be called, also affects the estimate of the expected count. [Default=1000000]
- **`sex_chroms_to_process`**: either X,Y,XY or NA. [Default = NA]
    - This specifies which (if any) sex chromosomes the user wants to run MAPS on: X = X chr only, Y = Y chr only, XY = both X and Y chroms, NA = none (just autosomal).
---

## DeepTools (plotCoverage,multiBigwigSummary) 
- **`--numberOfProcessors`**: Number of processors to use [Default = 1]
- **`--binSize`**: Size (in bases) of the windows sampled from the genome. [Default: 10000]
---

## HOMER(findMotifsGenome.pl)
- **`-p`**: Number of CPUs to use [Default 1]
- **`-S`**: Number of motifs to find [Default 25)
---

## Stripe Calling (gStripe)
- **`--max_workers`**: Number of CPUs to use [Default 1]
---

## Open2C (pairtools)
- **`--nproc`**: Number of processes
---

## Open2C (cooler Zomify)
- **`--nproc`**: Number of processes to use for batch processing chunks of pixels [Default: 1]
- **`--resolutions`**: Comma-separated list of target resolutions.
    - Use suffixes `B` or `N` to specify a progression:
    - `B` for binary (geometric steps of factor 2)
    - `N` for nice (geometric steps of factor 10 interleaved with steps of 2 and 5)
      - **Examples**:
        - `1000B` → 1000, 2000, 4000, 8000, …
        - `1000N` → 1000, 2000, 5000, 10000, …
        - `5000N` → 5000, 10000, 25000, 50000, …
---

## Open2C (cooltools insulation)
- **`--nproc`**: Number of processes to split the work between. [Default: 1]
- **`<res>`**: resolution [Default: 1000 2000 5000 10000 20000 50000 100000 200000 500000]
---

## Open2C (cooltools insulation)
- **`--nproc`**: Number of processes to split the work between. [Default: 1]
- **`<res>`**: resolution [Default: 100000]
---

## CALDER
- **`--ncores`**: Number of processes to split the work between. [Default: 1]
- **`<res>`**: resolution [Default: 100000]
---

## 3Dmodelling (MultiMM)
- **`--platform`**: name of the platform. Available choices: CPU, OpenCL, CUDA. [Default: CPU]
- **`--cpu_threads`**: Number of CPU threads (in case that CPU is chosen as platform).
- **`--modelling_level`**: Helping function to specify parameters of simulation.
    - Choose 'GENE', 'REGION', 'CHROM', or 'GW' depending on the resolution of interest.
- **`--gene_name`**: The name of the gene of interest. {Only valid when --modelling_level = GENE} 
- **`--chrom`**: Chromosome that corresponds the the modelling region of interest
    - case-sensivite for chr21 (not accepted : CHR21, Chr21, 21) {Only valid when --modelling_level = REGION/CHROM }
- **`--loc_start`**: Starting region coordinate {Only valid when --modelling_level = REGION}
- **`--loc_end`**: Ending region coordinate {Only valid when --modelling_level = REGION}
---




