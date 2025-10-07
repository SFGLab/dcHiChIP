---
title: Parameters for the Modules
contributors:
description: Detailed overview of configurable parameters for each module in the dcHiChIP pipeline.
toc: false
---
---
The dcHiChIP pipeline provides flexible module-level control through a comprehensive set of parameters. Below is an overview of key parameters categorized by module, with default values and usage notes.

## General Options

General pipeline execution options (output, resources, etc.)

**`--input`**: Input file-path
- Path to a comma-separated file containing information about the samples in the experiment.
- You will need to create a design file with information about the samples in your experiment before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 3 columns and a header row.

**`--output`**: Output file-path
- Path to the main output directory where all pipeline results, logs, and intermediate files will be stored
- Provide either a relative or absolute path. The directory will be created automatically if it does not exist. Example: `--outdir /data/dcHiChIP_results`.
- *Default flags:* `result`

**`--output`**: Output file-path
- Specifies the default number of CPU threads to allocate for each process, unless a specific module overrides it.
- Use this to control parallel execution and optimize runtime based on available cores. Example: `--threads 16` for a 16-core machine. Increasing threads improves speed but also increases memory usage.
- *Default flags:* `8`

**`--mem`**: Memory Allocation (GB)"
- Specifies the default amount of memory (in gigabytes) to allocate per process, unless overridden by module-specific settings.
- Set this according to the available system memory and dataset size. For example, `--mem 16` allocates 16 GB per task. Insufficient memory may cause process failures on large datasets.
- *Default flags:* `4`

## Mapping & Filtering
Parameters controlling read alignment, pair filtering, and mapping quality thresholds for dcHiChIP.

<code style="color:red; font-weight:bold;"> --bwa_index</code>: Path to BWA Index

- Path to the pre-built BWA index of the reference genome.
- This is required for aligning paired-end reads. Ensure all index files (.bwt, .sa, .ann, .amb, .pac) are in the same directory.
- Example: --bwa_index /data/genomes/hg38/bwa_index/genome.fa
- Type: string
- Default value: None










































## Mandatory Input Parameters

- **`input`**: Location of the samplesheet  
- **`outdir`**: Location of the output directory  
- **`fasta`**: Location of the BWA index file  
- **`genomics_features`**: Location of the genomic features file  
- **`profile`**: Execution environment (choose one: `docker`, `singularity`, `podman`, `shifter`, `charliecloud`, or `conda`)  
- **`max_cpus`**: Maximum allocated CPU cores  
- **`max_memory`**: Maximum allocated memory  
- **`max_time`**: Maximum execution time  
- **`jaspar_motif`**: JASPAR motif file (`.tsv` format)  
- **`mapq`**: Minimum mapping quality threshold  
- **`peak_quality`**: P-value cutoff for peak calling  
- **`genome_size`**: Reference genome (e.g., `hs` for human, `mm` for mouse)  

## Other Parameters

- **FILTER_QUALITY**  
  *Default flags:* `-F 0x04 -q ${params.mapq}`  
  *Override param:* `params.filter_quality_args`  
  *Example:* `--filter_quality_args "-F 0x904 -q 30"`

- **BWA_MEM**  
  *Default flags (aligner):* `-M -v 0` → `params.bwa_mem_args`  
  *Default flags (samtools view):* `-bh` → `params.bwa_mem_samtools_args`  
  *Examples:* `--bwa_mem_args "-M -t 16"`, `--bwa_mem_samtools_args "-bh -@ 8"`

- **REMOVE_DUPLICATES**  
  *Default flags:* `-n`  
  *Override param:* `params.remove_duplicates_args`  
  *Example:* `--remove_duplicates_args "-r -s"`

- **MACS3_CALLPEAK**  
  *Default flags:* `--nomodel -q ${params.peak_quality} -B --format BAMPE`  
  *Override param:* `params.macs3_callpeak_args`  
  *Example:* `--macs3_callpeak_args "--nomodel -q 0.01 -B --format BAMPE --keep-dup auto"`

- **MAPS**  
  *Default flags:* `-n 1 -a ${params.mapq}`  
  *Override param:* `params.maps_args`  
  *Example:* `--maps_args "-n 2 -a 10 --res 5000"`

- **JUICERTOOLS**  
  *Default flags:* `${params.ref_short}`  
  *Override param:* `params.juicertools_args`  
  *Example:* `--juicertools_args "pre -r 10000 -k KR"`

- **DEEPTOOLS_PLOTCOVERAGE**  
  *Default flags:* `--skipZeros`  
  *Override param:* `params.deeptools_plotcoverage_args`  
  *Example:* `--deeptools_plotcoverage_args "--bins 5000 --skipZeros"`

- **GSTRIPE**  
  *Default flags:* `--fix_bin_start`  
  *Override param:* `params.gstripe_args`  
  *Example:* `--gstripe_args "--chromsizes path/to.sizes --threads 8 --out results/stripes"`

- **DEEPTOOLS_PLOTCORRELATION**  
  *Default flags:* `--cmap RdYlBu --plotNumbers`  
  *Override param:* `params.deeptools_plotcorrelation_args`  
  *Example:* `--deeptools_plotcorrelation_args "--corMethod spearman --cmap RdBu --what pairwise"`

- **PAIRTOOLS_PARSE2**  
  *Default flags:* `parse --output-type pairs`  
  *Override param:* `params.pairtools_parse2_args`  
  *Example:* `--pairtools_parse2_args "parse --min-mapq 10 --output-type pairs"`

- **COOLER_CLOAD**  
  *Default flags:* `pairs --zero-based`  
  *Override param:* `params.cooler_cload_args`  
  *Example:* `--cooler_cload_args "pairs --zero-based -r 5000"`

- **COOLER_ZOOMIFY**  
  *Default flags:* *(empty)*  
  *Override param:* `params.cooler_zoomify_args`  
  *Example:* `--cooler_zoomify_args "--balance --resolutions 5k,10k,25k"`

- **COOLTOOLS_EIGSCIS**  
  *Default flags:* `--n-eigs 1`  
  *Override param:* `params.cooltools_eigscis_args`  
  *Example:* `--cooltools_eigscis_args "--n-eigs 3 --phasing-track GC"`

- **MULTIMM**  
  *Default flags:* *(project-defined; override if needed)*  
  *Override param:* `params.multimm_args`  
  *Example:* `--multimm_args "--fdr 0.05 --min-dist 1000"`

- **SAMTOOLS_FIXMATE**  
  *Default flags:* `-m`  
  *Override param:* `params.samtools_fixmate_args`  
  *Example:* `--samtools_fixmate_args "-m -@ 8"`

- **FILTER_PAIRES**  
  *Default flags:* `--mapq ${params.mapq}`  
  *Override param:* `params.filter_paires_args`  
  *Example:* `--filter_paires_args "--mapq 20 --max-dist 2000000"`

- **SAMTOOLS_MARKDUP**  
  *Default flags:* `-r -d ${params.optical_duplicate_distance}`  
  *Override param:* `params.samtools_markdup_args`  
  *Output prefix:* `markdup`  
  *Example:* `--samtools_markdup_args "-r -@ 8 -d 2500"`

- **SAMTOOLS_SORT_2**  
  *Default flags:* `-n`  
  *Override param:* `params.samtools_sort_2_args`  
  *Example:* `--samtools_sort_2_args "-n -@ 8"`

---

## Quick Copy‑Paste: params.config template

```groovy
params {
  max_cpus   = 32
  max_memory = '216.GB'
  max_time   = '120.h'

  input  = "/mnt/raid/test_case/samplesheet.csv"
  outdir = "/mnt/raid/test_case/results"

  fasta = "/mnt/raid/workspace/reference_genome/bwa_index/bwa1/hg38/Homo_sapiens_assembly38.fasta"
  genomics_features = "/mnt/raid/dcHiChIP/MAPS/MAPS_data_files/hg38/genomic_features/F_GC_M_MboI_10Kb_el.GRCh38.txt"
  jaspar_motif = "/mnt/raid/workspace/JASPAR/MA0139.1.tsv.gz"
  gtf = "/mnt/raid/workspace/UCSC/genes.gtf"
    
  mapq = 30
  peak_quality = 0.05
  genome_size = "hs"

  fastqc_args = null
  filter_quality_args = null
  bwa_mem_args = null
  bwa_mem_samtools_args = null
  remove_duplicates_args = null
  macs3_callpeak_args = null
  maps_args = null
  juicertools_args = null
  deeptools_plotcoverage_args = null
  deeptools_plotcorrelation_args = null
  gstripe_args = null
  pairtools_parse2_args = null
  cooler_cload_args = null
  cooler_zoomify_args = null
  cooltools_eigscis_args = null
  multimm_args = null
  samtools_fixmate_args = null
  filter_paires_args = null
  samtools_markdup_args = null
  samtools_sort_2_args = null
}
```

---

## FAQs
- **Do I need to set all of these?**  
  No. Most users will never touch them. Only override if you know you need different behavior.

- **What if I want to keep defaults and add a flag?**  
  Copy the default shown above and append your extra options.

- **How can I check what actually ran?**  
  Look inside the `.command.sh` file in the Nextflow work directory for any task.
