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

---

<code style="color:red; font-weight:bold;">--ref_short</code>: Reference Genome Short Name  
- A short identifier for the reference genome build used in the workflow. This is mainly used for labeling outputs and maintaining consistency across pipeline steps. Typical values are `hg38` (human) or `mm10` (mouse).  
- Choose a short, standard genome code matching your input reference (e.g., `hg38`, `mm10`). It should correspond to the genome files (FASTA, GTF, chrom sizes) you provide.  
- *Type*: string  
- *Default flags:* hg38  

---

<code style="color:red; font-weight:bold;">--outdir</code>: Output Directory  
- Path to the main output directory where all pipeline results, logs, and intermediate files will be stored.  
- Provide either a relative or absolute path. The directory will be created automatically if it does not exist. Example: `--outdir /data/dcHiChIP_results`.  
- *Type*: string  
- *Default flags:* results  

---

<code style="color:red; font-weight:bold;">--threads</code>: Number of CPU Threads  
- Specifies the default number of CPU threads to allocate for each process, unless a specific module overrides it.  
- Use this to control parallel execution and optimize runtime based on available cores. Example: `--threads 16` for a 16-core machine. Increasing threads improves speed but also increases memory usage.  
- *Type*: integer  
- *Default flags:* 8  

---

<code style="color:red; font-weight:bold;">--mem</code>: Memory Allocation (GB)  
- Specifies the default amount of memory (in gigabytes) to allocate per process, unless overridden by module-specific settings.  
- Set this according to the available system memory and dataset size. For example, `--mem 16` allocates 16 GB per task. Insufficient memory may cause process failures on large datasets.  
- *Type*: integer  
- *Default flags:* 4  

## Alignment & Filtering

Configuration options for read alignment, mapping quality filtering, duplicate removal, and sorting.  
These parameters control how raw HiChIP reads are aligned to the reference genome and preprocessed before downstream analyses.  

*Help:* Adjust these settings to balance mapping stringency, computational cost, and data quality.  
For example, lowering the mapping quality (`mapq`) may include more reads but increase background noise, while increasing it improves precision.

---

<code style="color:red; font-weight:bold;">--mapq</code>: Minimum Mapping Quality (MAPQ)  
- Sets the minimum MAPQ threshold for retaining aligned reads. Reads with mapping quality below this value are filtered out before downstream analysis.  
- Use higher MAPQ values (e.g., 30) to include only confidently aligned reads and reduce background noise. Lower values may retain more reads but can increase false positives. Typical range: 10–60.  
- *Type*: integer  
- *Default flags:* 30  

---

<code style="color:red; font-weight:bold;">--se_samtools_args</code>: SAMtools Arguments (Single-End)  
- Custom command-line flags passed to SAMtools when processing single-end read alignments.  
- Use this to fine-tune how SAMtools handles single-end BAM/SAM files — for example, adjusting compression, threading, or output format. Example: `--se_samtools_args "-@ 8 -bh"`. Leave empty (null) to use default SAMtools behavior.  
- *Type*: string  
- *Default flags:* null  

---

<code style="color:red; font-weight:bold;">--se_bwa_mem_args</code>: BWA-MEM Arguments (Single-End)  
- Additional command-line parameters for BWA-MEM when aligning single-end reads.  
- Use this to customize BWA-MEM behavior for single-end reads — for example, seed length (`-k`), mismatch penalties, or output verbosity. Example: `--se_bwa_mem_args "-k 19 -B 4 -O 6,6 -E 1,1"`. Leave empty (null) to use default alignment settings.  
- *Type*: string  
- *Default flags:* null  

---

<code style="color:red; font-weight:bold;">--bwa_mem_args</code>: BWA-MEM Arguments (Paired-End)  
- Custom command-line options for the BWA-MEM aligner used for paired-end read mapping.  
- Use this to modify BWA-MEM parameters such as alignment scoring, read group tagging, or reporting options. The default `-M -v 0` marks shorter split hits as secondary and suppresses verbose output. Example: `--bwa_mem_args "-M -K 100000000 -Y -R '@RG\tID:sample\tSM:sample'"`.  
- *Type*: string  
- *Default flags:* -M -v 0  

---

<code style="color:red; font-weight:bold;">--samtools_fixmate_args</code>: SAMtools Fixmate Arguments  
- Specifies additional command-line options for the `samtools fixmate` command, which ensures that read-pair information is consistent in the alignment file.  
- The default `-m` option marks missing mate reads and adds mate coordinate information to each read pair. Adjust this if you need to control how mate tags or secondary alignments are handled. Example: `--samtools_fixmate_args "-m -O bam"`.  
- *Type*: string  
- *Default flags:* -m  

---

<code style="color:red; font-weight:bold;">--remove_duplicates_args</code>: Duplicate Removal Arguments  
- Specifies custom command-line flags for the duplicate removal step, controlling how PCR or optical duplicates are identified and filtered.  
- The default `-n` flag skips duplicate removal but still counts duplicates for statistics. Modify this if you want to fully remove duplicates or change behavior depending on your data. Example: `--remove_duplicates_args "-n --stats dupstats.txt"`.  
- *Type*: string  
- *Default flags:* -n  

---

<code style="color:red; font-weight:bold;">--filter_quality_args</code>: Read Quality Filtering Arguments  
- Additional command-line options for filtering reads based on mapping or sequence quality before downstream analysis.  
- Use this to apply extra filtering thresholds beyond MAPQ, such as minimum alignment length or mismatch rate. Leave empty (null) to use default filtering behavior. Example: `--filter_quality_args "--min-MAPQ 10 --min-len 30"`.  
- *Type*: string  
- *Default flags:* null  

---

<code style="color:red; font-weight:bold;">--filter_paires_args</code>: Pair Filtering Arguments  
- Custom options for filtering valid read pairs based on distance, orientation, or pairing criteria during HiChIP read preprocessing.  
- Use this to refine which read pairs are retained for contact map generation. For example, you can exclude pairs beyond a distance threshold or with invalid orientations. Example: `--filter_paires_args "--max-dist 2000 --min-dist 100"`. Leave empty (null) to keep default pair filtering behavior.  
- *Type*: string  
- *Default flags:* null  

---

<code style="color:red; font-weight:bold;">--samtools_markdup_args</code>: SAMtools Markdup Arguments  
- Specifies additional options for the `samtools markdup` command, which marks or removes duplicate reads in BAM files.  
- Use this to control duplicate marking behavior, threading, or reporting options. For example, adding `-r` removes duplicates instead of marking them, and `-@ 8` sets the number of threads. Example: `--samtools_markdup_args "-r -@ 8"`. Leave empty (null) to use SAMtools default settings.  
- *Type*: string  
- *Default flags:* null  

---

<code style="color:red; font-weight:bold;">--samtools_sort_2_args</code>: SAMtools Sort (Second Pass) Arguments  
- Specifies additional parameters for the second sorting step performed by `samtools sort`, typically used for name-sorting or coordinate-sorting read pairs.  
- The default `-n` flag sorts alignments by read name, which is often required for pairwise operations. You can modify this to coordinate-sort (`-o`) or add threading options. Example: `--samtools_sort_2_args "-@ 8 -T tmp -o sorted.bam"`.  
- *Type*: string  
- *Default flags:* -n  










































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
