---
title: Parameters for the Modules
contributors:
description: Detailed overview of configurable parameters for each module in the dcHiChIP pipeline.
toc: false
---

---

The dcHiChIP pipeline provides flexible module-level control through a comprehensive set of parameters. Below is an overview of key parameters categorized by module, with default values and usage notes.

---

## General Options

General pipeline execution options (output, resources, etc.)

---

<code style="color:red; font-weight:bold;">--input</code>: Input Samplesheet (CSV)

- Path to a comma-separated file containing information about the samples in the experiment.
- You must create a design file with sample information before running the pipeline.
- This file should contain three columns with a header row, typically representing sample ID, read 1 path, and read 2 path.
- The file must have a .csv extension and follow the required structure for successful parsing.
- If you need assistance preparing your sample sheet, you can refer to the `example templates` available here:
  - [dcHiChIP Working Cases](https://sfglab.github.io/dcHiChIP/working_cases) - Select the case that matches your experimental setup.
- _Example_: `--input /data/samplesheet/samplesheet.csv`
- _Type_: string
- _Format_: file-path (CSV)
- _Default flags_: None

---

<code style="color:red; font-weight:bold;">--fasta</code>: Reference Genome FASTA (with BWA Index)

- Specifies the path to the reference genome FASTA file used for read alignment and genome indexing.
- The provided FASTA file must be accompanied by its corresponding BWA index files, typically including:

```
Homo_sapiens_assembly38.fasta
Homo_sapiens_assembly38.fasta.amb
Homo_sapiens_assembly38.fasta.ann
Homo_sapiens_assembly38.fasta.bwt
Homo_sapiens_assembly38.fasta.pac
Homo_sapiens_assembly38.fasta.sa
Homo_sapiens_assembly38.fasta.fai
```

- Ensure all index files are present in the same directory as the FASTA file before running the pipeline.
- Missing or mismatched index files can cause alignment failures.
- _Example_: `--fasta "/data/bwa_index/hg38/Homo_sapiens_assembly38.fasta"`
- This parameter is required for the alignment module and should match the genome build specified by --ref_short (e.g., hg38, mm10).
- Type: string
- Format: file-path
- Default flags: None

---

<code style="color:red; font-weight:bold;"> --ref_short</code>: Reference Genome Short Name

- A short identifier for the reference genome build used in the workflow. This is mainly used for labeling outputs and maintaining consistency across pipeline steps. Typical values are `hg38` (human) or `mm10` (mouse).
- Choose a short, standard genome code matching your input reference (e.g., `hg38`, `mm10`). It should correspond to the genome files (FASTA, GTF, chrom sizes) you provide.
- _Type_: string
- _Default flags_: hg38

---

<code style="color:red; font-weight:bold;"> --outdir</code>: Output Directory

- Path to the main output directory where all pipeline results, logs, and intermediate files will be stored.
- Provide either a relative or absolute path. The directory will be created automatically if it does not exist.
- _Example_: `--outdir "/data/dcHiChIP_results"`.
- _Type_: string
- _Default flags_: results

---

## Alignment & Filtering

Configuration options for read alignment, mapping quality filtering, duplicate removal, and sorting.

These parameters control how raw HiChIP reads are aligned to the reference genome and preprocessed before downstream analyses.

_Help_: Adjust these settings to balance mapping stringency, computational cost, and data quality.

For example, lowering the mapping quality (`mapq`) may include more reads but increase background noise, while increasing it improves precision.

---

<code style="color:red; font-weight:bold;"> --mapq</code>: Minimum Mapping Quality (MAPQ)

- Sets the minimum MAPQ threshold for retaining aligned reads. Reads with mapping quality below this value are filtered out before downstream analysis.
- Use higher MAPQ values (e.g., 30) to include only confidently aligned reads and reduce background noise. Lower values may retain more reads but can increase false positives. Typical range: 10–60.
- _Type_: integer
- _Default flags_: 30

---

<code style="color:red; font-weight:bold;"> --se_samtools_args</code>: SAMtools Arguments (Single-End)

- Custom command-line flags passed to SAMtools when processing single-end read alignments.
- Use this to fine-tune how SAMtools handles single-end BAM/SAM files — for example, adjusting compression, threading, or output format.
- _Example_: `--se_samtools_args "-@ 8 -bh"`. Leave empty (null) to use default SAMtools behavior.
- _Type_: string
- _Default flags_: null

---

<code style="color:red; font-weight:bold;"> --se_bwa_mem_args</code>: BWA-MEM Arguments (Single-End)

- Additional command-line parameters for BWA-MEM when aligning single-end reads.
- Use this to customize BWA-MEM behavior for single-end reads — for example, seed length (`-k`), mismatch penalties, or output verbosity.
- _Example_: `--se_bwa_mem_args "-k 19 -B 4 -O 6,6 -E 1,1"`. Leave empty (null) to use default alignment settings.
- _Type_: string
- _Default flags_: null

---

<code style="color:red; font-weight:bold;"> --bwa_mem_args</code>: BWA-MEM Arguments (Paired-End)

- Custom command-line options for the BWA-MEM aligner used for paired-end read mapping.
- Use this to modify BWA-MEM parameters such as alignment scoring, read group tagging, or reporting options. The default `-M -v 0` marks shorter split hits as secondary and suppresses verbose output.
- _Example_: `--bwa_mem_args "-M -K 100000000 -Y -R '@RG\tID:sample\tSM:sample'"`.
- _Type_: string
- _Default flags_: -M -v 0

---

<code style="color:red; font-weight:bold;"> --samtools_fixmate_args</code>: SAMtools Fixmate Arguments

- Specifies additional command-line options for the `samtools fixmate` command, which ensures that read-pair information is consistent in the alignment file.
- The default `-m` option marks missing mate reads and adds mate coordinate information to each read pair. Adjust this if you need to control how mate tags or secondary alignments are handled.
- _Example_: `--samtools_fixmate_args "-m -O bam"`.
- _Type_: string
- _Default flags_: -m

---

<code style="color:red; font-weight:bold;">--optical_duplicate_distance</code>: Optical Duplicate Distance

- Defines the maximum pixel distance between clusters on the flowcell that are considered optical duplicates during duplicate marking.
- Set this value to detect and optionally remove optical duplicates generated by sequencing instruments.
- A value of `0` disables optical duplicate detection. Typical values range between `100` and `2500`, depending on the sequencing platform.
- _Example_: `--optical_duplicate_distance 2500`.
- _Type_: integer
- _Default flags_: 0

---

<code style="color:red; font-weight:bold;"> --remove_duplicates_args</code>: Duplicate Removal Arguments

- Specifies custom command-line flags for the duplicate removal step, controlling how PCR or optical duplicates are identified and filtered.
- The default `-n` flag skips duplicate removal but still counts duplicates for statistics. Modify this if you want to fully remove duplicates or change behavior depending on your data.
- _Example_: `--remove_duplicates_args "-n --stats dupstats.txt"`.
- _Type_: string
- _Default flags_: -n

---

<code style="color:red; font-weight:bold;"> --filter_quality_args</code>: Read Quality Filtering Arguments

- Additional command-line options for filtering reads based on mapping or sequence quality before downstream analysis.
- Use this to apply extra filtering thresholds beyond MAPQ, such as minimum alignment length or mismatch rate. Leave empty (null) to use the default filtering behavior.
- _Example_: `--filter_quality_args "--min-MAPQ 10 --min-len 30"`.
- _Type_: string
- _Default flags_: null

---

<code style="color:red; font-weight:bold;"> --filter_paires_args</code>: Pair Filtering Arguments

- Custom options for filtering valid read pairs based on distance, orientation, or pairing criteria during HiChIP read preprocessing.
- Use this to refine which read pairs are retained for contact map generation. For example, you can exclude pairs beyond a distance threshold or with invalid orientations.
- _Example_: `--filter_paires_args "--max-dist 2000 --min-dist 100"`. Leave empty (null) to keep default pair filtering behavior.
- _Type_: string
- _Default flags_: null

---

<code style="color:red; font-weight:bold;"> --samtools_markdup_args</code>: SAMtools Markdup Arguments

- Specifies additional options for the `samtools markdup` command, which marks or removes duplicate reads in BAM files.
- Use this to control duplicate marking behavior, threading, or reporting options. For example, adding `-r` removes duplicates instead of marking them.
- _Example_: `--samtools_markdup_args "-r"`. Leave empty (null) to use SAMtools default settings.
- _Type_: string
- _Default flags_: null

---

<code style="color:red; font-weight:bold;"> --samtools_sort_2_args</code>: SAMtools Sort (Second Pass) Arguments

- Specifies additional parameters for the second sorting step performed by `samtools sort`, typically used for name-sorting or coordinate-sorting read pairs.
- The default `-n` flag sorts alignments by read name, which is often required for pairwise operations. You can modify this to coordinate-sort (`-o`) or add threading options.
- _Example_: `--samtools_sort_2_args "-@ 8 -T tmp -o sorted.bam"`.
- _Type_: string
- _Default flags_: -n

---

<code style="color:red; font-weight:bold;"> --bwa_mem_samtools_args</code>: BWA-MEM + SAMtools Combined Arguments

- Specifies additional command-line options applied jointly to the BWA-MEM alignment output and subsequent SAMtools processing steps.
- Use this to adjust output conversion and compression parameters after BWA-MEM alignment.
- For example, `-bh` converts to BAM format with the header included. Leave `null` to use the default behavior.
- _Example_: `--bwa_mem_samtools_args "-bh -@ 8"`
- _Type_: string
- _Default flags_: null

---

## Reference & Annotation

Configuration for reference genome files and annotation resources used throughout the pipeline.

These parameters ensure consistent genome build usage across mapping, feature annotation, and downstream analyses.

_Help_: Provide correct and consistent genome resources (FASTA index, GTF, chromosome sizes, and blacklist) matching your chosen reference build (e.g., hg38, mm10). Mismatched files can lead to coordinate errors or missing annotations.

---

<code style="color:red; font-weight:bold;"> --jaspar_motif</code>: JASPAR Motif File (TSV)

- Specifies the URL or local path to the JASPAR motif file (TSV format) used for motif scanning and peak annotation — typically representing transcription factor binding motifs such as CTCF.
- Provide the path or direct download link to a JASPAR motif file compatible with your reference genome. The default points to the CTCF motif (MA0139.1) for hg38.
- _Example_: `--jaspar_motif http://expdata.cmmt.ubc.ca/JASPAR/downloads/UCSC_tracks/2022/hg38/MA0139.1.tsv.gz`.
- _Type_: string
- _Default flags_: http://expdata.cmmt.ubc.ca/JASPAR/downloads/UCSC_tracks/2022/hg38/MA0139.1.tsv.gz

---

<code style="color:red; font-weight:bold;"> --blacklist</code>: ENCODE Blacklist Regions (BED)

- Path or URL to a BED file containing ENCODE blacklist regions that should be excluded from peak calling, loop detection, and coverage calculations. These regions are known to produce artificially high signal or mapping artifacts.
- Use the appropriate blacklist file for your genome build (e.g., hg19, hg38, mm10). The default points to the ENCODE hg38 blacklist (`ENCFF356LFX`).
- _Example_: `--blacklist /data/blacklist_region/hg38/hg38-blacklist.bed.gz`.
- _Type_: string
- _Default flags_: https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz

---

<code style="color:red; font-weight:bold;"> --gtf</code>: Gene Annotation File (GTF)

- Specifies the path or URL to the GTF file containing gene annotations for the reference genome. This file is used to assign peaks, loops, and other genomic features to known genes and transcripts.
- Provide a GTF file compatible with your chosen reference genome (e.g., Ensembl or GENCODE format). The default points to the GRCh38 annotation from the Illumina iGenomes collection.
- _Example_: `--gtf /data/gtf_file/hg38/genes.gtf`.
- _Type_: string
- _Default flags_: s3://ngi-igenomes/igenomes/Homo_sapiens/NCBI/GRCh38/Annotation/Genes/genes.gtf

---

<code style="color:red; font-weight:bold;"> --chrom_size</code>: Chromosome Sizes File

- Path to a two-column file listing chromosome names and their corresponding lengths, used to define genomic bounds during binning and matrix construction.
- This file ensures that cooler and related modules apply consistent chromosome boundaries. It must match the same reference genome build as the input FASTA and GTF files.
- _Example_: `--chrom_size /data/chrom_size/hg38/hg38.chrom.sizes`.
- _Type_: string
- _Default flags_: $projectDir/assets/hg38.chrom.sizes

---

<code style="color:red; font-weight:bold;"> --genomics_features</code>: Genomic Features File (MAPS)

- Specifies the path or URL to the genomic features file required by MAPS for loop calling. This file contains mappability, GC content, and restriction enzyme fragment information used for bias correction.
- Use the appropriate MAPS genomic features file matching your genome build and restriction enzyme (e.g., MboI or DpnII). The default points to the hg38 MboI 10 kb resolution file.
- _Example_: `--genomics_features /data/MAPS/features/MAPS_features_hg38_MboI_10kb.txt`.
- _Type_: string
- _Default flags_: https://raw.githubusercontent.com/HuMingLab/MAPS/refs/heads/master/MAPS_data_files/hg38/genomic_features/F_GC_M_MboI_10Kb_el.GRCh38.txt

## Peaks & Loops

Configuration of parameters related to peak calling and chromatin loop detection.

These settings control how enriched regions (peaks) and chromatin interactions (loops) are identified from HiChIP data.

_Help_: Adjust these options to fine-tune sensitivity and specificity in peak and loop detection.

For example, you can modify the p-value threshold for MACS3 or enable MAPS loop calling for specific experimental setups.

---

<code style="color:red; font-weight:bold;"> --peak_quality</code>: Peak Significance Threshold

- Sets the significance cutoff (p-value or q-value, depending on MACS3 configuration) for peak calling to identify enriched regions.
- Lower values (e.g., 0.01 or 1e-5) increase stringency, reducing false positives but possibly missing weaker peaks. The default of 0.05 is a balanced threshold.
- _Example_: `--peak_quality 0.05`.
- _Type_: number
- _Default flags_: 0.05

---

<code style="color:red; font-weight:bold;"> --genome_size</code>: Genome Size (MACS3)

- Specifies the genome size shortcut for MACS3 peak calling, corresponding to the total mappable genome length.
- Use MACS3-supported short codes like `hs` (human), `mm` (mouse), or provide an exact value (e.g., 2.7e9). Must match your reference genome build.
- _Example_: `--genome_size hs`.
- _Type_: string
- _Default flags_: hs

---

<code style="color:red; font-weight:bold;"> --macs3_callpeak_args</code>: MACS3 Additional Arguments

- Optional custom arguments passed directly to the MACS3 `callpeak` command for advanced configuration.
- Use this to customize peak calling — for example, enabling model-free mode, adjusting shift/extension sizes, or specifying control input.
- _Example_: `--macs3_callpeak_args "-q 0.01 --nomodel --shift -75 --extsize 150"`.
- _Type_: string
- _Default flags_: null

---

<code style="color:red; font-weight:bold;"> --maps_args</code>: MAPS Loop Calling Arguments

- Additional command-line flags for the MAPS (Model-based Analysis of PLAC-seq/HiChIP) tool, used for chromatin loop detection.
- Adjust MAPS behavior such as bin size, distance range, or resolution.
- _Example_: `--maps_args "--bin-size 5000 --cis-only"`. Leave empty (null) to use the default MAPS settings.
- _Type_: string
- _Default flags_: null

---

<code style="color:red; font-weight:bold;"> --skip_maps</code>: Skip MAPS Module

- Determines whether to skip the MAPS loop calling module during the pipeline execution.
- Set this to `true` to disable loop calling (useful for QC or peak-only runs). Set to `false` to perform full MAPS-based loop analysis.
- _Example_: `--skip_maps false`.
- _Type_: boolean
- _Default flags_: false

---

## Contact Matrices & Binning

Settings controlling the generation, binning, and normalization of contact matrices from HiChIP data.

These parameters determine how raw read pairs are converted into multi-resolution .cool files for downstream analysis.

_Help_: Adjust these options to tune the resolution and structure of the resulting chromatin contact maps.

For example, smaller bin sizes (e.g., 1 kb) provide higher resolution but require more memory and processing time, while larger bins (e.g., 10 kb or 25 kb) yield smoother matrices for large-scale analyses.

---

<code style="color:red; font-weight:bold;"> --cool_bin</code>: Base Bin Size (bp)

- Defines the base bin size, in base pairs, used for generating contact matrices in the .cool format.
- Smaller bin sizes (e.g., 1000 bp) provide higher resolution but increase memory and runtime requirements.
- Larger bins (e.g., 5000–25000 bp) are recommended for lower-depth datasets.
- _Example_: `--cool_bin 5000`.
- _Type_: integer
- _Default flags_: 1000

---

<code style="color:red; font-weight:bold;"> --cooler_cload_args</code>: Cooler Cload Arguments

- Specifies command-line options passed to the `cooler cload` command for loading read pairs into a .cool file.
- Use this to modify how columns from the pairs file are interpreted when generating contact matrices. The default is configured for standard pairtools output.
- _Example_: `--cooler_cload_args "pairs --zero-based -c1 1 -p1 2 -c2 3 -p2 4"`.
- _Type_: string
- _Default flags_: pairs --zero-based -c1 2 -p1 3 -c2 4 -p2 5

---

<code style="color:red; font-weight:bold;"> --cooler_zoomify_res</code>: Zoomify Resolution Preset

- Specifies the resolution preset key used for multi-resolution cooler files (`cooler zoomify`).
- Each key corresponds to a predefined set of bin sizes defined in `insulation_resultions`.
- Use this to select which resolution set to apply during matrix zoomification. For example, `1000N` corresponds to 1 kb–500 kb bins by default.
- _Example_: `--cooler_zoomify_res 1000N`.
- _Type_: string
- _Default flags_: 1000N

---

<code style="color:red; font-weight:bold;"> --cooler_zoomify_args</code>: Cooler Zoomify Additional Arguments

- Extra command-line options for the `cooler zoomify` command that generate multi-resolution cooler files.
- Use this to control balancing, overwrite behavior, or other Zoomify parameters.
- _Example_: `--cooler_zoomify_args "--balance --force"`. Leave empty ("") to use defaults.
- _Type_: string
- _Default flags_: null

---

<code style="color:red; font-weight:bold;"> --insulation_resultions</code>: Insulation Score Resolutions

- Specifies the resolution presets (in base pairs) used for insulation score and TAD boundary calculations.
- Each key defines a named preset containing multiple window sizes.
- Provide one or more resolution sets, where each key (e.g., 1000N) maps to a list of window sizes that correspond to the resolutions used for cooler zoomify. These determine the granularity of TAD detection and insulation profiling..
- _Example_: `--insulation_resultions {"1000N": "1000 2000 5000 10000 20000 50000 100000 200000 500000"}'
- _Type_: object
- _Default flags_: {"1000N": "1000 2000 5000 10000 20000 50000 100000 200000 500000"}

---

<code style="color:red; font-weight:bold;"> --cooltools_eigscis_args</code>: Cooltools Eigs-Cis Arguments

- Additional command-line flags for the `cooltools eigs-cis` command, which computes eigenvectors for A/B compartment analysis.
- Modify this to control the number of eigenvectors or specify contact type (cis/trans).
- _Example_: `--cooltools_eigscis_args "--n-eigs 2 --contact-type cis"`.
- _Type_: string
- _Default flags_: --n-eigs 1

---

<code style="color:red; font-weight:bold;"> --cooler_eigscis_resultion</code>: Eigenvector Resolution (bp)

- Specifies the resolution, in base pairs, used for eigenvector computation in A/B compartment analysis.
- Smaller values provide higher resolution but require denser contact maps.
- _Example_: `--cooler_eigscis_resultion 100000`.
- _Type_: integer
- _Default flags_: 100000

---

<code style="color:red; font-weight:bold;"> --calder_bin</code>: CALDER Bin Size

- Sets the bin size used for subcompartment calling by CALDER.
- Provide bin size as a number (e.g., 10000). Higher bin sizes reduce resolution but speed up computation.
- _Example_: `--calder_bin 25000`.
- _Type_: string
- _Default flags_: 100000

---

<code style="color:red; font-weight:bold;"> --calder_chrom</code>: CALDER Chromosomes

- Defines the chromosomes to be analyzed by the CALDER subcompartment detection module.
- You can specify a single chromosome or multiple chromosomes separated by commas.
- When providing multiple entries, do not include spaces between chromosome numbers.
- _Example_:
  - For a single chromosome: `--calder_chrom "1"`
  - For multiple chromosomes: `--calder_chrom "1,2,3,19,20"`
- _Type_: string
- _Default flags_: 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22

---

<code style="color:red; font-weight:bold;"> --gstripe_args</code>: gStripe Detection Arguments

- Extra command-line options for the `gStripe` tool, which detects stripe-like chromatin interaction patterns.
- Modify detection thresholds or bin correction behavior. The default `--fix_bin_start` ensures stripe detection starts from bin-aligned positions.
- _Example_: `--gstripe_args "--fix_bin_start --minlen 3 --qval 0.05"`.
- _Type_: string
- _Default flags_: --fix_bin_start

---

## Visualization & QC

Configuration for plotting and quality-control steps (e.g., FastQC, deepTools, Juicer Tools, pairtools).

These options affect coverage plots, sample correlations, and other diagnostic visualizations.

_Help_: Tune these to control plot styles, correlation methods, and output formats.

For large cohorts, consider lighter settings (e.g., skipping numbers on heatmaps) to speed up rendering.

---

<code style="color:red; font-weight:bold;"> --plot_method</code>: Correlation Method for Plots

- Defines the statistical method used to compute sample-to-sample correlations in deepTools visualizations.
- Choose between `spearman` (rank-based, robust to outliers) or `pearson` (linear correlation).
- The default `spearman` is typically used for read-count correlation heatmaps.
- _Example_: `--plot_method pearson`.
- _Type_: string
- _Default flags_: spearman

---

<code style="color:red; font-weight:bold;"> --plot_type</code>: Plot Type for Correlation Visualization

- Specifies the visualization format for correlation results — typically a heatmap or scatter plot.
- Select the style best suited to your comparison. `heatmap` provides an overview of pairwise correlations, while `scatter` highlights relationships between specific samples.
- _Example_: `--plot_type scatter`.
- _Type_: string
- _Default flags_: heatmap

---

<code style="color:red; font-weight:bold;"> --fastqc_args</code>: FastQC Arguments

- Specifies additional command-line options for FastQC quality control analysis of raw sequencing reads.
- Use this to control FastQC behavior, such as verbosity or output compression. The default `--quiet` suppresses detailed logging.
- _Example_: `--fastqc_args "--quiet --nogroup"`.
- _Type_: string
- _Default flags_: --quiet

---

<code style="color:red; font-weight:bold;"> --deeptools_plotcoverage_args</code>: deepTools PlotCoverage Arguments

- Additional flags for the deepTools `plotCoverage` module, which visualizes genome-wide read coverage across samples.
- Customize output options such as file format, scaling, or whether to skip empty bins. The default `--skipZeros` ignores regions with zero coverage. Example: `--deeptools_plotcoverage_args "--skipZeros --plotFileFormat pdf"`.
- _Type_: string
- _Default flags_: --skipZeros

---

<code style="color:red; font-weight:bold;"> --deeptools_plotcorrelation_args</code>: deepTools PlotCorrelation Arguments

- Custom flags for the deepTools `plotCorrelation` command, which computes and visualizes sample correlation matrices.
- Use this to modify correlation type, color scheme, and plot annotations. The default parameters produce a Spearman correlation heatmap with labeled values.
- _Example_: `--deeptools_plotcorrelation_args "--corMethod pearson --colorMap RdBu --what pairwise"`.
- _Type_: string
- _Default flags_: --skipZeros --plotTitle "Spearman Correlation of Read Counts" --colorMap RdYlBu --plotNumbers

---

<code style="color:red; font-weight:bold;"> --juicertools_args</code>: Juicer Tools Arguments

- Specifies extra command-line parameters for Juicer Tools utilities used for Hi-C/HiChIP contact map analysis.
- Set parameters for operations like normalization, balancing, or matrix extraction. Leave null to use the default Juicer settings.
- _Example_: `--juicertools_args "pre --resolutions 5000"`
- _Type_: string
- _Default flags_: null

---

<code style="color:red; font-weight:bold;"> --pairtools_parse2_args</code>: Pairtools Parse2 Arguments

- Additional options for the `pairtools parse2` command used to convert alignment files into a standardized pairs format.
- Modify column selection, filtering, or threading for parsing steps. Leave null for default settings.
- _Example_: `--pairtools_parse2_args "--add-columns chr1,chr2,pos1,pos2"`.
- _Type_: string
- _Default flags_: null

---

## 3D Genome Modelling (MultiMM)

Configuration options for 3D genome reconstruction and visualization using the MultiMM module.

These parameters define which genomic regions are modeled, the computational platform used (CPU or GPU), and optional fine-grained settings for chromosomal or locus-level simulations.

_Help_: Adjust these parameters to specify the modeling scope (whole chromosome, gene region, or specific locus) and computational mode.

MultiMM integrates chromatin contact data to predict spatial genome structures. For large models, prefer GPU acceleration if available.

---

<code style="color:red; font-weight:bold;">--multimm_platform</code>: MultiMM Computational Platform

- Specifies the compute platform used for MultiMM 3D genome modeling — either CPU or GPU-based execution, depending on system availability and dataset size.
- The following flags can be used:
  - `CPU` — standard desktop or cluster execution (default).
  - `OpenCL` — enables GPU acceleration using the OpenCL framework for compatible hardware.
  - `CUDA` — enables GPU acceleration using NVIDIA CUDA, providing faster simulations for large or high-resolution models.
- _Example_: `--multimm_platform CUDA`
- _Type_: string
- _Default flags_: CPU

---

<code style="color:red; font-weight:bold;"> --multimm_modelling_level</code>: Modelling Granularity Level

- Specifies the genomic scale at which MultiMM performs 3D genome modeling — ranging from a single gene to the entire genome.
- Choose the modeling depth according to your analysis goal and computational resources.
- The following modelling levels are available:
  - **GENE** — Provide a gene of interest (with an associated .bedpe file path). MultiMM models the gene with a default ±100 kb window around it. Compartment forces are _not_ considered at this level.
  - **REGION** — Specify a chromosome and genomic coordinates (start–end). Compartment interactions can optionally be included. Only the selected region is modeled.
  - **CHROMOSOME** — Provide a chromosome name; MultiMM automatically determines start and end coordinates. Compartments of data have to be imported.
  - **GW (Genome-Wide)** — Models the entire genome. No input for chromosome or coordinates is needed. Compartments of data have to be imported. This is the most computationally intensive mode, potentially taking minutes to hours depending on system performance.
- _Example_: `--multimm_modelling_level region`
- _Type_: string
- _Default flags_: chrom

---

<code style="color:red; font-weight:bold;"> --multimm_gene_name</code>: Gene Name for Modelling

- Specifies the gene symbol or identifier to be modeled in 3D when the modelling level is set to `gene`.
- Provide a valid gene symbol that exists within your reference annotation (e.g., `PIRAT1`, `RUNX1`). This parameter is ignored if the modelling level is set to `chrom` or `region`.
- _Example_: `--multimm_gene_name RUNX1`.
- _Type_: string
- _Default flags_: null

---

<code style="color:red; font-weight:bold;"> --multimm_chrom</code>: Chromosome to Model

- Defines which chromosome should be used for 3D genome modeling. Required for both _chromosome_ level and _region_ level modeling.
- Use standard chromosome naming consistent with your reference genome (e.g., `chr1`, `chr21`, `chrX`).
- _Example_: `--multimm_chrom chr10`.
- _Type_: string
- _Default flags_: chr21

---

<code style="color:red; font-weight:bold;"> --multimm_loc_start</code>: Locus Start Coordinate (bp)

- Specifies the genomic start coordinate (in base pairs) for _region_ level 3D modeling.
- Used only when the modelling level is set to `locus`. Must be a valid coordinate within the selected chromosome.
- _Example_: `--multimm_loc_start 28000000`.
- _Type_: integer
- _Default flags_: null

---

<code style="color:red; font-weight:bold;"> --multimm_loc_end</code>: Locus End Coordinate (bp)

- Specifies the genomic end coordinate (in base pairs) for _region_ level 3D modeling.
- Used together with `multimm_loc_start` to define the genomic window for 3D reconstruction.
- _Example_: `--multimm_loc_end 32000000`.
- _Type_: integer
- _Default flags_: null

---

<code style="color:red; font-weight:bold;"> --multimm_args</code>: MultiMM Additional Arguments

- Custom command-line arguments passed directly to MultiMM for fine-tuning modeling behavior and output.
- Use this for advanced control — e.g., number of models, iteration limits, or output directory.
- _Example_: `--multimm_args "--n-models 50 --max-iters 2000"`. Leave empty (`null`) to use defaults.
- _Type_: string
- _Default flags_: null

## **Notes on MultiMM Modelling Levels**

- _GENE level_

  - The user must provide the gene name using `--multimm_gene_name`.
  - Parameters such as `--multimm_chrom`, `--multimm_loc_start`, and `--multimm_loc_end` should remain **null**, as they are not required at this level.
  - Suitable for detailed, local reconstruction of a single gene and its regulatory environment.

- _REGION level_

  - The user must provide the chromosome name (`--multimm_chrom`), chromosome start (`--multimm_loc_start`), and chromosome end (`--multimm_loc_end`) coordinates of the region.
  - The gene name (`--multimm_gene_name`) should remain **null**.
  - Used when focusing on a defined genomic segment, such as a TAD or regulatory region.

- _CHROMOSOME level_

  - Only the chromosome name (`--multimm_chrom`) must be provided.
  - Parameters `--multimm_gene_name`, `--multimm_loc_start`, and `--multimm_loc_end` should be **null**.
  - Ideal for analyzing entire chromosome architecture and its folding pattern.

- _GW (GENOME-WIDE) level_

  - No chromosome, gene, or coordinate inputs are required.
  - Parameters `--multimm_gene_name`, `--multimm_chrom`, `--multimm_loc_start`, and `--multimm_loc_end` should all be **null**.
  - Used for complete genome-wide modeling to capture inter-chromosomal and large-scale spatial organization.

- _Automatic bead assignment_
  - The modelling level automatically determines the number of simulation beads.
  - Regardless of the user-defined `N_BEADS` value, specifying `--multimm_modelling_level` overrides it with the following defaults:
    - _GENE_: 1,000 beads
    - _REGION_: 5,000 beads
    - _CHROMOSOME_: 20,000 beads
    - _GW (Genome-Wide)_: 200,000 beads

---

## Chromatin Contact Domain (CCD) Calling

Configuration options for calling Chromatin Contact Domains (CCDs) from significant chromatin loops (e.g. MAPS .bedpe output).

These parameters control how strictly loop-dense genomic regions are identified as CCDs, based on loop coverage per bin and minimum domain length.

Help:
Adjust these parameters to balance sensitivity vs specificity.
By default, the pipeline reports compact, well-supported CCDs defined by ≥ 2 overlapping loops across ≥ 15 kb (~3 contiguous 5 kb bins), which performs robustly across most HiChIP/MAPS datasets.

CCD calling in this pipeline is coverage-driven (not diff-driven): domains are reported only where loops are present, ensuring stable and biologically interpretable results.

---

<code style="color:red; font-weight:bold;">--summits_only</code>: Use Only MAPS Cluster Summits

- Restricts CCD calling to loop anchors where `ClusterSummit == 1` in MAPS output.
- This focuses domain detection on the most confident interactions (cluster peaks) and reduces noise from weaker loop edges.
- Recommended for most analyses unless you explicitly want to include all loop bin-pairs.
- If not set, all intra-chromosomal MAPS loops are used.
- _Example_: `--summits_only`
- _Type_ : flag (boolean)
- _Default flags_: enabled (via pipeline default)

---

<code style="color:red; font-weight:bold;">--min_loops</code>: Minimum Supporting Loops per CCD

- Specifies the minimum number of overlapping loops required for a genomic region to be considered part of a CCD.
- At 5 kb resolution, this corresponds to the minimum loop coverage per 5 kb bin.
- Higher values increase stringency and reduce false positives; lower values are more permissive and capture weaker domain structures.
- Typical values:
  - 1 → very permissive (exploratory)
  - 2 → balanced (default)
  - 3+ → high-confidence domains only
- _Example_: `--min_loops 3`
- _Type_: integer
- _Default flags_: 2

---

<code style="color:red; font-weight:bold;">--min_length</code>: Minimum CCD Length (bp)

- Specifies the minimum genomic span (in base pairs) required for a CCD to be reported.
- Shorter domains are filtered out after CCD calling.
- With 5 kb bins, common values translate as:
  - 10000 bp ≈ 2 bins
  - 15000 bp ≈ 3 bins (default)
  - 25000 bp ≈ 5 bins (stricter)
- Increase this value to focus on larger structural domains; decrease it to capture compact/local CCDs.
- _Example_: `--min_length 15000`
- _Type_: integer
- _Default flags_: 15000

---

<code style="color:red; font-weight:bold;">--ccd_caller_args</code>: Custom CCD Caller Arguments

- Allows full manual control over the CCD caller command-line options.
- When provided, this string overrides the pipeline’s automatically constructed arguments.
- Leave as null to use pipeline defaults:
`--summits_only --min_loops <min_loops> --min_length <min_length>`
- Useful for advanced users who want to explicitly control all CCD parameters.
- _Example (no summit filtering)_: `--min_loops 3 --min_length 20000`
- _Example (summits only)_: `--summits_only --min_loops 3 --min_length 20000`
- _Type_: string
- _Default flags_: null

---

## **Notes on CCD Logic**

- CCDs are defined **strictly by loop coverage** - regions with no loops will never form domains.
- Domain boundaries reflect **contiguous loop-supported bins**, not local noise or signal fluctuations.
- This design is especially robust at **high resolution (5 kb bins)**, where diff-based methods tend to produce artificial micro-domains.

---

## **Maximum Job Request Options**

- These parameters define the upper resource limits available to the Nextflow scheduler during pipeline execution.
- They ensure that individual processes do not exceed the maximum computational resources defined by the user or system administrator.

---

<code style="color:red; font-weight:bold;">--max_cpus</code>: Maximum Number of CPU Cores

- Specifies the maximum number of CPU cores that can be allocated to any process in the pipeline.
- Increasing this value allows more parallel execution but may require a high-performance computing (HPC) environment.
- _Example_: `--max_cpus 50`
- _Type_: integer
- Default flags: null

---

<code style="color:red; font-weight:bold;">--max_memory</code>: Maximum Memory Allocation

- Defines the maximum amount of memory that can be used by any single process in the pipeline.
- The value must include units (MB, GB, or TB). It is recommended to allocate sufficient memory for alignment and matrix generation steps.
- _Example_: `--max_memory '216.GB'`
- Type: string
- Default flags: null

---

<code style="color:red; font-weight:bold;">--max_time</code>: Maximum Execution Time

- Specifies the maximum wall-clock time allowed for any process.
- The value must include time units (e.g., h for hours, m for minutes). It helps prevent indefinitely running jobs.
- _Example_: `--max_time '1200.h'`
- _Type_: string
- _Default flags_: null

---


## Quick Copy‑Paste: params.config template

```groovy
params {

  /* ========= GENERAL OPTIONS ========= */
  ref_short = "hg38"
  outdir    = "/mnt/raid/test_case/results"
  mem       = 4

  /* ========= INPUT FILES ========= */
  input  = "/mnt/raid/test_case/samplesheet.csv"

  /* ========= REFERENCES & ANNOTATION ========= */
  jaspar_motif       = "http://expdata.cmmt.ubc.ca/JASPAR/downloads/UCSC_tracks/2022/hg38/MA0139.1.tsv.gz"
  blacklist          = "https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz"
  gtf                = "s3://ngi-igenomes/igenomes/Homo_sapiens/NCBI/GRCh38/Annotation/Genes/genes.gtf"
  chrom_size         = "$projectDir/assets/hg38.chrom.sizes"
  genomics_features  = "https://raw.githubusercontent.com/HuMingLab/MAPS/refs/heads/master/MAPS_data_files/hg38/genomic_features/F_GC_M_MboI_10Kb_el.GRCh38.txt"

  /* ========= MAPPING & FILTERING ========= */
  mapq                   = 30
  se_samtools_args        = null
  se_bwa_mem_args         = null
  bwa_mem_args            = "-M -v 0"
  samtools_fixmate_args   = "-m"
  optical_duplicate_distance = 0
  remove_duplicates_args  = "-n"
  filter_quality_args     = null
  filter_paires_args      = null
  samtools_markdup_args   = null
  samtools_sort_2_args    = "-n"
  bwa_mem_samtools_args   = null

  /* ========= PEAKS & LOOPS ========= */
  peak_quality         = 0.05
  genome_size          = "hs"
  macs3_callpeak_args  = null
  maps_args            = null
  skip_maps            = false

  /* ========= CONTACT MATRICES & BINNING ========= */
  cool_bin                   = 1000
  cooler_cload_args          = "pairs --zero-based  -c1 2 -p1 3 -c2 4 -p2 5"
  cooler_zoomify_res         = 1000N
  cooler_zoomify_args        = null
  insulation_resultions      = '{"1000N": "1000 2000 5000 10000 20000 50000 100000 200000 500000"}'
  cooltools_eigscis_args     = "--n-eigs 1"
  cooler_eigscis_resultion   = 100000
  calder_bin                 = 100000
  gstripe_args               = "--fix_bin_start"

/* ========= Max job request options ========= */
  max_cpus: 32
  max_memory: `216.GB`
  max_time: `120.h`


  /* ========= VISUALIZATION & QC ========= */
  plot_method                 = "spearman"
  plot_type                   = "heatmap"
  fastqc_args                 = "--quiet"
  deeptools_plotcoverage_args = "--skipZeros"
  deeptools_plotcorrelation_args = '--skipZeros --plotTitle "Spearman Correlation of Read Counts" --colorMap RdYlBu --plotNumbers'
  juicertools_args            = null
  pairtools_parse2_args       = null

  /* ========= 3D GENOME MODELLING (MULTIMM) ========= */
  multimm_platform         = "CPU"
  multimm_modelling_level  = "chrom"
  multimm_gene_name        = null
  multimm_chrom            = "chr21"
  multimm_loc_start        = null
  multimm_loc_end          = null
  multimm_args             = null

  /* ========= PIPELINE EXECUTION LIMITS ========= */
  max_cpus   = 32
  max_memory = '216.GB'
  max_time   = '120.h'
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
