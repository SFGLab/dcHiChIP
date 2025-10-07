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

<code style="color:red; font-weight:bold;"> --ref_short</code>: Reference Genome Short Name  
- A short identifier for the reference genome build used in the workflow. This is mainly used for labeling outputs and maintaining consistency across pipeline steps. Typical values are `hg38` (human) or `mm10` (mouse).  
- Choose a short, standard genome code matching your input reference (e.g., `hg38`, `mm10`). It should correspond to the genome files (FASTA, GTF, chrom sizes) you provide.  
- *Type*: string  
- *Default flags:* hg38  

---

<code style="color:red; font-weight:bold;"> --outdir</code>: Output Directory  
- Path to the main output directory where all pipeline results, logs, and intermediate files will be stored.  
- Provide either a relative or absolute path. The directory will be created automatically if it does not exist.
- *Example*: `--outdir /data/dcHiChIP_results`.  
- *Type*: string  
- *Default flags:* results  

---

<code style="color:red; font-weight:bold;"> --mem</code>: Memory Allocation (GB)  
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

<code style="color:red; font-weight:bold;"> --mapq</code>: Minimum Mapping Quality (MAPQ)  
- Sets the minimum MAPQ threshold for retaining aligned reads. Reads with mapping quality below this value are filtered out before downstream analysis.  
- Use higher MAPQ values (e.g., 30) to include only confidently aligned reads and reduce background noise. Lower values may retain more reads but can increase false positives. Typical range: 10–60.  
- *Type*: integer  
- *Default flags:* 30  

---

<code style="color:red; font-weight:bold;"> --se_samtools_args</code>: SAMtools Arguments (Single-End)  
- Custom command-line flags passed to SAMtools when processing single-end read alignments.  
- Use this to fine-tune how SAMtools handles single-end BAM/SAM files — for example, adjusting compression, threading, or output format.
- *Example*: `--se_samtools_args "-@ 8 -bh"`. Leave empty (null) to use default SAMtools behavior.  
- *Type*: string  
- *Default flags:* null  

---

<code style="color:red; font-weight:bold;"> --se_bwa_mem_args</code>: BWA-MEM Arguments (Single-End)  
- Additional command-line parameters for BWA-MEM when aligning single-end reads.  
- Use this to customize BWA-MEM behavior for single-end reads — for example, seed length (`-k`), mismatch penalties, or output verbosity.
- *Example*: `--se_bwa_mem_args "-k 19 -B 4 -O 6,6 -E 1,1"`. Leave empty (null) to use default alignment settings.  
- *Type*: string  
- *Default flags:* null  

---

<code style="color:red; font-weight:bold;"> --bwa_mem_args</code>: BWA-MEM Arguments (Paired-End)  
- Custom command-line options for the BWA-MEM aligner used for paired-end read mapping.  
- Use this to modify BWA-MEM parameters such as alignment scoring, read group tagging, or reporting options. The default `-M -v 0` marks shorter split hits as secondary and suppresses verbose output.
- *Example*: `--bwa_mem_args "-M -K 100000000 -Y -R '@RG\tID:sample\tSM:sample'"`.  
- *Type*: string  
- *Default flags:* -M -v 0  

---

<code style="color:red; font-weight:bold;"> --samtools_fixmate_args</code>: SAMtools Fixmate Arguments  
- Specifies additional command-line options for the `samtools fixmate` command, which ensures that read-pair information is consistent in the alignment file.  
- The default `-m` option marks missing mate reads and adds mate coordinate information to each read pair. Adjust this if you need to control how mate tags or secondary alignments are handled.
- *Example*: `--samtools_fixmate_args "-m -O bam"`.  
- *Type*: string  
- *Default flags:* -m  

---

<code style="color:red; font-weight:bold;">--optical_duplicate_distance</code>: Optical Duplicate Distance  
- Defines the maximum pixel distance between clusters on the flowcell that are considered optical duplicates during duplicate marking.  
- Set this value to detect and optionally remove optical duplicates generated by sequencing instruments.
- A value of `0` disables optical duplicate detection. Typical values range between `100` and `2500`, depending on the sequencing platform.
- *Example*: `--optical_duplicate_distance 2500`.  
- *Type*: integer  
- *Default flags:* 0
  
---

<code style="color:red; font-weight:bold;"> --remove_duplicates_args</code>: Duplicate Removal Arguments  
- Specifies custom command-line flags for the duplicate removal step, controlling how PCR or optical duplicates are identified and filtered.  
- The default `-n` flag skips duplicate removal but still counts duplicates for statistics. Modify this if you want to fully remove duplicates or change behavior depending on your data.
- *Example*: `--remove_duplicates_args "-n --stats dupstats.txt"`.  
- *Type*: string  
- *Default flags:* -n  

---

<code style="color:red; font-weight:bold;"> --filter_quality_args</code>: Read Quality Filtering Arguments  
- Additional command-line options for filtering reads based on mapping or sequence quality before downstream analysis.  
- Use this to apply extra filtering thresholds beyond MAPQ, such as minimum alignment length or mismatch rate. Leave empty (null) to use the default filtering behavior.
- *Example*: `--filter_quality_args "--min-MAPQ 10 --min-len 30"`.  
- *Type*: string  
- *Default flags:* null  

---

<code style="color:red; font-weight:bold;"> --filter_paires_args</code>: Pair Filtering Arguments  
- Custom options for filtering valid read pairs based on distance, orientation, or pairing criteria during HiChIP read preprocessing.  
- Use this to refine which read pairs are retained for contact map generation. For example, you can exclude pairs beyond a distance threshold or with invalid orientations.
- *Example*: `--filter_paires_args "--max-dist 2000 --min-dist 100"`. Leave empty (null) to keep default pair filtering behavior.  
- *Type*: string  
- *Default flags:* null  

---

<code style="color:red; font-weight:bold;"> --samtools_markdup_args</code>: SAMtools Markdup Arguments  
- Specifies additional options for the `samtools markdup` command, which marks or removes duplicate reads in BAM files.  
- Use this to control duplicate marking behavior, threading, or reporting options. For example, adding `-r` removes duplicates instead of marking them, and `-@ 8` sets the number of threads.
- *Example*: `--samtools_markdup_args "-r -@ 8"`. Leave empty (null) to use SAMtools default settings.  
- *Type*: string  
- *Default flags:* null  

---

<code style="color:red; font-weight:bold;"> --samtools_sort_2_args</code>: SAMtools Sort (Second Pass) Arguments  
- Specifies additional parameters for the second sorting step performed by `samtools sort`, typically used for name-sorting or coordinate-sorting read pairs.  
- The default `-n` flag sorts alignments by read name, which is often required for pairwise operations. You can modify this to coordinate-sort (`-o`) or add threading options.
- *Example*: `--samtools_sort_2_args "-@ 8 -T tmp -o sorted.bam"`.  
- *Type*: string  
- *Default flags:* -n

---

<code style="color:red; font-weight:bold;"> --bwa_mem_samtools_args</code>: BWA-MEM + SAMtools Combined Arguments  
- Specifies additional command-line options applied jointly to the BWA-MEM alignment output and subsequent SAMtools processing steps.  
- Use this to adjust output conversion and compression parameters after BWA-MEM alignment.
- For example, `-bh` converts to BAM format with the header included. Leave `null` to use the default behavior.  
- *Example*: `--bwa_mem_samtools_args "-bh -@ 8"`  
- *Type*: string  
- *Default flags:* null

---

## Reference & Annotation

Configuration for reference genome files and annotation resources used throughout the pipeline. 

These parameters ensure consistent genome build usage across mapping, feature annotation, and downstream analyses.  

*Help:* Provide correct and consistent genome resources (FASTA index, GTF, chromosome sizes, and blacklist) matching your chosen reference build (e.g., hg38, mm10). Mismatched files can lead to coordinate errors or missing annotations.

---

<code style="color:red; font-weight:bold;"> --jaspar_motif</code>: JASPAR Motif File (TSV)  
- Specifies the URL or local path to the JASPAR motif file (TSV format) used for motif scanning and peak annotation — typically representing transcription factor binding motifs such as CTCF.  
- Provide the path or direct download link to a JASPAR motif file compatible with your reference genome. The default points to the CTCF motif (MA0139.1) for hg38.
- *Example*: `--jaspar_motif http://expdata.cmmt.ubc.ca/JASPAR/downloads/UCSC_tracks/2022/hg38/MA0139.1.tsv.gz`.  
- *Type*: string  
- *Default flags:* http://expdata.cmmt.ubc.ca/JASPAR/downloads/UCSC_tracks/2022/hg38/MA0139.1.tsv.gz  

---

<code style="color:red; font-weight:bold;"> --blacklist</code>: ENCODE Blacklist Regions (BED)  
- Path or URL to a BED file containing ENCODE blacklist regions that should be excluded from peak calling, loop detection, and coverage calculations. These regions are known to produce artificially high signal or mapping artifacts.  
- Use the appropriate blacklist file for your genome build (e.g., hg19, hg38, mm10). The default points to the ENCODE hg38 blacklist (`ENCFF356LFX`).
- *Example*: `--blacklist /refs/hg38-blacklist.bed.gz`.  
- *Type*: string  
- *Default flags:* https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz  

---

<code style="color:red; font-weight:bold;"> --gtf</code>: Gene Annotation File (GTF)  
- Specifies the path or URL to the GTF file containing gene annotations for the reference genome. This file is used to assign peaks, loops, and other genomic features to known genes and transcripts.  
- Provide a GTF file compatible with your chosen reference genome (e.g., Ensembl or GENCODE format). The default points to the GRCh38 annotation from the Illumina iGenomes collection.
- *Example*: `--gtf /data/genomes/GRCh38/genes.gtf`.  
- *Type*: string  
- *Default flags:* s3://ngi-igenomes/igenomes/Homo_sapiens/NCBI/GRCh38/Annotation/Genes/genes.gtf  

---

<code style="color:red; font-weight:bold;"> --chrom_size</code>: Chromosome Sizes File  
- Path to a two-column file listing chromosome names and their corresponding lengths, used to define genomic bounds during binning and matrix construction.  
- This file ensures that cooler and related modules apply consistent chromosome boundaries. It must match the same reference genome build as the input FASTA and GTF files.
- *Example*: `--chrom_size /refs/hg38.chrom.sizes`.  
- *Type*: string  
- *Default flags:* $projectDir/assets/hg38.chrom.sizes  

---

<code style="color:red; font-weight:bold;"> --genomics_features</code>: Genomic Features File (MAPS)  
- Specifies the path or URL to the genomic features file required by MAPS for loop calling. This file contains mappability, GC content, and restriction enzyme fragment information used for bias correction.  
- Use the appropriate MAPS genomic features file matching your genome build and restriction enzyme (e.g., MboI or DpnII). The default points to the hg38 MboI 10 kb resolution file.
- *Example*: `--genomics_features /refs/MAPS_features_hg38_MboI_10kb.txt`.  
- *Type*: string  
- *Default flags:* https://raw.githubusercontent.com/HuMingLab/MAPS/refs/heads/master/MAPS_data_files/hg38/genomic_features/F_GC_M_MboI_10Kb_el.GRCh38.txt

## Peaks & Loops

Configuration of parameters related to peak calling and chromatin loop detection.  

These settings control how enriched regions (peaks) and chromatin interactions (loops) are identified from HiChIP data.  

*Help:* Adjust these options to fine-tune sensitivity and specificity in peak and loop detection.  

For example, you can modify the p-value threshold for MACS3 or enable MAPS loop calling for specific experimental setups.

---

<code style="color:red; font-weight:bold;"> --peak_quality</code>: Peak Significance Threshold  
- Sets the significance cutoff (p-value or q-value, depending on MACS3 configuration) for peak calling to identify enriched regions.  
- Lower values (e.g., 0.01 or 1e-5) increase stringency, reducing false positives but possibly missing weaker peaks. The default of 0.05 is a balanced threshold.
- *Example*: `--peak_quality 0.05`.  
- *Type*: number  
- *Default flags:* 0.05  

---

<code style="color:red; font-weight:bold;"> --genome_size</code>: Genome Size (MACS3)  
- Specifies the genome size shortcut for MACS3 peak calling, corresponding to the total mappable genome length.  
- Use MACS3-supported short codes like `hs` (human), `mm` (mouse), or provide an exact value (e.g., 2.7e9). Must match your reference genome build.
- *Example*: `--genome_size hs`.  
- *Type*: string  
- *Default flags:* hs  

---

<code style="color:red; font-weight:bold;"> --macs3_callpeak_args</code>: MACS3 Additional Arguments  
- Optional custom arguments passed directly to the MACS3 `callpeak` command for advanced configuration.  
- Use this to customize peak calling — for example, enabling model-free mode, adjusting shift/extension sizes, or specifying control input.
- *Example*: `--macs3_callpeak_args "-q 0.01 --nomodel --shift -75 --extsize 150"`.  
- *Type*: string  
- *Default flags:* null  

---

<code style="color:red; font-weight:bold;"> --maps_args</code>: MAPS Loop Calling Arguments  
- Additional command-line flags for the MAPS (Model-based Analysis of PLAC-seq/HiChIP) tool, used for chromatin loop detection.  
- Adjust MAPS behavior such as bin size, distance range, or resolution.
- *Example*: `--maps_args "--bin-size 5000 --cis-only"`. Leave empty (null) to use the default MAPS settings.  
- *Type*: string  
- *Default flags:* null  

---

<code style="color:red; font-weight:bold;"> --skip_maps</code>: Skip MAPS Module  
- Determines whether to skip the MAPS loop calling module during the pipeline execution.  
- Set this to `true` to disable loop calling (useful for QC or peak-only runs). Set to `false` to perform full MAPS-based loop analysis.
- *Example*: `--skip_maps false`.  
- *Type*: boolean  
- *Default flags:* true

## Contact Matrices & Binning

Settings controlling the generation, binning, and normalization of contact matrices from HiChIP data.  

These parameters determine how raw read pairs are converted into multi-resolution .cool files for downstream analysis.  

*Help:* Adjust these options to tune the resolution and structure of the resulting chromatin contact maps.  

For example, smaller bin sizes (e.g., 1 kb) provide higher resolution but require more memory and processing time, while larger bins (e.g., 10 kb or 25 kb) yield smoother matrices for large-scale analyses.

---

<code style="color:red; font-weight:bold;"> --cool_bin</code>: Base Bin Size (bp)
- Defines the base bin size, in base pairs, used for generating contact matrices in the .cool format.  
- Smaller bin sizes (e.g., 1000 bp) provide higher resolution but increase memory and runtime requirements.
- Larger bins (e.g., 5000–25000 bp) are recommended for lower-depth datasets.
- *Example*: `--cool_bin 5000`.  
- *Type*: integer  
- *Default flags:* 1000  

---

<code style="color:red; font-weight:bold;"> --cooler_cload_args</code>: Cooler Cload Arguments  
- Specifies command-line options passed to the `cooler cload` command for loading read pairs into a .cool file.  
- Use this to modify how columns from the pairs file are interpreted when generating contact matrices. The default is configured for standard pairtools output.
- *Example*: `--cooler_cload_args "pairs --zero-based -c1 1 -p1 2 -c2 3 -p2 4"`.  
- *Type*: string  
- *Default flags:* pairs --zero-based  -c1 2 -p1 3 -c2 4 -p2 5  

---

<code style="color:red; font-weight:bold;"> --cooler_zoomify_res</code>: Zoomify Resolution Preset  
- Specifies the resolution preset key used for multi-resolution cooler files (`cooler zoomify`).
- Each key corresponds to a predefined set of bin sizes defined in `insulation_resultions`.  
- Use this to select which resolution set to apply during matrix zoomification. For example, `1000N` corresponds to 1 kb–500 kb bins by default.
- *Example*: `--cooler_zoomify_res 1000N`.  
- *Type*: string  
- *Default flags:* 1000N  

---

<code style="color:red; font-weight:bold;"> --cooler_zoomify_args</code>: Cooler Zoomify Additional Arguments  
- Extra command-line options for the `cooler zoomify` command that generate multi-resolution cooler files.  
- Use this to control balancing, overwrite behavior, or other Zoomify parameters.
- *Example*: `--cooler_zoomify_args "--balance --force"`. Leave empty ("") to use defaults.  
- *Type*: string  
- *Default flags:* null  

---

<code style="color:red; font-weight:bold;"> --insulation_resultions</code>: Insulation Score Resolutions  
- Specifies the resolution presets (in base pairs) used for insulation score and TAD boundary calculations.
- Each key defines a named preset containing multiple window sizes.  
- Provide one or more resolution sets, where each key (e.g., 1000N) maps to a list of window sizes that correspond to the resolutions used for cooler zoomify. These determine the granularity of TAD detection and insulation profiling..
- *Example*:  
  ```
  --insulation_resultions '{"1000N": "1000 2000 5000 10000 20000 50000 100000 200000 500000"}'
  ```  
- *Type*: object  
- *Default flags:* {"1000N": "1000 2000 5000 10000 20000 50000 100000 200000 500000"}  

---

<code style="color:red; font-weight:bold;">--cooltools_eigscis_args</code>: Cooltools Eigs-Cis Arguments  
- Additional command-line flags for the `cooltools eigs-cis` command, which computes eigenvectors for A/B compartment analysis.  
- Modify this to control the number of eigenvectors or specify contact type (cis/trans).
- *Example*: `--cooltools_eigscis_args "--n-eigs 2 --contact-type cis"`.  
- *Type*: string  
- *Default flags:* --n-eigs 1  

---

<code style="color:red; font-weight:bold;">--cooler_eigscis_resultion</code>: Eigenvector Resolution (bp)  
- Specifies the resolution, in base pairs, used for eigenvector computation in A/B compartment analysis.  
- Smaller values provide higher resolution but require denser contact maps.
- *Example*: `--cooler_eigscis_resultion 100000`.  
- *Type*: integer  
- *Default flags:* 100000  

---

<code style="color:red; font-weight:bold;">--calder_bin</code>: CALDER Bin Size  
- Sets the bin size used for subcompartment calling by CALDER.  
- Provide bin size as a number (e.g., 10000) or scientific notation (e.g., `10E3`). Higher bin sizes reduce resolution but speed up computation.
- *Example*: `--calder_bin 25000`.  
- *Type*: string  
- *Default flags:* 100000  

---

<code style="color:red; font-weight:bold;">--gstripe_args</code>: gStripe Detection Arguments  
- Extra command-line options for the `gStripe` tool, which detects stripe-like chromatin interaction patterns.  
- Modify detection thresholds or bin correction behavior. The default `--fix_bin_start` ensures stripe detection starts from bin-aligned positions.
- *Example*: `--gstripe_args "--fix_bin_start --minlen 3 --qval 0.05"`.  
- *Type*: string  
- *Default flags:* --fix_bin_start  

## Visualization & QC

Configuration for plotting and quality-control steps (e.g., FastQC, deepTools, Juicer Tools, pairtools).  

These options affect coverage plots, sample correlations, and other diagnostic visualizations.  

*Help:* Tune these to control plot styles, correlation methods, and output formats.  

For large cohorts, consider lighter settings (e.g., skipping numbers on heatmaps) to speed up rendering.

---

<code style="color:red; font-weight:bold;"> --plot_method</code>: Correlation Method for Plots  
- Defines the statistical method used to compute sample-to-sample correlations in deepTools visualizations.  
- Choose between `spearman` (rank-based, robust to outliers) or `pearson` (linear correlation).
- The default `spearman` is typically used for read-count correlation heatmaps.
- *Example*: `--plot_method pearson`.  
- *Type*: string  
- *Default flags:* spearman  

---

<code style="color:red; font-weight:bold;"> --plot_type</code>: Plot Type for Correlation Visualization  
- Specifies the visualization format for correlation results — typically a heatmap or scatter plot.  
- Select the style best suited to your comparison. `heatmap` provides an overview of pairwise correlations, while `scatter` highlights relationships between specific samples.
- *Example*: `--plot_type scatter`.  
- *Type*: string  
- *Default flags:* heatmap  

---

<code style="color:red; font-weight:bold;"> --fastqc_args</code>: FastQC Arguments  
- Specifies additional command-line options for FastQC quality control analysis of raw sequencing reads.  
- Use this to control FastQC behavior, such as verbosity or output compression. The default `--quiet` suppresses detailed logging.
- *Example*: `--fastqc_args "--quiet --nogroup"`.  
- *Type*: string  
- *Default flags:* --quiet  

---

<code style="color:red; font-weight:bold;"> --deeptools_plotcoverage_args</code>: deepTools PlotCoverage Arguments  
- Additional flags for the deepTools `plotCoverage` module, which visualizes genome-wide read coverage across samples.  
- Customize output options such as file format, scaling, or whether to skip empty bins. The default `--skipZeros` ignores regions with zero coverage. Example: `--deeptools_plotcoverage_args "--skipZeros --plotFileFormat pdf"`.  
- *Type*: string  
- *Default flags:* --skipZeros  

---

<code style="color:red; font-weight:bold;"> --deeptools_plotcorrelation_args</code>: deepTools PlotCorrelation Arguments  
- Custom flags for the deepTools `plotCorrelation` command, which computes and visualizes sample correlation matrices.  
- Use this to modify correlation type, color scheme, and plot annotations. The default parameters produce a Spearman correlation heatmap with labeled values.
- *Example*: `--deeptools_plotcorrelation_args "--corMethod pearson --colorMap RdBu --what pairwise"`.  
- *Type*: string  
- *Default flags:* --skipZeros --plotTitle "Spearman Correlation of Read Counts"  --colorMap RdYlBu --plotNumbers  

---

<code style="color:red; font-weight:bold;"> --juicertools_args</code>: Juicer Tools Arguments  
- Specifies extra command-line parameters for Juicer Tools utilities used for Hi-C/HiChIP contact map analysis.  
- Set parameters for operations like normalization, balancing, or matrix extraction. Leave null to use the default Juicer settings.
- *Example*: `--juicertools_args "pre --threads 8 --resolutions 5000"`  
- *Type*: string  
- *Default flags:* null  

---

<code style="color:red; font-weight:bold;"> --pairtools_parse2_args</code>: Pairtools Parse2 Arguments  
- Additional options for the `pairtools parse2` command used to convert alignment files into a standardized pairs format.  
- Modify column selection, filtering, or threading for parsing steps. Leave null for default settings.  
- *Example*: `--pairtools_parse2_args "--add-columns chr1,chr2,pos1,pos2 --nproc 8"`. 
- *Type*: string  
- *Default flags:* null

 
## 3D Genome Modelling (MultiMM)

Configuration options for 3D genome reconstruction and visualization using the MultiMM module.  

These parameters define which genomic regions are modeled, the computational platform used (CPU or GPU), and optional fine-grained settings for chromosomal or locus-level simulations.  

*Help:* Adjust these parameters to specify the modeling scope (whole chromosome, gene region, or specific locus) and computational mode.  

MultiMM integrates chromatin contact data to predict spatial genome structures. For large models, prefer GPU acceleration if available.

---

<code style="color:red; font-weight:bold;">--multimm_platform</code>: MultiMM Computational Platform  
- Specifies the compute platform used for MultiMM 3D genome modeling — either CPU or GPU-based execution, depending on system availability and dataset size.  
- The following flags can be used:  
  - `CPU` — standard desktop or cluster execution (default).  
  - `OpenCL` — enables GPU acceleration using the OpenCL framework for compatible hardware.  
  - `CUDA` — enables GPU acceleration using NVIDIA CUDA, providing faster simulations for large or high-resolution models.  
- *Example*: `--multimm_platform CUDA`  
- *Type*: string  
- *Default flags:* CPU   

---

<code style="color:red; font-weight:bold;"> --multimm_modelling_level</code>: Modelling Granularity Level  
- Specifies the genomic scale at which MultiMM performs 3D genome modeling — ranging from a single gene to the entire genome.  
- Choose the modeling depth according to your analysis goal and computational resources.
- The following modelling levels are available:
  - **GENE** — Provide a gene of interest (with an associated .bedpe file path). MultiMM models the gene with a default ±100 kb window around it. Compartment forces are *not* considered at this level.
  - **REGION** — Specify a chromosome and genomic coordinates (start–end). Compartment interactions can optionally be included. Only the selected region is modeled.
  - **CHROMOSOME** — Provide a chromosome name; MultiMM automatically determines start and end coordinates. Compartments of data have to be imported.
  - **GW (Genome-Wide)** — Models the entire genome. No input for chromosome or coordinates is needed. Compartments of data have to be imported. This is the most computationally intensive mode, potentially taking minutes to hours depending on system performance.  
- *Example*: `--multimm_modelling_level region`  
- *Type*: string  
- *Default flags:* chrom  

---

<code style="color:red; font-weight:bold;"> --multimm_gene_name</code>: Gene Name for Modelling  
- Specifies the gene symbol or identifier to be modeled in 3D when the modelling level is set to `gene`.  
- Provide a valid gene symbol that exists within your reference annotation (e.g., `PIRAT1`, `RUNX1`). This parameter is ignored if the modelling level is set to `chrom` or `region`.
- *Example*: `--multimm_gene_name RUNX1`.  
- *Type*: string  
- *Default flags:* null  

---

<code style="color:red; font-weight:bold;"> --multimm_chrom</code>: Chromosome to Model  
- Defines which chromosome should be used for 3D genome modeling. Required for both *chromosome* level and *region* level modeling.  
- Use standard chromosome naming consistent with your reference genome (e.g., `chr1`, `chr21`, `chrX`).
- *Example*: `--multimm_chrom chr10`.  
- *Type*: string  
- *Default flags:* chr21  

---

<code style="color:red; font-weight:bold;"> --multimm_loc_start</code>: Locus Start Coordinate (bp)  
- Specifies the genomic start coordinate (in base pairs) for *region* level 3D modeling.  
- Used only when the modelling level is set to `locus`. Must be a valid coordinate within the selected chromosome.
- *Example*: `--multimm_loc_start 28000000`.  
- *Type*: integer  
- *Default flags:* null  

---

<code style="color:red; font-weight:bold;"> --multimm_loc_end</code>: Locus End Coordinate (bp)  
- Specifies the genomic end coordinate (in base pairs) for *region* level 3D modeling.  
- Used together with `multimm_loc_start` to define the genomic window for 3D reconstruction.
- *Example*: `--multimm_loc_end 32000000`.  
- *Type*: integer  
- *Default flags:* null  

---

<code style="color:red; font-weight:bold;"> --multimm_args</code>: MultiMM Additional Arguments  
- Custom command-line arguments passed directly to MultiMM for fine-tuning modeling behavior and output.  
- Use this for advanced control — e.g., number of models, iteration limits, or output directory.
- *Example*: `--multimm_args "--n-models 50 --max-iters 2000"`. Leave empty (`null`) to use defaults.  
- *Type*: string  
- *Default flags:* null  


---

## Quick Copy‑Paste: params.config template

```groovy
params {

  /* ========= GENERAL OPTIONS ========= */
  ref_short = "hg38"
  outdir    = "/mnt/raid/test_case/results"
  threads   = 8
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
  cooler_zoomify_res         = "1000N"
  cooler_zoomify_args        = null
  insulation_resultions      = '{"1000N": "1000 2000 5000 10000 20000 50000 100000 200000 500000"}'
  cooltools_eigscis_args     = "--n-eigs 1"
  cooler_eigscis_resultion   = 100000
  calder_bin                 = "10E3"
  gstripe_args               = "--fix_bin_start"

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
