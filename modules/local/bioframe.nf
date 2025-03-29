process BIOFRAME {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::bioframe=0.7.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioframe:0.7.2--pyhdfd78af_0' :
        'biocontainers/bioframe:0.7.2--pyhdfd78af_0' }"

    input:
    tuple val(meta) , path(peaks)
    path(jaspar_motif)
    path (blacklist)

    output:
    tuple val(meta), path("*_peak_with_CTCFmotif.tsv"),  emit: tsv
    path  "versions.yml",            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    #!/usr/bin/env python3
    # taken from https://bioframe.readthedocs.io/en/latest/tutorials/tutorial_assign_motifs_to_peaks.html
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import bioframe
    import sys

    #base_dir = "/tmp/bioframe_tutorial_data/"
    #assembly = "GRCh38"
    ctcf_peaks = bioframe.read_table("${peaks}",schema="narrowPeak",)
    first_col = ctcf_peaks.columns[0]
    if not str(ctcf_peaks[first_col].iloc[0]).startswith("chr"):
        ctcf_peaks[first_col] = "chr" + ctcf_peaks[first_col].astype(str)

    ctcf_motifs = bioframe.read_table("${jaspar_motif}", schema="jaspar", skiprows=1)

    ### Overlap peaks & motifs
    df_peaks_motifs = bioframe.overlap(ctcf_peaks, ctcf_motifs, suffixes=("_1", "_2"), return_index=True)
    df_peaks_motifs["pval_2"] = df_peaks_motifs["pval_2"].fillna(-1)
    idxmax_peaks_motifs = (df_peaks_motifs.groupby(["chrom_1", "start_1", "end_1"])["pval_2"].idxmax().values)
    df_peaks_maxmotif = df_peaks_motifs.loc[idxmax_peaks_motifs]
    df_peaks_maxmotif["pval_2"].replace(-1, np.nan, inplace=True)

    ### filter peaks overlapping blacklisted regions for hg38
    blacklist = bioframe.read_table("${blacklist}", schema="bed3",)
    closest_to_blacklist = bioframe.closest(ctcf_peaks, blacklist)

    # first let's select the columns we want for our final dataframe of peaks
    # with motifs
    df_peaks_maxmotif = df_peaks_maxmotif[
        [
            "chrom_1",
            "start_1",
            "end_1",
            "fc_1",
            "chrom_2",
            "start_2",
            "end_2",
            "pval_2",
            "strand_2",
        ]
    ]


    #  then rename columns for convenience when subtracting
    for i in df_peaks_maxmotif.keys():
        if "_1" in i:
            df_peaks_maxmotif.rename(columns={i: i.split("_")[0]}, inplace=True)

    # now subtract, expanding the blacklist by 1kb
    df_peaks_maxmotif_clean = bioframe.subtract(df_peaks_maxmotif, bioframe.expand(blacklist, 1000))

    # dataframe containing positions of CTCF ChIP peaks, including the strongest motif underlying that peak, and after conservative filtering for proximity to blacklisted regions

    df_peaks_maxmotif_clean.to_csv("${prefix}_peak_with_CTCFmotif.tsv", sep="\\t", index=False)

    with open ("versions.yml", "w") as version_file:
	    version_file.write("\\"${task.process}\\":\\n    python: {}\\n".format(sys.version.split()[0].strip()))

    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    touch ${prefix}_peak_with_CTCFmotif.tsv
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
