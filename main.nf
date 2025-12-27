
process DownloadCounts {

    publishDir "${params.outdir}/sample_counts", mode: 'copy'

    input:
    val gse_id 

    output:
        path "expression_data.pkl"
        path "samples.tsv"

    script:
    """
    python3 "${projectDir}/scripts/download_counts.py" --gse_id ${gse_id}
    """
}


process BuildMatrices {
    publishDir "${params.outdir}/matrices_TPM_CPM_meta", mode: 'copy'

    input:
        path expression_data
        path samples

    output:
        path "counts.tsv"
        path "cpm.tsv"
        path "tpm.tsv"
        path "meta.tsv"

    script:
    """
    python3 "${projectDir}/scripts/build_matrices.py" "${expression_data}" "${samples}"
    """
}


process QualityControlPlots {

    publishDir "${params.outdir}/qc_plots", mode: 'copy'

    input:
        path counts_tsv
        path cpm_tsv
        path tpm_tsv
        path meta_tsv

    output:
        path "*.png"
        path "*.pdf"

    script:
    """
    python3 "${projectDir}/scripts/qc_plots.py" \\
           "${counts_tsv}" \\
           "${cpm_tsv}" \\
           "${tpm_tsv}" \\
           "${meta_tsv}"
    """
}

process DifferentialExpression {
    publishDir "${params.outdir}/differential_expression", mode: 'copy'

    input:
        path counts_tsv
        path meta_tsv

    output:
        path "deseq_results.csv"
        path "top_genes.csv"
        path "volcano_plot.png"
        path "volcano_plot.pdf"

    script:
    """
    python3 "${projectDir}/scripts/differential_expression.py" \\
           "${counts_tsv}" \\
           "${meta_tsv}"
    """
}

workflow {

    ch_gse = Channel.value(params.gse_id)

    ch_download = DownloadCounts(ch_gse)

    ch_matrices = BuildMatrices(ch_download)

    QualityControlPlots(ch_matrices)

    DifferentialExpression(
    ch_matrices[0],   
    ch_matrices[3]   
)
}
