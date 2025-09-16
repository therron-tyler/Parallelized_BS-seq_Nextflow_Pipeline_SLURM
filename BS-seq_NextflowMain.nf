nextflow.enable.dsl=2

workflow {
    reads_ch       = Channel.fromPath(params.input_pattern)
    metrics_ch     = MethylSeq(reads_ch)

    metrics_all_ch = metrics_ch.collect()
    aggregate_ch   = AggregateMetrics(metrics_all_ch)
    beta_ch        = CalcBeta(metrics_all_ch)
    dss_ch         = Cov2DSS(metrics_all_ch)
    bsseq_ch       = ConvertBSseq(dss_ch)
    dml_dmr_ch     = CallDMLsDMRs(bsseq_ch)
    flat_dmr_ch    = dml_dmr_ch.flatten()
    sig_dmr_ch     = flat_dmr_ch.filter { file -> file.name.endsWith('_significant_DMRs.csv') }
    annotated_ch   = AnnotateDMRs(sig_dmr_ch)
    Visualize(beta_ch)
}

process MethylSeq {
    publishDir "Results/BSseq_Sample_Processing", mode: 'copy', pattern: 'metrics_*'
    tag { read.baseName.replaceFirst(/_L\d+.*$/, '') }
    input:
        path read
    output:
        path "metrics_${read.baseName.replaceFirst(/_L\\d+.*$/, '')}", emit: metrics_ch
    script:
    """
    module load ${params.moduleLoad}
    SAMPLE=${read.baseName.replaceFirst(/_L\\d+.*$/, '')}
    bash ${params.pipeline_script} \
         -i \$SAMPLE \
         -1 ${read} \
         -g ${params.genome_dir} \
         -o metrics_\${SAMPLE}
    """
}

process AggregateMetrics {
    publishDir "Results/Aggregate_QC_Metrics", mode: 'copy'
    input:
        path metrics_dirs
    output:
        path 'aggregate_metrics.tsv', emit: aggregate_ch
    script:
    """
    module load python/3.8.4
    python ${params.aggregate_script} \
           --input ${metrics_dirs.join(' ')} \
           --output aggregate_metrics.tsv
    """
}

process CalcBeta {
    publishDir "Results/BetaScore_Calculation", mode: 'copy'
    input:
        path metrics_dirs
    output:
        path 'per_cpg_methylation_matrix.tsv', emit: beta_ch
    script:
    """
    module load python/3.8.4
    python ${params.beta_script}
    """
}

process Cov2DSS {
    publishDir "Results/DSS_Input", mode: 'copy'
    tag 'cov2dss'
    input:
        path metrics_dirs
    output:
        path "${params.dss_outdir}", emit: dss_ch
    script:
    """
    module load python/3.8.4
    mkdir -p ${params.dss_outdir}
    for dir in ${metrics_dirs.join(' ')}; do
    python ${params.cov2dss_script} \
           -p "\$dir/**/*.cov*" \
           -o ${params.dss_outdir}
    done
    """
}

process ConvertBSseq {
    publishDir "Results/BSseq", mode: 'copy'
    input:
        path dss_dir
    output:
        path "${params.bsseq_outfile}", emit: bsseq_ch
    script:
    """
    module load R/4.2.0
    Rscript ${params.bsseq_script} \
           --input_dir ${dss_dir} \
           --output ${params.bsseq_outfile}
    """
}

process CallDMLsDMRs {
    publishDir "Results/DMR_Calling", mode: 'copy'
    input:
        path bsseq_file
    output:
        path "${params.dml_out_prefix}*.csv", emit: dml_dmr_ch
 
    script:
    """
    module load R/4.2.0
    Rscript ${params.dml_dmr_script} \
           --groups ${params.groups_csv} \
           --grp1 ${params.grp1} \
           --grp2 ${params.grp2} \
           --bsobj ${bsseq_file} \
           --out ${params.dml_out_prefix}
    """
}

process AnnotateDMRs {
    publishDir "Results/GeneAnnotation_to_DMR", mode: 'copy'
    tag { dmr_csv.baseName }
    input:
        path dmr_csv
    output:
        path "${params.annotation_outfile}.csv", emit: annotated_ch
    script:
    """
    module load R/4.2.0
    Rscript ${params.annotation_script} \
           --dmr_csv ${dmr_csv} \
           --ref_seq ${params.RefSeq} \
           --out ${params.annotation_outfile}
    """
}

process Visualize {
    publishDir "Results/Bscore_Visualization", mode: 'copy'
    input:
        path matrix_tsv
    output:
        path 'bscore_distribution.png'
    script:
    """
    module load python/3.8.4
    python ${params.viz_script} \
           $matrix_tsv \
           --cdf \
           --density \
           --bins ${params.viz_bins} \
           --out bscore_distribution.png
    """
}
