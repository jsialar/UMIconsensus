process CALL_CONSENSUS {
    //publishDir "${params.output}/${sample}/align_consensus", mode: 'copy'    
    
    input:
        tuple val( sample ), val(target), path ( bam ), path ( bam_bai )        
        
        output:
        tuple val( "${sample}" ), val("${target}"), path("*.fastq.gz"), emit: consensus_fastq, optional:true
        tuple val( "${sample}" ), val("${target}"), path("*.fastq.gz"), path("*.csv"), emit: fastq_w_table
        tuple val( "${sample}" ), path("*.csv"), emit: consensus_table
    
    script:
    def keepundetermined = params.keep_nonconsensus ? "--keep" : ""
    """
    call_consensus.py --inputpath ./ \
    --nonpolysnp ${params.snplist} \
    --samplename ${sample} \
    ${keepundetermined}
    """
}