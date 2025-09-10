process PROCESS_CONSENSUS {
   
    input:
        tuple val( sample ), val(target), path(fastq), path(consensustable)
    
    output:
        tuple val( "${sample}" ), path("qual.csv"), emit: qual
        tuple val( "${sample}" ), path("hap.csv"), emit: hap
    
    script:
    """
    process_consensus.py \
    --sample ${sample} \
    --target ${target} \
    --fastqpath ${fastq} \
    --tablepath ${consensustable}    
    """
}