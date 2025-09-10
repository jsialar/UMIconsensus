process SPLIT_GROUP {
    //publishDir "${params.output}/${sample}/align_split/", mode: 'copy'    
    
    input:
        tuple val( sample ), val( target ),path (input) 
        
    
    output:
        tuple val( "${sample}" ), val("${target}"), path("*.bam"), path("*.bai"), optional: true,  emit: splitted_bam        
    
    script:

    """
    samtools split \
    -d UG \
    -M -1 \
    ${input}
    samtools index -M P*.bam
    """
}
