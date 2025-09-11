process PARSE_BED {
   
    input:
    tuple val( target ), val( line )

    output:
    tuple val( target ), path("*.bed"), emit: bed_channel

    script:
    """
    echo '$line' > ${target}.bed
    """
}

process GROUP_UMI {
   
    input:
        tuple val( sample ), path (input_bam) , path (input_bai), val(target), path(bed)
        
    
    output:
        tuple val( "${sample}" ), val("${target}"), path("${target}.bam"), emit: grouped_bam, optional: true      
        tuple val( "${sample}" ), path("*.tsv"), emit: umi_sizes, optional: true
    
    script:

    """
    group.py -I ${input_bam} \
    --size-out ${sample}_sizedistribution.tsv \
    --output-bam \
    --log=group.log \
    --bed=${bed} \
    --sizethreshold=${params.familysizethreshold} \
    --umi-separator=":" \
    --method=${params.method}
    """
}