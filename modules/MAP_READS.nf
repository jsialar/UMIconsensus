process MERGE_FASTQ {
    //publishDir "${params.output}/${samplename}/merged_fastq", mode: 'copy', pattern: "*.fastq"
    
    input:
        tuple val(samplename), path(fastqs)
    
    output:
        tuple val( "${samplename}" ), path( "${samplename}.merged.fastq.gz"), emit: merged_fastq
    
    script:
    """
        cat ${fastqs} > ${samplename}.merged.fastq.gz
        
    """
}


process TRIM_FASTQ {
    //publishDir "${params.output}/${samplename}/trimmed_fastq", mode: 'copy'    
    
    input:
        tuple val(samplename), path(fastq)
    
    output:
        tuple val( "${samplename}" ), path("${samplename}.trimmed.fastq.gz"),  optional: true, emit: trimmed_fastq
        tuple val( "${samplename}" ), path("${samplename}.json"), emit: cutadapt_info
        tuple val( "${samplename}" ), path("${samplename}.untrimmed.fastq.gz"), optional: true, emit: untrimmed_fastq
    
    script:
    def adapter_full="^AATGATACGGCGACCACCGAGATCTACACNNNNNNNNACACTCTTTCCCTACACGACGCTCTTCCGATCT...AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG\$"
    def adapter_short="ACACTCTTTCCCTACACGACGCTCTTCCGATCT...AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
    def save_untrimmed = params.save_untrimmed ? "--untrimmed-output ${samplename}.untrimmed.fastq.gz" : "--discard-untrimmed"
    """
        cutadapt -g F=${adapter_full} \
            -e 0.15 \
            --revcomp \
            --cores ${params.threads} \
            --minimum-length ${params.minlength} \
            --maximum-length ${params.maxlength} \
            --json ${samplename}.json \
            ${save_untrimmed} \
            ${fastq} | \
        cutadapt --cut 8 \
            --cut -8 \
            --cores ${params.threads} \
            --rename '{id}:{cut_prefix}{cut_suffix}' \
            -o ${samplename}.trimmed.fastq.gz - 
    """
}


process MAP_READS {
    //publishDir "${params.output}/${sample}/align/", mode: 'copy'

    input:
        tuple val( samplename ), path( trimmed_fastq )
        path reference
    output:
        tuple val( "${samplename}"), path ( "${samplename}.bam" ), path ( "${samplename}.bam.bai" ), emit: bamfile

    script:
    """
        minimap2 \
          ${params.minimap2_param} \
          -t ${params.threads} \
          ${reference} \
          ${trimmed_fastq} | 
        samtools sort \
          -@ ${params.threads} \
          -o ${samplename}.bam - && \
        samtools index \
          -@ ${params.threads} \
          ${samplename}.bam
    """
}