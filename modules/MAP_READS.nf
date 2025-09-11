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
        val samplecount
    
    output:
        tuple val( "${samplename}" ), path("${samplename}.trimmed.fastq.gz"),  optional: true, emit: trimmed_fastq
        tuple val( "${samplename}" ), path("${samplename}.json"), emit: cutadapt_info
        tuple val( "${samplename}" ), path("${samplename}.untrimmed.fastq.gz"), optional: true, emit: untrimmed_fastq
    
    script:
    def adapter_full="^AATGATACGGCGACCACCGAGATCTACACNNNNNNNNACACTCTTTCCCTACACGACGCTCTTCCGATCT...AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG\$"
    def adapter_short="ACACTCTTTCCCTACACGACGCTCTTCCGATCT...AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
    def save_untrimmed = params.save_untrimmed ? "--untrimmed-output ${samplename}.untrimmed.fastq.gz" : "--discard-untrimmed"
    def numberOfCores = Runtime.getRuntime().availableProcessors()
    def threads = (numberOfCores / samplecount) as int 
    """
        cutadapt -g F=${adapter_full} \
            -e 0.15 \
            --revcomp \
            --cores ${threads} \
            --minimum-length ${params.minlength} \
            --maximum-length ${params.maxlength} \
            --json ${samplename}.json \
            ${save_untrimmed} \
            ${fastq} | \
        cutadapt --cut 8 \
            --cut -8 \
            --cores ${threads} \
            --rename '{id}:{cut_prefix}{cut_suffix}' \
            -o ${samplename}.trimmed.fastq.gz - 
    """
}


process MAP_READS {
    //publishDir "${params.output}/${sample}/align/", mode: 'copy'
    maxForks 4    

    input:
        tuple val( samplename ), path( trimmed_fastq )
        path reference
        val samplecount

    output:
        tuple val( "${samplename}"), path ( "${samplename}.bam" ),  path ( "${samplename}.bam.bai" ), emit: bamfile        
        tuple val( "${samplename}"), path ( "${samplename}.bam" ), emit: bamonly
        tuple val( "${samplename}"), path ( "${samplename}.bam.bai" ), emit: bamindexonly

    script:
    def numberOfCores = Runtime.getRuntime().availableProcessors()
    def forks = samplecount < 4 ? samplecount : 4
    def threads = (numberOfCores / forks) as int
    """
        minimap2 \
          ${params.minimap2_param} \
          -t ${threads} \
          ${reference} \
          ${trimmed_fastq} | 
        samtools sort \
          -@ ${threads} \
          -o ${samplename}.bam - && \
        samtools index \
          -@ ${threads} \
          ${samplename}.bam
    """
}