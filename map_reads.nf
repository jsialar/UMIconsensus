// nextflow.preview.output = true

include { MERGE_FASTQ; TRIM_FASTQ; MAP_READS } from './modules/MAP_READS.nf'
include { GENERATE_CONSENSUS } from './workflows/GENERATE_CONSENSUS.nf'
include { publish_artifact } from './modules/common.nf'

workflow {

    Channel
    .fromPath("${params.input}/barcode*/*.fastq.gz")
    .map{ 
        fastqs -> 
        def barcode = fastqs.parent.name
        tuple(barcode, fastqs)
        }
    .groupTuple( by: 0 ) 
    .set{ existing_fastqs }
    
    existing_fastqs.count().set{sample_count}
    //sample_count.view()

    MERGE_FASTQ(existing_fastqs)
    TRIM_FASTQ(MERGE_FASTQ.out.merged_fastq, sample_count)
    MAP_READS(TRIM_FASTQ.out.trimmed_fastq,  file("${params.reference}", checkIfExists: true), sample_count)

    TRIM_FASTQ.out.cutadapt_info
    .combine(["trim"])
    .mix(
        TRIM_FASTQ.out.untrimmed_fastq
            .combine(["trim"]),

        MAP_READS.out.bamonly
            .combine(["aligned"]),

        MAP_READS.out.bamindexonly
            .combine(["aligned"]),
   
    )
    .set{output_ch}


    publish_artifact(output_ch) 

}
