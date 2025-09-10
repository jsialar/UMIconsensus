// nextflow.preview.output = true

include { MERGE_FASTQ; TRIM_FASTQ; MAP_READS } from './modules/MAP_READS.nf'
include { GENERATE_CONSENSUS } from './workflows/GENERATE_CONSENSUS.nf'

process publish_artifact {
    publishDir "${params.output}/${sample}/${module}", mode: 'copy', pattern: "*"
    input:
        tuple val(sample), file(outfile), val(module)
        
    output:
        file outfile

    script:
    """
    echo "Writing output files"
    """
    }

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

    MERGE_FASTQ(existing_fastqs)
    TRIM_FASTQ(MERGE_FASTQ.out.merged_fastq)
    MAP_READS(TRIM_FASTQ.out.trimmed_fastq,  file("${params.reference}", checkIfExists: true))

    GENERATE_CONSENSUS(MAP_READS.out.bamfile)

    TRIM_FASTQ.out.cutadapt_info
    .combine(["trim"])
    .mix(
        TRIM_FASTQ.out.untrimmed_fastq
            .combine(["trim"]),

        MAP_READS.out.bamonly
            .combine(["aligned"]),

        MAP_READS.out.bamindexonly
            .combine(["aligned"]),
        
        GENERATE_CONSENSUS.out.umi_sizes_ch
            .map{file -> tuple(file.getSimpleName().tokenize('_')[0], file)}
            .combine(["group_umi"]),

        GENERATE_CONSENSUS.out.consensustable_ch
            .map{file -> tuple(file.getSimpleName().tokenize('_')[0], file)}
            .combine(["call_consensus"]),

        GENERATE_CONSENSUS.out.consensus_quals
            .map{file -> tuple(file.getSimpleName().tokenize('_')[0], file)}
            .combine(["call_consensus"]),

        GENERATE_CONSENSUS.out.consensus_haps
            .map{file -> tuple(file.getSimpleName().tokenize('_')[0], file)}
            .combine(["call_consensus"])
    )
    .set{output_ch}


    publish_artifact(output_ch)
 

    // publish:
    // cutadapt_info = TRIM_FASTQ.out.cutadapt_info
    // untrimmed_fastq = TRIM_FASTQ.out.untrimmed_fastq
    // mapped_reads = MAP_READS.out.bamfile
    // umi_sizes = GENERATE_CONSENSUS.out.umi_sizes_ch
    // consensustable = GENERATE_CONSENSUS.out.consensustable_ch
    // consensus_fastq = GENERATE_CONSENSUS.out.consensus_fastq
    // consensus_quals = GENERATE_CONSENSUS.out.consensus_quals
    // consensus_haps = GENERATE_CONSENSUS.out.consensus_haps


}

// output {

//     cutadapt_info {
//         mode 'copy'
//         path { item -> 
//             def sample = item[0]
//             "${params.output}/${sample}/trim/"
//         }
//     }

//     untrimmed_fastq {
//         mode 'copy'
//         path { item -> 
//             def sample = item[0]
//             "${params.output}/${sample}/trim/"
//         }
//     }

//     mapped_reads {
//         mode 'copy'
//         path { item -> 
//             def sample = item[0]
//             "${params.output}/${sample}/aligned/"
//         }
//     }

//     umi_sizes {
//         mode 'copy'
//         path { file -> 
//             def sample = file.getSimpleName().tokenize('_')[0]
//             "${params.output}/${sample}/group_umi/"
//         }
//     }

//     consensustable {
//         mode 'copy'
//         path { file -> 
//             def sample = file.getSimpleName().tokenize('_')[0]
//             "${params.output}/${sample}/call_consensus/"
//         }
//     }

//     consensus_fastq {
//         mode 'copy'
//         path { item -> 
//             def sample = item[0]
//             "${params.output}/${sample}/call_consensus/"
//         }
//     }

//     consensus_quals {
//         mode 'copy'
//         path { file -> 
//             def sample = file.getSimpleName().tokenize('_')[0]
//             "${params.output}/${sample}/call_consensus/"
//         }
//     }

//     consensus_haps {
//         mode 'copy'
//         path { file -> 
//             def sample = file.getSimpleName().tokenize('_')[0]
//             "${params.output}/${sample}/call_consensus/"
//         }
//     }
// }