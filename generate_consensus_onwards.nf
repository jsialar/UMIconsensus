// nextflow.preview.output = true
// Removed: Nextflow 25.10 no longer accepts this assignment in a script file
// ("No such variable: output"). This workflow publishes via publish_artifact
// and has no output {} block, so the flag was never doing anything. main.nf
// already has the same line commented out.

include { GENERATE_CONSENSUS } from './workflows/GENERATE_CONSENSUS.nf'
include { publish_artifact } from './modules/common.nf'


workflow {

    Channel.fromFilePairs("${params.bam}{,.bai}", flat: true)
    .set{bam_ch}

    GENERATE_CONSENSUS(bam_ch)

    GENERATE_CONSENSUS.out.umi_sizes_ch
        .map{file -> tuple(file.getSimpleName().tokenize('_')[0], file)}
        .combine(["group_umi"])
        .mix(

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



    if (params.additional_outputs) {
        output_ch.mix(
            GENERATE_CONSENSUS.out.consensus_fastq_ch
                .combine(["fastq_consensus"]),

            GENERATE_CONSENSUS.out.splitted_bam                
                .combine(["splitted_bam"]),

            GENERATE_CONSENSUS.out.splitted_bai
                .combine(["splitted_bam"])

            )
            .set{output_ch2}

    } else {
        output_ch2 = output_ch
    }

    publish_artifact(output_ch2)
}
