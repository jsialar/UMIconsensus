include { PARSE_BED; GROUP_UMI } from '../modules/GROUP_UMI.nf'
include { SPLIT_GROUP } from '../modules/SPLIT_GROUP.nf'
include { CALL_CONSENSUS } from '../modules/CALL_CONSENSUS.nf'
include { PROCESS_CONSENSUS } from '../modules/PROCESS_CONSENSUS.nf'

workflow GENERATE_CONSENSUS{
    take:
    bamfile

    main:
    // Use Nextflow to split input bed file
    Channel.fromPath( file("${params.bed}", checkIfExists: true) )
    .splitText()
    .map { line ->
        def fields = line.tokenize('\t')
        def target = fields[3].trim()  // Trim to remove trailing \n or spaces
        return tuple(target, line)
    }
    .set{ bed_channel }
    
    PARSE_BED( bed_channel )

    bamfile
    .combine(PARSE_BED.out.bed_channel)
    .set{BAM_ch}

    GROUP_UMI(BAM_ch)
    
    GROUP_UMI.out.umi_sizes
    .collectFile{ file -> ["${file[0]}_umisizes.tsv", file[1]] }
    .set{umi_sizes_ch}

    GROUP_UMI.out.grouped_bam
    .filter { tuple -> 
        def bam_file = tuple[2] 
        bam_file.size() > 0
    }
    .set{filtered_grouped_bam}

    SPLIT_GROUP(filtered_grouped_bam)

    CALL_CONSENSUS(SPLIT_GROUP.out.splitted_bam)

    CALL_CONSENSUS.out.consensus_table
    .collectFile{ file -> ["${file[0]}_consensustable.csv", file[1]] }
    .set{consensustable_ch}

    PROCESS_CONSENSUS(CALL_CONSENSUS.out.fastq_w_table)

    PROCESS_CONSENSUS.out.qual
    .collectFile{ file -> ["${file[0]}_quals.csv", file[1]] }
    .set{quals_ch}

    PROCESS_CONSENSUS.out.hap
    .collectFile{ file -> ["${file[0]}_haps.csv", file[1]] }
    .set{haps_ch}

    emit:
    umi_sizes_ch
    consensustable_ch
    consensus_fastq = CALL_CONSENSUS.out.consensus_fastq
    consensus_quals = quals_ch
    consensus_haps = haps_ch

}