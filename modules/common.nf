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