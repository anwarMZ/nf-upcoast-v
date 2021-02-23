process genome_annotation {
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.*", mode: "copy"

    label 'prokka annotation'

    input:
    tuple(sampleName, path(forward), path(reverse))

    output:
    tuple sampleName, path("${sampleName}.sorted.bam"), emit: bam
    path("${sampleName}.flagstats.txt")

    script:
    """
    bwa mem -t 8 ${composite_reference} ${forward} ${reverse} | samtools sort --threads 6 -T "temp" -O BAM -o ${sampleName}.sorted.bam
    samtools flagstat ${sampleName}.sorted.bam > ${sampleName}.flagstats.txt
    """
}

process pangenome_analysis {
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.*", mode: "copy"

    label 'shovill assembly'

    input:
    tuple(sampleName, path(forward), path(reverse))

    output:
    tuple sampleName, path("${sampleName}.sorted.bam"), emit: bam
    path("${sampleName}.flagstats.txt")

    script:
    """
    bwa mem -t 8 ${composite_reference} ${forward} ${reverse} | samtools sort --threads 6 -T "temp" -O BAM -o ${sampleName}.sorted.bam
    samtools flagstat ${sampleName}.sorted.bam > ${sampleName}.flagstats.txt
    """
}
