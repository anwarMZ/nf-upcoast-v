process fastQC{

    tag { "${params.prefix}/${sampleName}" }
    publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":","_")}", pattern: "*fastqc*", mode: "copy"

    cpus 2

    input:
    tuple(sampleName, path(forward), path(reverse))

    output:
    path("*.{zip,html}"), emit: fastqc_files

    """
    fastqc -q \
    -t ${task.cpus} ${forward} ${reverse}
    """
}

process multiQC{

    tag { "${params.prefix}" }
    publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":","_")}", pattern: "*.html", mode: "copy", overwrite: false

    //cpu 2

    input:
    path("*")

    output:
    file "${params.prefix}_multiqc_report.html"

    """
    multiqc . -s -m fastqc -n ${params.prefix}_multiqc_report.html
    """
}

process grabIndex {

    //label 'smallcpu'

    input:
    path(index_folder)

    output:
    file("*.fna.*")

    script:
    """
    ln -sf $index_folder/*.fna.* ./
    """
}

process indexReference {

    tag { "GRC38" }

    input:
    path(human_ref)

    output:
    file("*.fa*")

    script:
    """
    bwa index -a bwtsw $human_ref
    """
}

process readmapping {
    tag { "${params.prefix}/${sampleName}" }
    publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.*", mode: "copy"

    cpus 2

    input:
    tuple(sampleName, path(forward), path(reverse), path(reference))
    path(indexed_reference)

    output:
    tuple sampleName, path("${sampleName}.sorted.bam"), emit: bam
    path("${sampleName}.flagstats.txt")

    script:
    """
    bwa mem -t ${task.cpus} ${reference} ${forward} ${reverse} | samtools sort --threads 6 -T "temp" -O BAM -o ${sampleName}.sorted.bam
    samtools flagstat ${sampleName}.sorted.bam > ${sampleName}.flagstats.txt
    """
}


process dehostBamFiles {
    tag { "${params.prefix}/${sampleName}" }
    publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.dehosted.bam", mode: "copy"
    publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}*.csv", mode: "copy"
    //label 'mediumcpu'

    input:
    tuple(sampleName, path(bam))

    output:
    tuple sampleName, path("${sampleName}.dehosted.bam"), emit: dehosted_bam
    path("${sampleName}*.csv"), emit: csv

    script:

    def rev = workflow.commitId ?: workflow.revision ?: workflow.scriptId

    """
    samtools index ${bam}
    dehost.py --file ${bam} \
    -q ${params.keep_min_map_quality} \
    -Q ${params.remove_min_map_quality} \
    -o ${sampleName}.dehosted.bam \
    -R ${rev}
    """
}

process generateDehostedReads {

    tag { "${params.prefix}/${sampleName}" }
    publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}-dehosted_R*", mode: "copy"

    //label 'mediumcpu'

    input:
    tuple(sampleName, path(dehosted_bam))

    output:
    //tuple(sampleName, path("${sampleName}-dehosted_R1"), path("${sampleName}-dehosted_R2*")), emit: dehosted
    tuple sampleName, path("${sampleName}-dehosted_R1*"), path("${sampleName}-dehosted_R2*") , emit: dehosted

    script:
    """
    samtools fastq -1 ${sampleName}-dehosted_R1.fastq -2 ${sampleName}-dehosted_R2.fastq ${dehosted_bam}
    gzip ${sampleName}-dehosted_R*.fastq
    """
}

process combineDehostedCSVs {

    tag { "${params.prefix}/${sampleName}" }
    publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":","_")}", pattern: "dehosting_summary.csv", mode: "copy"
    //label 'smallcpu'

    input:
    path("*.csv")

    output:
    path("removal_summary.csv")

    script:
    """
    csvtk concat *_stats.csv > removal_summary.csv
    """
}

process readTrimming {
    /**
    * Trims paired fastq using trim_galore (https://github.com/FelixKrueger/TrimGalore)
    * @input tuple(sampleName, path(forward), path(reverse))
    * @output trimgalore_out tuple(sampleName, path("*_val_1.fq.gz"), path("*_val_2.fq.gz"))
    */

    //tag { sampleName }

    publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":","_")}", pattern: '*_val_{1,2}.fq.gz', mode: 'copy'

    cpus 8

    input:
    tuple(sampleName, path(forward), path(reverse))

    output:
    tuple(sampleName, path("*_val_1.fq.gz"), path("*_val_2.fq.gz")) optional true

    script:
    """
    if [[ \$(gunzip -c ${forward} | head -n4 | wc -l) -eq 0 ]]; then
      exit 0
    else
      trim_galore --paired $forward $reverse
    fi
    """
}

process assembly {
    publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.*", mode: "copy"

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

process assemblyQC {
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
