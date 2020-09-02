#!/usr/local/bin nextflow


// Using checkIfExists & followLinks to confirm if the files are avaialable and
// also for IRIDA configuration to follow symlinks
Channel.fromFilePairs('data/*_R{1,2}*.fastq.gz',
                                  checkIfExists: true, followLinks: true)
                                  .into {samples_fastqc_ch; samples_assembly_ch}


process fastqc {

  input:
  tuple val(sample_id), file(reads_file) from samples_fastqc_ch

  output:
  file("fastqc_${sample_id}_logs") into fastqc_ch

  script:
  """
  mkdir fastqc_${sample_id}_logs
  fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads_file}
  """

}


process assembly {

  input:
  tuple val(sample_id), file(reads_file) from samples_assembly_ch

  output:
  path('shovill_${sample_id}') into assembly_ch_result

  script:
  """
  shovill --outdir shovill_${sample_id} --R1 ${reads_file[0]} --R2 ${reads_file[1]} --depth 0 --gsize 5.1 --cpus 2 --ram 12 --trim
  """
}
