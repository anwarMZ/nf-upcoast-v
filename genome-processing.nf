#!/usr/local/bin nextflow


// Using checkIfExists & followLinks to confirm if the files are avaialable and
// also for IRIDA configuration to follow symlinks
samples_ch = Channel.fromFilePairs('data/*{1,2}.fastq.gz', checkIfExists: true,
                                followLinks: true)
//samples_ch.view()

process fastqc {

  input:
  tuple val(sample_id), file(reads_file) from samples_ch

  output:
  file("fastqc_${sample_id}_logs") into fastqc_ch

  script:
  """
  mkdir fastqc_${sample_id}_logs
  fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads_file}
  """

}
