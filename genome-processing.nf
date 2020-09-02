#!/usr/local/bin nextflow


VpResolver_version = '1.0.0'
nf_required_version = '20.07.1'

log.info "".center(60, "=")
log.info "VpResolver".center(60)
log.info "Version ${VpResolver_version}".center(60)
log.info "Unified Pathogen Control Onehealth Approach Specifically Targeting \
Vibrio (UPCOAST-V) - VpResolver".center(60)
log.info "".center(60, "=")


params.help = false
def printHelp() {
    log.info """
    Example usage:
    nextflow run VpResolver --reads '*_R{1,2}.fastq.gz'
    Mandatory arguments:
      --reads               Path to input data (must be surrounded with single quotes).
    Output options:
      --output_dir           Output directory, where results will be saved
                             (default: ${params.output_dir}).
    Refer to the online manual for more information on available options:
    """
}

def printSettings() {
    log.info "Running with the following settings:".center(60)
    for (option in params) {
        if (option.key in ['cluster-options', 'help']) {
            continue
        }
        log.info "${option.key}: ".padLeft(30) + "${option.value}"
    }
    log.info "".center(60, "=")
}


try {
    if ( ! nextflow.version.matches(">= $nf_required_version") ){
        throw GroovyException('Nextflow version too old')
    }
} catch (all) {
    log.error "\n" +
              "".center(60, "=") + "\n" +
              "VpResolver requires Nextflow version $nf_required_version!".center(60) + "\n" +
              "You are running version $workflow.nextflow.version.".center(60) + "\n" +
              "Please run `nextflow self-update` to update Nextflow.".center(60) + "\n" +
              "".center(60, "=") + "\n"
    exit(1)
}

if ( params.help ) {
    printHelp()
    exit(0)
}
printSettings()


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
  shovill --outdir shovill_${sample_id} --R1 ${reads_file[0]} --R2 ${reads_file[1]} --depth 0 --gsize 5.1 --cpus ${task.cpus} --ram 12 --trim
  """
}
