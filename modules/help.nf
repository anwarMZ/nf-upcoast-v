def printHelp() {
  log.info"""
  Usage:
    nextflow run nf-upcoast-v -profile (conda) ( --illumina ) --prefix [prefix] [workflow-options]

  Description:
    Analyze vibrio parahaemolyticus Whole Genomic Sequencing (WGS) datasets from clinical and environmental samples


    All options set via CLI can be set in conf directory

  Nextflow arguments (single DASH):
    -profile                  Allowed values: conda

  Mandatory workflow arguments (mutually exclusive):
    --illumina                Run the Illumina workflow

  Illumina workflow options:
    Mandatory:
      --prefix                A (unique) string prefix for output files.
                              Sequencing run name is a good choice e.g DDMMYY_MACHINEID_RUN_FLOWCELLID.
      --directory             Path to a directory containing paired-end Illumina reads.
                              Reads will be found and paired RECURSIVELY beneath this directory.
    Optional:
      --outdir                Output directory (Default: ./results)
      --illuminaKeepLen       Length (bp) of reads to keep after primer trimming (Default: 20)
      --illuminaQualThreshold Sliding window quality threshold for keeping
                              reads after primer trimming (Default: 20)
      --ref                   Path to Human reference fasta file
                              
  """.stripIndent()
}
