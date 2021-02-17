#!/usr/bin/env nextflow

// enable dsl2
nextflow.preview.dsl = 2

//params.dbs = ".github/data/resources"

// include modules
include {printHelp} from './modules/help.nf'
include {makeFastqSearchPath} from './modules/util.nf'

// import subworkflows
include {vp_Illumina} from './workflows/illumina.nf'


if (params.help){
    printHelp()
    exit 0
}

if (params.profile){
    println("Profile should have a single dash: -profile")
    System.exit(1)
}

if ( params.illumina ) {
   if ( !params.directory ) {
       println("Please supply a directory containing fastqs with --directory.")
       println("Use --help to print help")
       System.exit(1)
   }

   //if ( (params.bed && ! params.ref) || (!params.bed && params.ref) ) {
    //   println("--bed and --ref must be supplied together")
    //   System.exit(1)
   //}
} /*else if ( params.nanopolish ) {
   if (! params.basecalled_fastq ) {
       println("Please supply a directory containing basecalled fastqs with --basecalled_fastq. This is the output directory from guppy_barcoder or guppy_basecaller - usually fastq_pass. This can optionally contain barcodeXX directories, which are auto-detected.")
   }
   if (! params.fast5_pass ) {
       println("Please supply a directory containing fast5 files with --fast5_pass (this is the fast5_pass directory)")
   }
   if (! params.sequencing_summary ) {
       println("Please supply the path to the sequencing_summary.txt file from your run with --sequencing_summary")
       System.exit(1)
   }
   if ( params.bed || params.ref ) {
       println("ivarBed and alignerRefPrefix only work in illumina mode")
       System.exit(1)
   }
} else if ( params.medaka ) {
   if (! params.basecalled_fastq ) {
       println("Please supply a directory containing basecalled fastqs with --basecalled_fastq. This is the output directory from guppy_barcoder or guppy_basecaller - usually fastq_pass. This can optionally contain barcodeXX directories, which are auto-detected.")
   }
} */else {
       println("Please select a workflow with --nanopolish, --illumina or --medaka, or use --help to print help")
       System.exit(1)
}


if ( ! params.prefix ) {
     println("Please supply a prefix for your output files with --prefix")
     println("Use --help to print help")
     System.exit(1)
} else {
     if ( params.prefix =~ /\// ){
         println("The --prefix that you supplied contains a \"/\", please replace it with another character")
         System.exit(1)
     }
}



// main workflow
workflow {
   if ( params.illumina ) {
     fastqSearchPath = makeFastqSearchPath( params.illuminaPrefixes, params.illuminaSuffixes, params.fastq_exts )
	   Channel.fromFilePairs( fastqSearchPath, flat: true)
	          .filter{ !( it[0] =~ /Undetermined/ ) }
	          .set{ ch_filePairs }
   }

   main:

     if ( params.illumina ) {
       //println("This will call Illumina workflow")
       Illumina(ch_filePairs)
     } else {
         println("Please select a workflow with --nanopolish, --illumina or --medaka")
     }

}
