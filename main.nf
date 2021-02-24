#!/usr/bin/env nextflow

// enable dsl2
nextflow.preview.dsl = 2

params.dbs = ".github/data/resources"

// include modules
include {printHelp} from './modules/help.nf'
include {makeFastqSearchPath} from './modules/util.nf'

// import subworkflows
include {vp_QualityControl} from './workflows/Illumina.nf'
include {vp_Dehosting} from './workflows/Illumina.nf'




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
}
else {
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
       vp_QualityControl(ch_filePairs)

       Channel.fromPath( "$params.dbs/*.fna", checkIfExists: true )
                          .set{ ch_HumanReference }


        Channel.fromPath( "$params.dbs/*.fasta", checkIfExists: true )
                          .set{ ch_vpReference }
        //vp_prepareReferenceFiles()
        vp_Dehosting(ch_filePairs, ch_vpReference, ch_HumanReference, )

     } else {
         println("Please select a workflow with --illumina or --nanopore")
     }

}
