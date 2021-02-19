#!/usr/bin/env nextflow

// enable dsl2
nextflow.preview.dsl = 2

// import modules
//include {articDownloadScheme } from '../modules/artic.nf'
include {fastQC} from '../modules/reads.nf'
include {multiQC} from '../modules/reads.nf'
//include {generateCompositeReference} from '../modules/reads.nf'
include {grabIndex} from '../modules/reads.nf'
include {indexReference} from '../modules/reads.nf'
include {readmapping} from '../modules/reads.nf'
include {dehostBamFiles} from '../modules/reads.nf'
include {generateDehostedReads} from '../modules/reads.nf'
include {combineDehostedCSVs} from '../modules/reads.nf'

/*include {readTrimming} from '../modules/reads.nf'
include {indexReference} from '../modules/reads.nf'
include {readMapping} from '../modules/reads.nf'
*/



// import subworkflows
//include {Genotyping} from './typing.nf'

/*
workflow prepareReferenceFiles {
    // Get reference fasta

    if (params.ref) {
      Channel.fromPath(params.ref)
              .set{ ch_refFasta }
    }


   Either get BWA aux files from reference
     location or make them fresh

    if (params.ref) {
      // Check if all BWA aux files exist, if not, make them
      bwaAuxFiles = []
      refPath = new File(params.ref).getAbsolutePath()
      new File(refPath).getParentFile().eachFileMatch( ~/.*.bwt|.*.pac|.*.ann|.*.amb|.*.sa/) { bwaAuxFiles << it }

      if ( bwaAuxFiles.size() == 5 ) {
        Channel.fromPath( bwaAuxFiles )
               .set{ ch_bwaAuxFiles }

        ch_refFasta.combine(ch_bwaAuxFiles.collect().toList())
                   .set{ ch_preparedRef }
      } else {
        indexReference(ch_refFasta)
        indexReference.out
                      .set{ ch_preparedRef }
      }
    }

    else {
      indexReference(ch_refFasta)
      indexReference.out
                    .set{ ch_preparedRef }
    }


    emit:
      bwaindex = ch_preparedRef
      //bedfile = ch_bedFile
      reffasta = ch_refFasta
}
*/

workflow vp_Dehosting {
    take:
      ch_filePairs
      ch_HumanReference

    main:

      //generateReference(ch_HumanReference)

      if ( params.bwa_index ){
        grabIndex("${params.bwa_index}")

        grabIndex.out
                .set{ ch_index }
      } else {
        indexReference(ch_HumanReference)
        indexReference.out
                .set{ ch_index }
      }

      readmapping(ch_filePairs.combine(ch_HumanReference),ch_index)

      dehostBamFiles(readmapping.out.bam)

      generateDehostedReads(dehostBamFiles.out.dehosted_bam)

      generateDehostedReads.out.dehosted.set{ ch_filePairs_dehosted }
      combineDehostedCSVs(dehostBamFiles.out.csv.collect())

    emit:
      dehosted = ch_filePairs_dehosted

}

/*workflow sequenceAnalysis {
    take:
      ch_filePairs_dehosted
      ch_preparedRef
      ch_bedFile


    main:
      readTrimming(ch_filePairs_dehosted)

      readMapping(readTrimming.out.combine(ch_preparedRef))

      trimPrimerSequences(readMapping.out.combine(ch_bedFile))

      callVariants(trimPrimerSequences.out.ptrim.combine(ch_preparedRef.map{ it[0] }))

      makeConsensus(trimPrimerSequences.out.ptrim)

      makeQCCSV(trimPrimerSequences.out.ptrim.join(makeConsensus.out, by: 0)
                                   .combine(ch_preparedRef.map{ it[0] }))

      makeQCCSV.out.csv.splitCsv()
                       .unique()
                       .branch {
                           header: it[-1] == 'qc_pass'
                           fail: it[-1] == 'FALSE'
                           pass: it[-1] == 'TRUE'
                       }
                       .set { qc }

      writeQCSummaryCSV(qc.header.concat(qc.pass).concat(qc.fail).toList())

      collateSamples(qc.pass.map{ it[0] }
                           .join(makeConsensus.out, by: 0)
                           .join(trimPrimerSequences.out.mapped))

      if (params.outCram) {
        bamToCram(trimPrimerSequences.out.mapped.map{it[0] }
                        .join (trimPrimerSequences.out.ptrim.combine(ch_preparedRef.map{ it[0] })) )

      }

    emit:
      qc_pass = collateSamples.out
      variants = callVariants.out.variants
}*/

workflow vp_QualityControl {
    take:
      ch_filePairs

    main:
      fastQC(ch_filePairs)
      multiQC(fastQC.out.collect())
      // Build or download fasta, index and bedfile as required
      //prepareReferenceFiles()

      // Dehost illumina sequence files for analysis
      //dehost(ch_filePairs, humanRef_ch)
      //Channel.fromPath( "$params.ref/*.fna", checkIfExists: true )
      //                    .set{ ch_HumanReference }

      dehosting(ch_filePairs, ch_HumanReference)



      // Actually do analysis
      //sequenceAnalysis(dehosting.out.dehosted, prepareReferenceFiles.out.bwaindex, prepareReferenceFiles.out.bedfile)


}
