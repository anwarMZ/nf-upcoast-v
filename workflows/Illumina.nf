#!/usr/bin/env nextflow

// enable dsl2
nextflow.preview.dsl = 2

// import modules
//include {articDownloadScheme } from '../modules/artic.nf'
include {fastQC} from '../modules/reads.nf'
include {multiQC} from '../modules/reads.nf'
include {generateCompositeReference} from '../modules/reads.nf'
include {grabCompositeIndex} from '../modules/reads.nf'
include {indexReference} from '../modules/reads.nf'
include {mapToCompositeIndex} from '../modules/reads.nf'
include {dehostBamFiles} from '../modules/reads.nf'
include {generateDehostedReads} from '../modules/reads.nf'
include {combineDehostedCSVs} from '../modules/reads.nf'

/*include {readTrimming} from '../modules/reads.nf'
include {indexReference} from '../modules/reads.nf'
include {readMapping} from '../modules/reads.nf'
*/



// import subworkflows
//include {Genotyping} from './typing.nf'

workflow vp_QualityControl {
    take:
      ch_filePairs

    main:
      fastQC(ch_filePairs)
      multiQC(fastQC.out.collect())


      // Actually do analysis
      //sequenceAnalysis(dehosting.out.dehosted, prepareReferenceFiles.out.bwaindex, prepareReferenceFiles.out.bedfile)


}

workflow vp_Dehosting {
    take:
      ch_filePairs
      ch_vpReference
      ch_HumanReference

    main:

      generateCompositeReference(ch_HumanReference, ch_vpReference)

      if ( params.bwa_index ){
        grabCompositeIndex("${params.bwa_index}")
        grabCompositeIndex.out
                .set{ ch_index }
      } else {
        indexReference(generateCompositeReference.out)
        indexReference.out
                .set{ ch_index }
      }

      mapToCompositeIndex(ch_filePairs.combine(generateCompositeReference.out.fasta),ch_index)

      dehostBamFiles(mapToCompositeIndex.out.bam)

      generateDehostedReads(dehostBamFiles.out.bam)

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
