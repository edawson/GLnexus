task genotypeTask{
    Array[File] inputVCFs
    Array[File] inputTBIs
    File bedFile
    Int threads
    Int memory
    String outbase

    Int memGB = memory + 11

    command {
        glnexus_cli --bed ${bedFile} --threads ${threads} --mem-gbytes ${memory} --list ${write_lines(inputVCFs)} > ${outbase}.glnexus.g.bcf
    }

    runtime{
        docker : "erictdawson/glnexus"
        memory : "${memGB}" + " GB"
        cpu : "${threads}"
        preemptible_tries : 1
    }

    output{
        File glnexusBCF = "${outbase}.glnexus.g.bcf"
    }
}

workflow GLNexusWorkflow{
    Array[File] inputVCFs
    Array[File] inputTBIs
    File bedFile
    Int threads
    Int memory

    call genotypeTask{
        input:
            inputVCFs=inputVCFs,
            inputTBIs=inputTBIs,
            bedFile=bedFile,
            threads=threads,
            memory=memory
        }

}
