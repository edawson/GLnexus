task genotypeTask{
    Array[File] inputVCFs
    Array[File] inputTBIs
    File bedFile
    Int threads
    Int memory
    Int diskGB
    String outbase

    Int memGB = memory + 11

    command {
        glnexus_cli --bed ${bedFile} --threads ${threads} --mem-gbytes ${memory} --list ${write_lines(inputVCFs)} > ${outbase}.glnexus.g.bcf
    }

    runtime{
        docker : "erictdawson/glnexus"
        memory : "${memGB}" + " GB"
        disks : "local-disk " + diskGB + " HDD"
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

    Int diskGB = ceil(size(inputVCFs, "GB") + size(inputTBIs, "GB")) + 20

    call genotypeTask{
        input:
            inputVCFs=inputVCFs,
            inputTBIs=inputTBIs,
            bedFile=bedFile,
            threads=threads,
            memory=memory,
            diskGB=diskGB
        }

}
