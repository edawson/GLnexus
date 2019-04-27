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
        glnexus_cli --config gatk --bed ${bedFile} --threads ${threads} --mem-gbytes ${memory} --list ${write_lines(inputVCFs)} > ${outbase}.glnexus.g.bcf
    }

    runtime{
        docker : "erictdawson/glnexus"
        memory : "${memGB}" + " GB"
        disks : "local-disk " + diskGB + " HDD"
        cpu : "${threads}"
        preemptible_tries : 3
    }

    output{
        File glnexusBCF = "${outbase}.glnexus.g.bcf"
    }
}

task getArrayChromosomeVCFTask{
    Array[File] inputVCFs
    Array[File] inputTBIs
    String chrom

    Int diskGB

    command<<<
      for i in ${sep=' ' inputVCFs}; do
        outbase=$( basename $(basename $i ".gz") ".vcf")
        echo "bcftools $i ${chrom} > $outbase.${chrom}.vcf && bgzip -c $outbase.${chrom}.vcf > $outbase.${chrom}.vcf.gz" >> jfile.txt && \
        python /usr/bin/launcher.py -i jfile.txt -c 1 -n 4
    >>>

    runtime{
        docker : "erictdawson/base"
        cpu : 4
        memory : "3 GB"
        disks : "local-disk " + diskGB + " HDD"
        preemptible : 3
    }

    output{
        Array[File] chromVCFs = glob("*.${chrom}.vcf.gz")
        Array[File] chromTBIs = glob("*.${chrom}.vcf.gz.tbi")
    }
}

workflow GLNexusWorkflow{
    Array[File] inputVCFs
    Array[File] inputTBIs
    File bedFile
    File chromFile
    Int threads
    Int memory

    ## We can't use dynamic sizing because
    ## it is not yet implemented in Firecloud for arrays
    Int diskGB = 3250

    Array[String] chroms = read_lines(chromFile)

    scatter ( chrom in chroms ){
        call getArrayChromosomeVCFTask{
            input:
                inputVCFs=inputVCFs,
                inputTBIs=inputTBIs,
                chrom=chrom,
                diskGB=diskGB
        }

        Int chromDiskGB = 500

        call genotypeTask{
            input:
                inputVCFs=getArrayChromosomeVCFTask.chromVCFs,
                inputTBIs=getArrayChromosomeVCFTask.chromTBIs,
                bedFile=bedFile,
                threads=threads,
                memory=memory,
                diskGB=chromDiskGB
        }
    }

}
