task genotypeTask{
    Array[File] inputVCFs
    File bedFile
    Int threads
    Int memory
    Int diskGB
    String outbase
    String chrom

    Int memGB = memory + 11

    command {
        export LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libjemalloc.so && \
        glnexus_cli --config gatk \
        --bed ${bedFile} --threads ${threads} \
        --mem-gbytes ${memory} \
        --list ${write_lines(inputVCFs)} > ${outbase}.glnexus.${chrom}.g.bcf
    }

    runtime{
        docker : "erictdawson/glnexus"
        memory : "${memGB}" + " GB"
        disks : "local-disk " + diskGB + " SSD"
        cpu : "${threads}"
        preemptible : 4
    }

    output{
        File glnexusBCF = "${outbase}.glnexus.${chrom}.g.bcf"
    }
}

task sliceVCFByChromTask{
    File inputVCF
    String chrom
    Int diskGB
    Int threads = 1

    String outbase = basename(basename(inputVCF, ".gz"), ".vcf")

    command{
        tabix -f ${inputVCF} && tabix -h ${inputVCF} ${chrom} | bgzip -c -@ ${threads} > ${outbase}.${chrom}.vcf.gz
    }

    runtime{
        docker : "erictdawson/samtools"
        cpu : "${threads}"
        memory : "1GB"
        disks : "local-disk " + diskGB + " HDD"
        preemptible : 4
    }

    output{
        File chromVCFGZ = "${outbase}.${chrom}.vcf.gz"
    }
}

task filterSingleChromVCFTask{
    String chrom
    String vcf

    command{
        if [ $(echo ${vcf} | grep "${chrom}.vcf") ]; then touch yes; fi
    }

    runtime{
        docker : "erictdawson/base"
        cpu : 1
        memory : "1 GB"
    }

    output{
        Array[File] isCorrectChrom = glob("yes")
    }
}

task filterChromVCFsTask{
    String chrom
    Array[String] vcfs
    Int? diskGB = 10

    command{
        grep "${chrom}.vcf" ${write_lines(vcfs)} > files.txt
    }

    runtime{
        docker : "erictdawson/base"
        cpu : 1
        memory : "1 GB"
        disks : "local-disk " + diskGB + " HDD"
        preemptible : 3
    }

    output{
        Array[String] chromFiles = read_lines("files.txt")
    }
}

task bcftoolsConcatTask{
    Array[File] chromLevelBCFs
    String outbase
    Int threads
    Int mem
    Int diskGB

    command{
        bcftools concat --threads ${threads} --output-type z --output ${outbase}.vcf.gz ${sep=" " chromLevelBCFs}
    }

    runtime{
        docker : "erictdawson/base"
        cpu : "${threads}"
        memory : mem + " GB"
        disks : "local-disk " + diskGB + " HDD"
        preemptible : 2
    }

    output{
        File concatVCFGZ = "${outbase}.vcf.gz"
    }
}

workflow GLNexusWorkflow{
    Array[File] inputVCFs
    Array[File] inputTBIs
    File bedFile
    File chromFile
    String outbase
    Int threads
    Int memory

    Int chromDiskGB = 500
    Int bcfDiskGB = 100
    
    
    Array[String] chroms = read_lines(chromFile)

    Array[Pair[String, File]] chrom_by_samp = cross(chroms, inputVCFs)

    scatter(chrom_sample in chrom_by_samp){
        Int diskGB = ceil(size(chrom_sample.right, "GB")) + 20
        call sliceVCFByChromTask{
            input:
                inputVCF=chrom_sample.right,
                chrom=chrom_sample.left,
                diskGB=diskGB,
                threads=4
        }
    }

    scatter (chrom in chroms){
        Array[String] c = [chrom]
        ## Get the chromosome-specific VCFs for every sample
        call filterChromVCFsTask{
            input:
                chrom=chrom,
                vcfs=sliceVCFByChromTask.chromVCFGZ
        }
        ## Trim to just the chrom-specific VCFs for this chromosome
        #Array[File] chrom_vcfs = select_all(filterChromVCFsTask.defined)

        ## Genotype only the selected VCFs for each chrom
 
        call genotypeTask{
            input:
                inputVCFs=filterChromVCFsTask.chromFiles,
                bedFile=bedFile,
                threads=threads,
                memory=memory,
                diskGB=chromDiskGB,
                outbase=outbase,
                chrom=chrom
            }
    }


    call bcftoolsConcatTask{
        input:
            chromLevelBCFs=genotypeTask.glnexusBCF,
            outbase=outbase,
            threads=4,
            mem=3,
            diskGB=bcfDiskGB

    }

}
