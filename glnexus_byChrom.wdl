task genotypeTask{
    Array[File] inputVCFs
    Array[File] inputTBIs
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

task getArrayChromosomeVCFTask{
    Array[File] inputVCFs
    Array[File] inputTBIs
    String chrom

    Int? threads = 4

    Int diskGB

    command<<<
      while read i
      do
        ln -s $i ./$(basename $i) && \
        ln -s $i.tbi ./$(basename $i).tbi && \
        shortname=$(basename $i) && \
        outbase=$( basename $(basename $i ".gz") ".vcf") && \
        echo "tabix -f $shortname && tabix -h $shortname ${chrom} > $outbase.${chrom}.vcf && bgzip -c $outbase.${chrom}.vcf > $outbase.${chrom}.vcf.gz && tabix $outbase.${chrom}.vcf.gz" >> ./jfile.txt
      done < ${write_lines(inputVCFs)} &&
      python /usr/bin/launcher.py -i ./jfile.txt -c 1 -n ${threads}
    >>>

    runtime{
        docker : "erictdawson/base"
        cpu : "${threads}"
        memory : "3 GB"
        disks : "local-disk " + diskGB + " SSD"
        preemptible : 4
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
    String outbase
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
                diskGB=diskGB,
                threads=threads
        }

        Int chromDiskGB = 500

        call genotypeTask{
            input:
                inputVCFs=getArrayChromosomeVCFTask.chromVCFs,
                inputTBIs=getArrayChromosomeVCFTask.chromTBIs,
                bedFile=bedFile,
                threads=threads,
                memory=memory,
                diskGB=chromDiskGB,
                outbase=outbase,
                chrom=chrom
        }
    }

}
