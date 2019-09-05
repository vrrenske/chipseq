##CHIP-SEQ pipeline based on Jasper Stinenbosch' R pipeline.
##Renske van Raaphorst, august 2019
##Unil, Lausanne

###########################################################

#Define variables

#the sequence output: gunzipped fastq files
IDS, = glob_wildcards("Test_folder/{id}.fastq.gz")

###########################################################

#Set rules
rule all:
    input:
        expand("Test_folder/{id}.fastq", id = IDS),
        "rename.txt",
        "map-log.txt",
        "macs2-commands.txt",
        "Output/done.txt",
        "Output/depth_files",
        "Output/plots"


#unzip the files
rule gunzip:
    input:
        "Test_folder/{id}.fastq.gz"
    output:
        "Test_folder/{id}.fastq"
    shell:
        "gunzip {input}"

#take each pair (R1/R2) and each input vs IP.
#give simplified new name to each.
#give back grouped lists of names.

rule relocate_trim:
    input:
        "Test_folder"
    output:
        "rename.txt",
        "list-of-files.txt"
    shell:
        "python python_orderfiles_snake.py"

#map reads using bowtie2 & compress using samtools
rule map_compress:
    input:
        filelist = rules.relocate_trim.output[1],
        reference = "D39V.FASTA"
    output:
        "map-log.txt"
    run:
        shell("bowtie2-build {input.reference} reference.index")
        test2 = open(input.filelist, "r")
        test2 = test2.readlines()
        for n in range(0, len(test2)):
            test2[n] = test2[n].partition("\n")[0]
        i = open("map-log.txt", "w+")
        samlist = []
        for n in range(0, len(test2), 2):
            bn = test2[n].partition("-R")[0]
            i.write(bn + ".bam\n")
            samlist.append(bn+".sam")
            shell("bowtie2 -x reference.index -1 " + test2[n] + " -2 " + test2[n+1] + " -S " + bn + ".sam -p 8")
        i.close()
        for n in samlist:
            print(n)
            shell("samtools view -Sb " + n + " > " + n.partition(".")[0] + ".bam")
            os.remove(n)

#macs2 peak calling: make the macs2 commands using python 3
rule peaks_makecommands:
    input:
        rules.map_compress.output
    output:
        "macs2-commands.txt"
    run:
        i = open("macs2-commands.txt", "w+")
        i.write("#!/bin/bash\n")
        import glob
        for n in glob.glob("Output/*"):
            i.write("macs2 callpeak -t " + " ".join(glob.glob(n + "/IP*.bam")) + " -c " +  " ".join(glob.glob(n + "/input*.bam")) + " -n " + n + " -f BAMPE -g 2.03e+6\n")
        i.close()


#macs2 peak calling 2: run the actual calling
#run with --use-conda
rule peaks_run:
    input:
        rules.peaks_makecommands.output
    output:
        "Output/done.txt"
    conda:
        "envs/macs2.yaml"
    shell:
        """
        ./macs2-commands.txt
        echo "Macs2-done" > {output}
        """

#R: bin the reads & save as .Rda
rule sam_depth:
    input:
        rules.map_compress.output
    output:
        directory("Output/depth_files")
    run:
        import os
        import glob
        os.mkdir("Output/depth_files")
        for n in glob.glob("Output/*/*.bam"):
            shell("samtools sort " + n + " -o " +  n)
            shell("samtools depth -a " + n + " > Output/depth_files/" + n.split("/")[2].partition(".bam")[0]+".txt")


#R: use plotly to plot chip over sequence
rule R_plots:
    input:
        rules.sam_depth.output
    output:
        directory("Output/plots")
    shell:
        """
        mkdir {output}
        Rscript convertBAMs.R
        """
