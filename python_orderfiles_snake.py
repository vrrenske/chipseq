
####python recognise each pair
###RvR

import glob
import os

folder = "Test_folder/*"

ids = glob.glob(folder)

ids.sort()

trimstart = "trimmomatic PE "
trimend = " ILLUMINACLIP:TruSeq3-PE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:100"

f = open("rename.txt", "w+")

f.write("#!/bin/bash\n")

if not os.path.exists("Output"):
    os.mkdir("Output")

for n in range(0, len(ids), 4):
    bn = ids[n].partition("-")[0].partition("/")[2]
    repnum = bn[-1:]
    bn = bn[:-1]

    if not os.path.exists("Output/" + bn):
        os.mkdir("Output/" + bn)

    os.rename(ids[n], "Output/" + bn + "/IP-" + bn + "-" + repnum + "-R1.fastq" )
    os.rename(ids[n+1], "Output/" + bn + "/IP-" + bn  + "-" + repnum + "-R2.fastq")
    os.rename(ids[n+2], "Output/" + bn + "/input-" + bn  + "-" + repnum + "-R1.fastq")
    os.rename(ids[n+3], "Output/" + bn + "/input-" + bn  + "-" + repnum + "-R2.fastq")

    f.write(trimstart + " ".join(glob.glob("Output/"+bn+"/IP*")) +
    " Output/"+bn+"/IP-"+bn + "-" + repnum + "-R1-paired.fastq" +
    " Output/"+bn+"/IP-"+bn + "-" + repnum + "-R1-unpaired.fastq" +
    " Output/"+bn+"/IP-"+bn + "-" + repnum +"-R2-paired.fastq" +
    " Output/"+bn+"/IP-"+bn + "-" + repnum +"-R2-unpaired.fastq" +
    trimend + "\n")

    f.write(trimstart + " ".join(glob.glob("Output/"+bn+"/input*")) +
    " Output/"+bn+"/input-"+bn + "-" + repnum+"-R1-paired.fastq" +
    " Output/"+bn+"/input-"+bn + "-" + repnum+"-R1-unpaired.fastq" +
    " Output/"+bn+"/input-"+bn + "-" + repnum+"-R2-paired.fastq" +
    " Output/"+bn+"/input-"+bn + "-" + repnum+"-R2-unpaired.fastq" +
    trimend + "\n")

f.close()

os.system("./rename.txt")

plantorem = glob.glob("Output/*/*unpaired.fastq")
plantorem = plantorem + glob.glob("Output/*/*R1.fastq") + glob.glob("Output/*/*R2.fastq")
for n in plantorem:
    os.remove(n)

z = open("list-of-files.txt", "w+")

for n in glob.glob("Output/*/*"):
    z.write(n)
    z.write("\n")

z.close()
