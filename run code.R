setwd("Y:/GROUPS/gr_Veening/Users/Renske/chip-seq/analysis/radA/6-BamFiles")
source("convertBams.R")


input_39121 <- "input-3912-1.bam"
IP_39121 <- "IP-3912-1.bam"

rel_39121 <- getBAMS(IP_39121, input_39121)
save(rel_39121, file = "read_3912_1.Rda")
gc() #clear memory

input_39122 <- "input-3912-2.bam"
IP_39122 <- "IP-3912-2.bam"

rel_39122 <- getBAMS(IP_39122, input_39122)
save(rel_39122, file = "read_3912_2.Rda")

gc() #clear memory

input_3912CSP_1 <- "input-CSP3912-1.bam"
IP_3912CSP_1 <- "IP-CSP3912-1.bam"

rel_3912CSP_1 <- getBAMS(IP_3912CSP_1,input_3912CSP_1)

save(rel_3912CSP_1, file="read_3912CSP_1.Rda")

gc()

input_3912CSP_2 <- "input-CSP3912-2.bam"
IP_3912CSP_2 <- "IP-CSP3912-2.bam"

rel_3912CSP_2 <- getBAMS(IP_3912CSP_2, input_3912CSP_2)
save(rel_3912CSP_2, file="read_3912CSP_2.Rda")
gc()

input_36041 <- "input-n3604-1.bam"
IP_36041 <- "IP-n3604-1.bam"

rel_36041 <- getBAMS(IP_36041, input_36041)
save(rel_36041, file="read_3604_1.Rda")

gc()

input_36042 <- "input-n3604-2.bam"
IP_36042 <- "IP-n3604-2.bam"

rel_36042 <- getBAMS(IP_36042, input_36042)
save(rel_36042, file="read_3604_2.Rda")
gc()

input_3604CSP_1 <- "input-3604MCSP-1.bam"
IP_3604CSP_1 <- "IP-3604MCSP-1.bam"

rel_3604CSP_1 <- getBAMS(IP_3604CSP_1, input_3604CSP_1)

save(rel_3604CSP_1, file="read_3604CSP_1.Rda")

gc()

input_3604CSP_2 <- "input-3604MCSP-2.bam"
IP_3604CSP_2 <- "IP-3604MCSP-2.bam"

rel_3604CSP_2 <- getBAMS(IP_3604CSP_2, input_3604CSP_2)
save(rel_3604CSP_2, file="read_3604CSP_2.Rda")
gc()

input_4000CSP_2 <- "input-4000CSP-2.bam"
IP_4000CSP_2 <- "IP-4000CSP-2.bam"

rel_4000CSP_2 <- getBAMS(IP_4000CSP_2, input_4000CSP_2)
save(rel_4000CSP_2, file="read_4000CSP_2.Rda")
gc()