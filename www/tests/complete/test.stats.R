outfile <- "test.stats.txt"

write("Test statistics", file = outfile, append = FALSE)
write("", file = outfile, append = TRUE)

write("=====================", file = outfile, append = TRUE)

write("Wx~Wy free, W free", file = outfile, append = TRUE)
source("tests/depmod.R")
source("write.output.R")

write("=====================", file = outfile, append = TRUE)

write("pCCA", file = outfile, append = TRUE)
source("tests/pcca.R")
source("write.output.R")

write("=====================", file = outfile, append = TRUE)

write("Wx = Wy, W free", file = outfile, append = TRUE)
source("tests/depmod.simcca.R")
source("write.output.R")

write("=====================", file = outfile, append = TRUE)

write("Wx~Wy free, W >=0", file = outfile, append = TRUE)
source("tests/pcca.nnW.R")
source("write.output.R")

write("=====================", file = outfile, append = TRUE)

write("Wx = Wy, W >=0", file = outfile, append = TRUE)
source("tests/psimcca.nnW.R")
source("write.output.R")

write("=====================", file = outfile, append = TRUE)

write("pFA", file = outfile, append = TRUE)
source("tests/pfa.R")
source("write.output.R")

write("=====================", file = outfile, append = TRUE)

write("pFA; non-negative W", file = outfile, append = TRUE)
source("tests/pfa.nnW.R")
source("write.output.R")

write("=====================", file = outfile, append = TRUE)

write("pPCA", file = outfile, append = TRUE)
source("tests/ppca.R")
source("write.output.R")

write("\n----------------------\n", file = outfile, append = TRUE)  

session.info <- sessionInfo()
save(session.info, file = "sessioninfo.RData")        

