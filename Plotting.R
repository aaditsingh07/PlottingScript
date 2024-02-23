library(biomaRt)
library(biomartr)

# Set working directory
setwd("./")

# Read CSV file
readfile <- read.csv("SRR765980.annovar.hg38_multianno1.csv")

# Extract chromosome-wise mutation counts
sums <- c(44690, 25012, 27034, 16598, 19892, 22812, 46090, 15646, 17882, 20546, 24918, 24798, 7568, 13797, 16212, 19972, 26358, 6666, 24338, 12984, 5078, 10346, 17136, 1320)
muts <- numeric()

for (i in 1:22) {
  chrnum <- paste("chr", i, sep = "")
  muts[i] <- length(which(readfile$Chr == chrnum))
}

# Special handling for chrX and chrY
muts[23] <- length(which(readfile$Chr == "chrX"))
muts[24] <- length(which(readfile$Chr == "chrY"))

# Calculate logarithms
logmuts <- log(muts, base = 10)
logsums <- log(sums, base = 10)

# Plot number of mutations by number of exons+introns+genes
plot(logsums, logmuts, main = "Number of Mutations by Number of Exons, Introns, & Genes",
     xlab = "Number of Exons, Introns, and Genes per Chromosome",
     ylab = "Number of Mutations")
abline(lm(logmuts ~ logsums), col = "red")

# Plot number of mutations per chromosome
chr_labs <- NULL
for (i in 1:22){
chr_labs <- append(chr_labs,paste("chr", i, sep = ""))
}
chr_labs <- append(chr_labs,c("chrX", "chrY"))
barplot(muts, names.arg = chr_labs, las=2, xlab = "Chromosome Number", ylab = "Number of Mutations", main="Number of Mutations per Chromosome")