# kamyr_digester.R
# Time-stamp: <22 Oct 2018 19:24:55 c:/x/rpack/nipals/old/kamyr_digester.R>

# Source: http://openmv.net/
# The data does not appear to be used in the book at the above link.

# This data has 352/6622 = 5% missing values

libs(asreml,dplyr,fs,kw,lattice,readxl,readr,reshape2,tibble)

setwd("c:/x/rpack/nipals/old/")
dat0 <- read_csv("kamyr_digester.csv")
dat <- dat0
head(dat)

kd <- as.matrix(as.data.frame(dat[,-1]))
rownames(kd) <- dat[[1]]
dim(kd) # 301 22

mod1 <- nipals(kd)
mod1$R2
