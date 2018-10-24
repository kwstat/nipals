# food.R
# Time-stamp: <23 Oct 2018 13:16:54 c:/x/rpack/nipals/old/food.R>

# Source
# John Hartigan, "Clustering Algorithms", page 289.
# The percentage of all households with various foods.
# Data has 3/320 missing values. Grung has a worked example, but
# only shows a plot, no values.

# Electronic version from (via GPL)
# https://people.sc.fsu.edu/~jburkardt/datasets/hartigan/file45.txt
# Note, 00 values are really NA.

# This data can also be found in cluster.datasets::european.foods
# but again the 0 values should be NA

# Also used by
# Bjørn Grung, Rolf Manne (1998).
# Missing values in principal component analysis
# Chemometrics and Intelligent Laboratory Systems, 42, 125-139.

libs(asreml,dplyr,fs,kw,lattice,readxl,readr,reshape2,tibble)

setwd("c:/x/rpack/nipals/old/")
dat0 <- read_excel("food.xlsx")
dat0 <- read_csv("food.csv")
dat <- dat0
head(dat)

food0 <- tribble(
 ~ code,  ~ name,      ~ Ger, ~ Ita, ~ Fra, ~ Net, ~ Bel, ~ Lux, ~ Eng, ~ Pol, ~ Aus, ~ Swi, ~ Swe, ~ Den, ~ Nor, ~ Fin, ~ Spa, ~ Ire,
"GC", "ground coffee",     90, 82, 88, 96, 94, 97, 27, 72, 55, 73, 97, 96, 92, 98, 70, 13,
"IC", "instant coffee",    49, 10, 42, 62, 38, 61, 86, 26, 31, 72, 13, 17, 17, 12, 40, 52,
"TB", "tea bags",          88, 60, 63, 98, 48, 86, 99, 77, 61, 85, 93, 92, 83, 84, 40, 99,
"SS", "sugarless sweets",  19, 02, 04, 32, 11, 28, 22, 02, 15, 25, 31, 35, 13, 20, NA, 11,
"BP", "packaged biscuits", 57, 55, 76, 62, 74, 79, 91, 22, 29, 31, NA, 66, 62, 64, 62, 80,
"SP", "packaged soup",     51, 41, 53, 67, 37, 73, 55, 34, 33, 69, 43, 32, 51, 27, 43, 75,
"ST", "tinned soup",       19, 03, 11, 43, 25, 12, 76, 01, 01, 10, 43, 17, 04, 10, 02, 18,
"IP", "instant potatoes",  21, 02, 23, 07, 09, 07, 17, 05, 05, 17, 39, 11, 17, 08, 14, 02,
"FF", "frozen fish",       27, 04, 11, 14, 13, 26, 20, 20, 15, 19, 54, 51, 30, 18, 23, 05,
"VF", "frozen vegetables", 21, 02, 05, 14, 12, 23, 24, 03, 11, 15, 45, 42, 15, 12, 07, 03,
"AF", "fresh apples",      81, 67, 87, 83, 76, 85, 76, 22, 49, 79, 56, 81, 61, 50, 59, 57,
"OF", "fresh oranges",     75, 71, 84, 89, 76, 94, 68, 51, 42, 70, 78, 72, 72, 57, 77, 52,
"FT", "tinned fruit",      44, 09, 40, 61, 42, 83, 89, 08, 14, 46, 53, 50, 34, 22, 30, 46,
"JS", "shop jam",          71, 46, 45, 81, 57, 20, 91, 16, 41, 61, 75, 64, 51, 37, 38, 89,
"CG", "garlic clove",      22, 80, 88, 15, 29, 91, 11, 89, 51, 64, 09, 11, 11, 15, 86, 05,
"BR", "butter",            91, 66, 94, 31, 84, 94, 95, 65, 51, 82, 68, 92, 63, 96, 44, 97,
"ME", "margarine",         85, 24, 47, 97, 80, 94, 94, 78, 72, 48, 32, 91, 94, 94, 51, 25,
"OO", "olive oil",         74, 94, 36, 13, 83, 84, 57, 92, 28, 61, 48, 30, 28, 17, 91, 31,
"YT", "yogurt",            30, 05, 57, 53, 20, 31, 11, 06, 13, 48, 02, 11, 02, NA, 16, 03,
"CD", "crispbread",        26, 18, 03, 15, 05, 24, 28, 09, 11, 30, 93, 34, 62, 64, 13, 09,
)
food <- as.matrix(as.data.frame(food0[,3:18]))
rownames(food) <- food0[[1]]

# Grung appears to not have scaled the data.
# This almost matches Grung's fig 2, but Ireland does not match.
# I checked the Ireland data and it is ok.
mod2 <- nipals(t(food), scale=FALSE)
xx=mod2$scores[,1] * mod2$eig[1] # just stretch axis to match Grung Fig 2
yy=mod2$scores[,2] * mod2$eig[2]
plot(xx,-yy)
text(xx,-yy,names(xx))

# Grung and Manne Table 1 
# 11.8, 7.63, 5.27, 3.92, 2.90
# Not sure how they calcualte Explained variance in Table 1.
# I can match the first value
round(sqrt(mod2$eig)/sum(sqrt(mod2$eig)),3)
# [1] 0.118 0.104 0.095 0.087 0.080 0.074 0.072 0.063 0.058 0.054 0.050 0.043
#[13] 0.039 0.027 0.021 0.016
