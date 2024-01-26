
# These data are analyzed by this paper:
# The impact of missing measurements on PCA and PLS prediction and monitoring applications
# But not very usefully for us.  (Why not?)

# Pulp quality is measured by the lignin content remaining in the pulp: the Kappa number. This data set is used to understand which variables in the process influence the Kappa number, and if it can be predicted accurately enough for an inferential sensor application.

dat <- read.csv("https://openmv.net/file/kamyr-digester.csv")
head(dat)

out <- nipals(dat[ , -1])

# data also available in R
library(VIM)
data(pulplignin)
dat <- pulplignin
str(dat) # 301 obs, 23 vars

VIM::aggr(dat)

m1 <- nipals::nipals(dat[,-1])
