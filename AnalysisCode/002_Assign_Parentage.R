# Parentage assignment with sequoia

library(sequoia)

# Create life history input
lh <- read.csv("SharedData/LifeHistoryData.csv")

LH_dolphins <- lh[, c("Dolphin.ID", "Birth.Date", "Sex")]

LH_dolphins$Sex[LH_dolphins$Sex == "FEMALE"] <- 1
LH_dolphins$Sex[LH_dolphins$Sex == "MALE"] <- 2
LH_dolphins$Sex[LH_dolphins$Sex == "UNKNOWN"] <- 3

LH_dolphins$birthyear <- substr(LH_dolphins$Birth.Date, 1, 4)

LH_dolphins$birthyear[LH_dolphins$birthyear == ""] <- -9

LH_dolphins <- LH_dolphins[, c("Dolphin.ID", "Sex", "birthyear")]

# Load in genotypes
goodloci <- GenoConvert(InFile = "SharedData/genotypes.txt", InFormat = "seq")

set.seed(1)

ped_with_sibs <- sequoia(GenoM = goodloci,
                         LifeHistData = LH_dolphins,
                         Err = 0.05,
                         Module = "par",
                         Tassign = 0.5, # default 0.5
                         Tfilter = -2) # default -2

# save(ped_with_sibs, file = "IntermediateData/sequoia_ped_with_sibs.RData")

est_conf <- EstConf(ped_with_sibs$PedigreePar, LH_dolphins, nSim = 100, nCores = 3,
                    args.sim = list(nSnp = 400, SnpError = 0.05, ParMis = c(0.4, 0.4)),
                    args.seq = list(Module = "ped", Err = 0.05, Tassign = 0.5, CalcLLR = FALSE))

# save(est_conf, file = "IntermediateData/pedigree_est_conf.RData")
