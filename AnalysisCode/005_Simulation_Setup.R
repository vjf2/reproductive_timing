## Set up data for simulation
source("AnalysisCode/get_age_function.R")

library(SocGen)

# Select subjects who are sexed, available after 1984 and at least 12 by 2016

lh <- read.csv("SharedData/LifeHistoryData.csv")

lh$Birth.Date <- as.Date(lh$Birth.Date)

parents <- lh$Dolphin.ID[which(lh$Sex %in% c("MALE", "FEMALE"),
                               lh$Birth.Date <= "2004-1-1")]

# Create availability list of all individuals

fast_avail <- data.frame(dolphin_id = parents,
                         entry = rep(NA, length(parents)),
                         depart = rep(NA, length(parents)))

fast_avail$entry <- lh$Birth.Date[match(fast_avail$dolphin_id, lh$Dolphin.ID)]

fast_avail$entry <- fast_avail$entry + (10 * 365.25)
fast_avail$depart <- lh$Depart.Study.Date[match(fast_avail$dolphin_id, lh$Dolphin.ID)] |> as.Date()

# Give the model a schedule with TRUE/FALSE for each dolphin's avaiability on each survey day

dates <- seq.Date(as.Date("1984-11-01"),
                  as.Date("2017-11-01"), by = "year")

schedule <- schedulize(fast_avail, id = "dolphin_id", start = "entry", end = "depart",
                       dates = dates, format = "sim")

rownames(schedule) <- as.character(dates)

rownames(schedule) <- format(as.Date(rownames(schedule)), "%Y")

age_mat <- schedule

# Matrix of ages
for (i in 1:ncol(age_mat)) {
  dol <- colnames(age_mat)[i]
  birthyear <- format(lh$Birth.Date[which(lh$Dolphin.ID == dol)], "%Y") |> as.numeric()

  age_mat[, i] <- as.numeric(rownames(age_mat)) - birthyear

}

schedule <- schedule * age_mat

females <- lh$Dolphin.ID[which(lh$Sex == "FEMALE")]
males <- lh$Dolphin.ID[which(lh$Sex == "MALE")]

moms <- schedule[, which(colnames(schedule) %in% females)]

dads <- schedule[, which(colnames(schedule) %in% males)]

sampled <- lh$Dolphin.ID[which(lh$Genotype.Available == "Y")]

offspring <- lh[which(lh$Genotype.Available == "Y"),
                c("Dolphin.ID", "Birth.Date")]

offspring$Conception.Date <- offspring$Birth.Date - (53 * 7) # assume 53 weeks gestation
offspring$Conception.Year <- format(offspring$Conception.Date, "%Y") |> as.numeric()
offspring <- offspring[which(offspring$Conception.Year >= 1984), ]

offspring$ranmom <- ""
offspring$randad <- ""

# save(offspring, moms, dads, get_age, lh, sampled, file = "IntermediateData/simulation_setup.RData")
