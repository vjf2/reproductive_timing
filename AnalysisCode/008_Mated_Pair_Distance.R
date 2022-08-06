# Compare mated pairs distances

library(SocGen)
library(coin)

# Take all births with known paternities

load("IntermediateData/sequoia_ped_with_sibs.RData")

ped <- ped_with_sibs$PedigreePar

with_dad <- ped[!is.na(ped$sire), ]

# Centroids of all individuals who were genotyped

cpts <- read.csv("SharedData/utm_centroid_positions.csv", row.names = 1)

centroids <- dist(cpts)

dist_mat <- as.matrix(centroids)
dist_mat <- dist_mat / 1000 # in km

# For each female, get a list of males that were alive and at least 12 years old at the conception date

with_both <- with_dad[!is.na(with_dad$dam), ]

lh <- read.csv("SharedData/LifeHistoryData.csv")

with_both$conception_date <- lh$Birth.Date[match(with_both$id, lh$Dolphin.ID)]

with_both$conception_date <- as.Date(with_both$conception_date)

with_both$conception_date <- with_both$conception_date - 372

# Make an availabilty matrix for all the biopsied males
lh$Birth.Date <- as.Date(lh$Birth.Date)

parents <- lh$Dolphin.ID[which(lh$Sex == "MALE",
                                 lh$Birth.Date <= "2004-1-1")]

fast_avail <- data.frame(dolphin_id = parents, entry = rep(NA, length(parents)), depart = rep(NA, length(parents)))

fast_avail$entry <- lh$Birth.Date[match(fast_avail$dolphin_id, lh$Dolphin.ID)]


fast_avail$entry <- fast_avail$entry + (12 * 365.25)
fast_avail$depart <- lh$Depart.Study.Date[match(fast_avail$dolphin_id,
                                                  lh$Dolphin.ID)] |> as.Date()

dates <- unique(with_both$conception_date)

schedule <- schedulize(fast_avail, id = "dolphin_id", start = "entry", end = "depart",
                         dates = dates, format = "sim")

rownames(schedule) <- as.character(dates)

distances <- list()
actual_distance <- c()
percentile <- c()

for (i in 1:length(with_both$id)) {

  calf <- with_both$id[i]
  mom <- with_both$dam[i]
  dad <- with_both$sire[i]
  conception_date <- with_both$conception_date[i]

  potential_dads <- schedule[as.character(conception_date), ]
  potential_dads <- names(potential_dads[which(potential_dads == TRUE)])

  distances[[i]] <- dist_mat[mom,
                             intersect(potential_dads, colnames(dist_mat))]

  actual_distance[i] <- dist_mat[mom, dad]

  percentile[i] <- ecdf(distances[[i]])(actual_distance)

}

with_both$actual_distance <- actual_distance

windows()
# pdf(file="TablesFigures/Fig6_mated_pair_distances.pdf", width = 7, height = 6)
hist(unlist(distances), breaks = 20,
     main = NA,
     probability = TRUE,
     xlab = "Distance between centroids (km)", ylab = "Proportion",
     border = "darkgrey",
     yaxt = "n",
     ylim = c(0, 0.085),
     xlim = c(0, 50))

axis(2, las = 1)

segments(x0 = mean(actual_distance), y0 = 0,
         x1 = mean(actual_distance), y1 = 0.085,
         col = "red", lwd = 2)

segments(x0 = quantile(actual_distance, 0.025), y0 = 0,
         x1 = quantile(actual_distance, 0.025), y1 = 0.085,
         col = "red", lty = 2, lwd = 1)

segments(x0 = quantile(actual_distance, 0.975), y0 = 0,
         x1 = quantile(actual_distance, 0.975), y1 = 0.085,
         col = "red", lty = 2, lwd = 1)

dev.off()

# Compare distances using coin
dat1 <- data.frame(distance = actual_distance, type = rep("actual_distance", length(actual_distance)))
dat2 <- data.frame(distance = mean_distances, type = rep("mean_distances", length(mean_distances)))

dat <- rbind(dat1, dat2)
dat$type <- as.factor(dat$type)

coin::wilcox_test(distance ~ type, data = dat)
