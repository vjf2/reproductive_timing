## Create Figure 1

# Frequency distributions of age ranges for all adult females and males (nf = 173, nm = 162) from which genetic data were obtained (left) and total observed individuals (nf = 273, nm = 254) during the sampling period from 2013 to 2019. The 1.07:1 female:male ratio in the sampled dataset is representative of the 1.07:1 ratio observed in the population, but substantially lower than the 2.01:1 obtained parentage assignment ratio 

lh <- read.csv("SharedData/LifeHistoryData.csv")

# Note - Depart Study Date is either Death Date, six months after last sighting, or end of biopsy sampling period, whichever is earliest

aggregate(Age.Depart ~ Sex,
          data = lh[which(lh$Genotype.Available == "Y"), ],
          mean)

genotyped_adults <- lh[which(lh$Age.Depart >= 10 &
                               lh$Genotype.Available == "Y"), ]

table(genotyped_adults$Sex)
table(genotyped_adults$Sex)[1] / table(genotyped_adults$Sex)[2] # ratio

range(genotyped_adults$Age.Depart)

fem_cuum_age <- sapply(10:49,
                       function(x) {
                         sum(x <= genotyped_adults$Age.Depart[which(genotyped_adults$Sex == "FEMALE")])})

mal_cuum_age <- sapply(10:49,
                       function(x) {
                         sum(x <= genotyped_adults$Age.Depart[which(genotyped_adults$Sex == "MALE")])})

# Now get total for population
# All adults seen during sampling period

all_adults <- lh[which(lh$Age.Depart >= 10 &
                         lh$Depart.Study.Date >= "2013-07-14" &
                         lh$Sex %in% c("MALE", "FEMALE")), ]

table(all_adults$Sex)
table(all_adults$Sex)[1] / table(all_adults$Sex)[2] # ratio

range(all_adults$Age.Depart)

fem_cuum_age_total <- sapply(10:52,
                             function(x) {
                               sum(x <= all_adults$Age.Depart[which(all_adults$Sex == "FEMALE")])})

mal_cuum_age_total <- sapply(10:52,
                             function(x) {
                               sum(x <= all_adults$Age.Depart[which(all_adults$Sex == "MALE")])})

windows()
# pdf("TablesFigures/Fig1_age_distributions_R1.pdf", width = 8, height = 4)
par(mfrow = c(1, 2), mar = c(5.1, 2.5, 4.1, 0.5), oma = c(0, 2, 0, 0))

barplot(fem_cuum_age, names.arg = 10:49, col = adjustcolor("purple", alpha.f = 0.5),
        xaxt = "n",
        yaxt = "n", main = "Genetic Data", border = NA, xlab = "Age (yrs)") -> locs1
barplot(mal_cuum_age, names.arg = 10:49, col = adjustcolor("darkgreen", alpha.f = 0.5),
        xaxt = "n",
        yaxt = "n", add = TRUE, border = NA)
mtext("Count", side = 2, line = 3)

axis(1, at = c(locs1[seq(1, length(locs1), 5)], 48.7), labels = seq(10, 50, 5))
axis(2, las = 1)

text(24, 141, "excess female counts", col = adjustcolor("purple", alpha.f = 0.8),
     cex = 0.75)

arrows(x0 = 11, y0 = 134, x1 = 14, y1 = 138,
       angle = 25, length = 0.10, code = 1,
       col = adjustcolor("purple", alpha.f = 0.8))

barplot(fem_cuum_age_total, names.arg = 10:52, col = adjustcolor("purple", alpha.f = 0.5),
        xaxt = "n",
        yaxt = "n", main = "Population Data", border = NA, xlab = "Age (yrs)") -> locs2
barplot(mal_cuum_age_total, names.arg = 10:52, col = adjustcolor("darkgreen", alpha.f = 0.5),
        xaxt = "n",
        yaxt = "n", add = TRUE, border = NA)

axis(1, at = c(locs2[seq(1, length(locs2), 5)]), labels = seq(10, 50, 5))
axis(2, las = 1)

text(41.5, 82, "excess male counts", col = adjustcolor("darkgreen", alpha.f = 0.8),
     cex = 0.75)

arrows(x0 = 28, y0 = 71, x1 = 31.4, y1 = 77,
       angle = 25, length = 0.10, code = 1,
       col = adjustcolor("darkgreen", alpha.f = 0.8))

legend("topright", pch = 15, col = c("#A020F080", "#00640080"),
       legend = c("Female", "Male"))

dev.off()
