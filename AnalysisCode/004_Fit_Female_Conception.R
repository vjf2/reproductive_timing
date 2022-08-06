# Create Figure 3

# Proportion of females conceiving a calf that survives to at least two years of age at each maternal age (n = 510 calves), compared with the proportion of just the calves that were assigned a mother in the genetic parentage assignment procedure (n = 141). Histogram represents all observed mother-calf pairs, with gaussian kernel density estimates (KDE) and a fitted gamma distribution overlaid 

source("AnalysisCode/get_age_function.R")

lh <- read.csv("SharedData/LifeHistoryData.csv")
lh$Birth.Date <- as.Date(lh$Birth.Date)

lh$momage <- get_age(data = lh,
                     column = "Mother.ID",
                     atDate = "Birth.Date",
                     lhdata = lh)

lh$Last.Sighting.Date <- as.Date(lh$Last.Sighting.Date)

lh$calf_age_lsd <- get_age(data = lh,
                           column = "Dolphin.ID",
                           atDate = "Last.Sighting.Date",
                           lhdata = lh)


sur_calf <- lh[-which(lh$calf_age_lsd <= 2 |
                        is.na(lh$calf_age_lsd)), ]
sur_calf <- sur_calf[which(sur_calf$Mother.ID != ""), ]

sur_calf$conception_age <- sur_calf$momage - (53 / 52) # 53 weeks gestation

ages <- c(na.omit(sur_calf$conception_age))

sampled_ages <- c(na.omit(sur_calf$conception_age[which(sur_calf$Genotype.Available == "Y")]))

t.test(ages, sampled_ages)

# Fit density
library(fitdistrplus)

descdist(ages, discrete = FALSE)

fit.gamma <- fitdist(ages, "gamma")

plot(fit.gamma)

# shape, rate 8.95 0.49

windows()
# pdf(file="TablesFigures/Fig3_observed_biopsied_maternities.pdf")
hist(sur_calf$conception_age, prob = TRUE, breaks = seq(5, 45, 2.5),
     ylab = "Proportion", xlab = "Maternal age at Conception",
     border = "darkgrey", main = NA, ylim = c(0, 0.1),
     # xlim = c(5, 45),
     yaxt = "n",
     xaxt = "n")

axis(1, at = seq(5, 45, 5))
axis(2, las = 1)

lines(density(na.omit(sur_calf$conception_age), adjust = 1.5), col = "black", lwd = 2)

lines(density(sampled_ages, adjust = 1.5), col = "purple", lwd = 2)

lines(7:41, dgamma(7:41, shape = 8.957, rate = 0.495), col = "red", lty = 2, lwd = 2)

legend("topright", legend = c("KDE Observed", "KDE Genetic Assignment", "Gamma Fit"),
       col = c("black", "purple", "red"),
       lty = c(1, 1, 2), lwd = c(2, 2, 2))


dev.off()

# Confirm gamma fit

x <- rgamma(n = length(ages), shape = 8.957, rate = 0.495)
t.test(x, ages)
