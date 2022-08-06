# Select best fitting paternal age
load("IntermediateData/simulation_results.RData")

# Extract best fits
# Constrain to only simulations which output a mean sampled paternal age within
# +/- 3 years of observed paternal age

sim_agg$mean_realdad_age_sampled[1]

sim_agg_realistic <- sim_agg[which(sim_agg$mean_randad_age_sampled >= 20.3 & sim_agg$mean_randad_age_sampled <= 26.3), ]

range(sim_agg_realistic$ranparent_ratio)

top <- 95

res <- sim_agg_realistic[which(sim_agg_realistic$parent_ratio_diff %in%
                                 sort(abs(sim_agg_realistic$parent_ratio_diff))[1:95]), ]

res1 <- sim_agg_realistic[which.min(abs(sim_agg_realistic$parent_ratio_diff)), ]

colr <- c("grey", rep("grey", nrow(res) - 1))

windows()
# pdf("TablesFigures/Fig5_relative_probabilities_werror.pdf")
plot(10:50, dgamma(10:50, shape = 9, rate = 0.5), type = "n", ylim = c(0, 0.08),
     col = "darkgreen",
     xlab = "Age at Conception",
     ylab = "Proportion", lwd = 2,
     yaxt = "n")

axis(2, las = 1)

for (i in 1:nrow(res)) {
  lines(10:50, dgamma(10:50, shape = res[i, 1], rate = res[i, 2]), col =
          adjustcolor(colr[i], alpha.f = 0.7), lwd = 1)

} # set of new male best fits

lines(10:50,
      dgamma(10:50, shape = res1[1, 1], rate = res1[1, 2]),
      col = "darkgreen", lwd = 2)

lines(10:50, dgamma(10:50, shape = 9.0, rate = 0.5), col = "purple", lwd = 2) # female best fit


legend("topright", legend = c("Female", "Male Best Fit"),
       col = c("purple", "darkgreen"), lty = 1, lwd = 2)

dev.off()

dgamma(1:50, shape = 22.4, rate = 0.71) |> which.max()

# Plot against observed data

sampled_ages <- offspring$realdad_age[!is.na(offspring$realdad_age)]

windows()
# pdf(file="TablesFigures/Fig4_predicted_biopsied_paternities.pdf")
hist(sampled_ages, prob = TRUE, breaks = seq(10, 40, 2.5),
     ylab = "Proportion", xlab = "Paternal Age at Conception",
     border = "darkgrey", main = NA, ylim = c(0, 0.08),
     yaxt = "n")

lines(density(sampled_ages, adjust = 1.5), col = "darkgreen", lwd = 2)

axis(2, las = 1)

lines(10:41, dgamma(10:41, shape = 22.4, rate = 0.71), col = "red", lty = 2, lwd = 2)

legend("topright", legend = c("KDE Genetic Assignment", "Gamma Fit from Simulation"),
       col = c("darkgreen", "red"),
       lty = c(1, 2), lwd = c(2, 2))

dev.off()

# Compare gamma fit

x <- rgamma(n = length(sampled_ages), shape = 22.4, rate = 0.71)
t.test(x, sampled_ages)

# write.csv(sim_agg, "SharedData/simulation_results.csv")
