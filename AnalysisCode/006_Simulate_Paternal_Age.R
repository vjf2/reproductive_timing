# Simulate random parentages for each offspring, where females are assigned according to
# their known probabilities, and male ages at conception are varied

load("IntermediateData/simulation_setup.RData")

# Generate 20 random iterations of each parameter set

iteration <- 1:20
shape.seq <- seq(1, 30, 0.1)
rate.seq <- seq(0.1, 1.2, 0.01)

params <- expand.grid(iteration, shape.seq, rate.seq)

n <- nrow(params)

sim_res <- data.frame(iteration = params[, 1],
                      shape = params[, 2],
                      rate = params[, 3],
                      ranmom_biopsied = numeric(n),
                      randad_biopsied = numeric(n),
                      ranparent_ratio = numeric(n),
                      mean_ranmom_age = numeric(n),
                      mean_randad_age = numeric(n),
                      mean_ranmom_age_sampled = numeric(n),
                      mean_randad_age_sampled = numeric(n))

# Constrain so at least 95% of conceptions are between ages 9 and 55
sim_res$pgamma <- pgamma(55, sim_res$shape, sim_res$rate) - pgamma(9, sim_res$shape, sim_res$rate)

sim_res <- sim_res[which(sim_res$pgamma >= 0.95), ]

# Execture on cluster, ~ 4 hour runtime on Intel(R) Core(TM) i7-8565U CPU @ 1.80GHz 

library(parallel)
library(foreach)
library(doParallel)
library(doRNG)

cl <- makeCluster(detectCores() - 1, outfile = "test_out.txt")
clusterExport(cl, c("offspring", "moms", "dads", "sampled", "sim_res"))
registerDoParallel(cl)
registerDoRNG(seed = 286567440)

starttime <- Sys.time()

trial1 <- foreach(j = 1:nrow(sim_res), .errorhandling = "pass") %dopar% {

  for (i in 1:nrow(offspring)) {

    baby <- offspring$Dolphin.ID[i]
    conception_year <- offspring$Conception.Year[i]

    p_mom <- moms[as.character(conception_year),
                  which(moms[as.character(conception_year), ] != 0)]

    mom_prob <- dgamma(p_mom, shape = 9.0, rate = 0.5)

    p_dad <- dads[as.character(conception_year),
                  which(dads[as.character(conception_year), ] != 0)]

    dad_prob <- dgamma(p_dad, shape = sim_res[j, 2], rate = sim_res[j, 3])

    offspring$ranmom[i] <- sample(names(p_mom), 1, prob = mom_prob)
    offspring$randad[i] <- sample(names(p_dad), 1, prob = dad_prob)

  }


  # ranmom and randad age at conception
  offspring$ranmom_age <- get_age(data = offspring,
                                  column = "ranmom",
                                  lhdata = lh)

  offspring$randad_age <- get_age(data = offspring,
                                  column = "randad",
                                  lhdata = lh)

  # sample based on ages of real samples
  offspring$ranmom_biopsied <- ifelse(offspring$ranmom %in% sampled, 1, 0)
  offspring$randad_biopsied <- ifelse(offspring$randad %in% sampled, 1, 0)

  sim_res[j, "ranmom_biopsied"] = sum(offspring$ranmom_biopsied)
  sim_res[j, "randad_biopsied"] = sum(offspring$randad_biopsied)
  sim_res[j, "ranparent_ratio"] = sum(offspring$ranmom_biopsied) / sum(offspring$randad_biopsied)
  sim_res[j, "mean_ranmom_age"] = mean(offspring$ranmom_age)
  sim_res[j, "mean_randad_age"] = mean(offspring$randad_age)
  sim_res[j, "mean_ranmom_age_sampled"] = mean(offspring$ranmom_age[offspring$ranmom_biopsied == 1])
  sim_res[j, "mean_randad_age_sampled"] = mean(offspring$randad_age[offspring$randad_biopsied == 1])

  if (j %% 1000 == 0) {print(j)}

  return(sim_res[j, ])

}

endtime <- Sys.time()

stopCluster(cl)

endtime - starttime

#############################
# Compare to observed data

sim_res <- do.call("rbind", trial1)

sim_agg <- aggregate(. ~ shape + rate, sim_res, mean)

load("IntermediateData/sequoia_ped_with_sibs2.RData")

offspring$realmom <- lh$Mother.ID[match(offspring$Dolphin.ID, lh$Dolphin.ID)]
offspring$realdad <- lh$Father.ID[match(offspring$Dolphin.ID, lh$Dolphin.ID)]

offspring$realdad <- ped_with_sibs$PedigreePar$sire[match(offspring$Dolphin.ID,
                                                          ped_with_sibs$PedigreePar$id)]

offspring$realmom_biopsied <- ifelse(offspring$realmom %in% sampled, 1, 0)
offspring$realdad_biopsied <- ifelse(offspring$realdad %in% sampled, 1, 0)

offspring$realmom_age <- get_age(data = offspring,
                                 column = "realmom",
                                 lhdata = lh)


offspring$realdad_age <- get_age(data = offspring,
                                 column = "realdad",
                                 lhdata = lh)

sim_agg[, "realmom_biopsied"] = sum(offspring$realmom_biopsied)
sim_agg[, "realdad_biopsied"] = sum(offspring$realdad_biopsied)
sim_agg[, "realparent_ratio"] = sum(offspring$realmom_biopsied) / sum(offspring$realdad_biopsied)
sim_agg[, "mean_realmom_age"] = mean(offspring$realmom_age, na.rm = TRUE)
sim_agg[, "mean_realdad_age"] = mean(offspring$realdad_age, na.rm = TRUE)
sim_agg[, "mean_realmom_age_sampled"] = mean(offspring$realmom_age[offspring$realmom_biopsied == 1])
sim_agg[, "mean_realdad_age_sampled"] = mean(offspring$realdad_age[offspring$realdad_biopsied == 1])

sim_agg$parent_ratio_diff <- sim_agg$realparent_ratio - sim_agg$ranparent_ratio

sim_agg$mean_age_diff <- sim_agg$mean_realdad_age_sampled - sim_agg$mean_randad_age_sampled

# save(sim_agg, offspring, file="IntermediateData/simulation_results.RData")

