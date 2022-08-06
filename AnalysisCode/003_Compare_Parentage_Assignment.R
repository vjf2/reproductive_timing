# Compare assignment rates
# Table 1 and Figure2

load("IntermediateData/sequoia_ped_with_sibs.RData")

parents <- ped_with_sibs$PedigreePar

parents$birthyear <- ped_with_sibs$LifeHist$BirthYear[match(parents$id,
                                                            ped_with_sibs$LifeHist$ID)]

lh <- read.csv("SharedData/LifeHistoryData.csv")

# Maternal interbirth interval is 4.25 overall, 4.97 surviving, so visualize 5-year cohorts

seq(1970, 2018, 5) |> length()

parents$birthcohort <- cut(parents$birthyear, 10)

dams <- aggregate(dam ~ birthcohort, parents, length)
sires <- aggregate(sire ~ birthcohort, parents, length)

births <- table(parents$birthcohort) |> as.data.frame()

names(births) <- c("birthcohort", "births")

parentage <- Reduce(merge, list(births, dams, sires))

parentage$damprop <- parentage$dam / parentage$births
parentage$sireprop <- parentage$sire / parentage$births

forplot <- t(parentage[, c("birthcohort", "damprop", "sireprop")])
colnames(forplot) <- forplot[1, ]
forplot <- forplot[-1, ]
mode(forplot) <- "numeric"

xlabnames <- gsub("\\(|\\]", "", colnames(forplot))
xlabnames <- gsub(",", "-", xlabnames)

windows()
# pdf("TablesFigures/Fig2_proportion_assigned.pdf", width = 9, height = 4)
x <- barplot(forplot, beside = TRUE,
             legend.text = c("Maternities", "Paternities"),
             col = adjustcolor(c("purple", "darkgreen"), alpha.f = 0.5),
             ylim = c(0, 1),
             xlab = "Offspring Birth Cohort",
             ylab = "Proportion Assigned",
             args.legend = list(x = 4, y = 0.98, bty = "n", cex = 1.25),
             xaxt = "n", yaxt = "n",
             space = c(0, 0.5))

axis(1, xlabnames, at = seq(1.5, 17.5, 2.5))
axis(2, las = 1)
text(x[1, ] + 0.5, y = forplot[1, ] + 0.07, parentage$births)

dev.off()

# Model
parents$sire01 <- ifelse(is.na(parents$sire), 0, 1)
parents$dam01 <- ifelse(is.na(parents$dam), 0, 1)

# Check for interaction (none found)
p2 <- parents[, c("dam01", "sire01", "birthyear", "id")]

p2 <- reshape(p2, direction = "long",
              varying = 1:2,
              v.names = "parent",
              times = c("dam", "sire"),
              timevar = "sex")

# Model individuals born after start of study (1984)
p2 <- p2[which(p2$birthyear >= 1984), ]

p2$sighting_rate <- lh$Sighting.Rate[match(p2$id, lh$Dolphin.ID)]

# No significant interaction between birth year and father assignment
# Or between sighting rate and father assignment
fullmod <- glm(parent ~ birthyear + sighting_rate + sex + sighting_rate:sex + birthyear:sex,
               data = p2,
               family = binomial)

summary(fullmod)

# No interactions so report simple model

pmod <- glm(parent ~ birthyear + sighting_rate + sex,
            data = p2,
            family = binomial)

summary(pmod)

library(sjPlot)

tab1 <- tab_model(pmod, show.ci = FALSE, show.stat = TRUE, show.se = TRUE,
                  show.r2 = TRUE, show.obs = FALSE, file = "TablesFigures/parent_assignment.html",
                  title = "Parentage Assignment", digits = 3)

tab1
