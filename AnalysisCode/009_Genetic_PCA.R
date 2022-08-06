library(adegenet)
library(vegan)

lh <- read.csv("SharedData/LifeHistoryData.csv")
rel <- read.csv("SharedData/female_relatedness_coefs.csv")

all_names <- c(rel$ID1, rel$ID2) |> unique()

good_mom <- all_names

# Remove one of each closely related pair

for (i in all_names) {

  my_rel <- rel[grep(i, rel$pair), ]

  if (any(my_rel$DyadML >= 0.1923)) {good_mom <- setdiff(good_mom, i)}

  rel <- rel[which(rel$ID1 %in% good_mom &
                                  rel$ID2 %in% good_mom), ]
}

# Select offspring

offspring <- lh$Dolphin.ID[which(lh$Mother.ID %in% good_mom &
                                   lh$Genotype.Available == "Y")]

# Load in genotypes and convert to genind object

loci <- read.table("SharedData/genotypes.txt", row.names = 1)

loci <- loci[which(rownames(loci) %in% c(good_mom, offspring)), ]

x.mat <- as.matrix(loci)
x.mat[x.mat == 0] <- "1/1" # homozygote reference
x.mat[x.mat == 1] <- "1/2" # heterozygote
x.mat[x.mat == 2] <- "2/2" # homozygote alternate
x.mat[x.mat == -9] <- "" # missing data
x.gid <- df2genind(x.mat, ind.names = rownames(x.mat), sep = "/", ploidy = 2)

mom.gid.X <- tab(x.gid[good_mom], freq = TRUE, NA.method = "mean")
off.gid.X <- tab(x.gid[offspring], freq = TRUE, NA.method = "mean")

x.gid.pca1 <- dudi.pca(mom.gid.X, scale = TRUE, scannf = FALSE, nf = 2)

ind.sup.coord <- suprow(x.gid.pca1, off.gid.X)$lisup

#Load in dad info
load("IntermediateData/sequoia_ped_with_sibs.RData")

ped <- ped_with_sibs$PedigreePar
ind.sup.coord$dad <- ped$sire[match(rownames(ind.sup.coord), ped$id)]

ind.sup.coord$good_dad <- ifelse(!is.na(ind.sup.coord$dad),
                                 "yes", "no")

head(ind.sup.coord)

ind.sup.coord$color <- ifelse(ind.sup.coord$good_dad == "yes", "green4", "red4")

# Plot
ind.sup.coord$pch <- ifelse(ind.sup.coord$good_dad == "yes", 16, 1)

windows()
# pdf(file="TablesFigures/Fig7_offspring_pca.pdf")
plot(x.gid.pca1$li$Axis1, x.gid.pca1$li$Axis2, type = "n", asp = 1,
     xlab = "PC1", ylab = "PC2", yaxt = "n")

axis(2, las = 1)
# moms
points(x.gid.pca1$li$Axis1, x.gid.pca1$li$Axis2, pch = 2, cex = 0.75, col = "grey")
# calves
points(ind.sup.coord$Axis1, ind.sup.coord$Axis2, col = ind.sup.coord$color, pch = ind.sup.coord$pch,
       cex = 1.25)

abline(h = 0)
abline(v = 0)

df <- ind.sup.coord
find_hull <- function(df) df[chull(df[, 1], df[, 2]), ]
hulls <- plyr::ddply(df, "good_dad", find_hull)

polygon(hulls[hulls$good_dad == "no", ], col = adjustcolor("red4", alpha.f = 0.2))
polygon(hulls[hulls$good_dad == "yes", ], col = adjustcolor("green4", alpha.f = 0.2))


legend("topright", pch = c(1, 16, 2), legend = c("No father assigned", "Father assigned", "Mother"),
       col = c("red4", "green4", "grey"),
       pt.cex = c(1.25, 1.25, 0.75))

dev.off()

# Statistical test

# Select just the offspring
off.gid <- x.gid[which(rownames(x.mat) %in% rownames(ind.sup.coord)), ]

dol_dist <- dist(off.gid)

dads <- ind.sup.coord$good_dad[match(rownames(x.mat), rownames(ind.sup.coord))]
dads <- dads[!is.na(dads)]

strata(off.gid) <- data.frame(dads)
dol_stra <- strata(off.gid)

set.seed(286567440)
res <- vegan::adonis(dol_dist ~ dads, data = dol_stra, permutations = 999)
