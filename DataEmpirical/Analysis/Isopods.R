library(nlme)
setwd("~/MutationSelectionDrift/DataEmpirical/Analysis/DataFrame")
data <- read.table("Ncat50_Isopods_rootedtree.nhx_cds.highcoverage.list_Sample12_Replicates6_merged_node_False_1_run.LogPopulationSize.tsv",h=TRUE, sep = "\t")
names(data)[1] <- "PopulationSize"

print(aov(PopulationSize~Ocular.structure+Rep, data))

mosaicplot(table(data$Habitat, data$Pigmentation), shade = TRUE)
mosaicplot(table(data$Habitat, data$Ocular.structure), shade = TRUE)
mosaicplot(table(data$Ocular.structure, data$Pigmentation), shade = TRUE)

out <- anova(lm(PopulationSize~Ocular.structure+Rep, data))
print(out)
cat(capture.output(out), file="anova_ocular_structure.txt", sep="\n", append=FALSE)

out <- anova(lm(PopulationSize~Pigmentation+Rep, data))
print(out)
cat(capture.output(out), file="anova_pigmentation.txt", sep="\n", append=FALSE)

print(summary(lm(PopulationSize~Habitat+Rep, data)))
print(summary(lm(PopulationSize~Habitat+Ocular.structure+Pigmentation, data)))
out <- anova(lm(PopulationSize~Habitat+Rep, data))
print(out)
cat(capture.output(out), file="anova_habitat.txt", sep="\n", append=FALSE)

print(summary(lm(PopulationSize~Habitat, data)))
print(summary(lm(PopulationSize~Pigmentation, data)))
print(summary(lm(PopulationSize~Ocular.structure, data)))

