#                 09.4.AgglomerativeHierarchical Clustering Statistical Study - LICA
#       Author: Alexperezm | Master's End of Degree Project - 2021-2022

#Objective: Statistical analysis of the clustering results developed by the
#hierarchical Clustering Model.

library(readxl)
results_HC <- read_excel("09.04.LICA_survival.xlsx", 
                                   na = "NA")

boxplot(survival ~ S, data=results_HC)
boxplot(survival ~ M, data=results_HC)
boxplot(survival ~ L, data=results_HC)

salida.aov = aov(survival ~ S, data=results_HC)
summary(salida.aov)

salida.aov = aov(survival ~ M, data=results_HC)
summary(salida.aov)

salida.aov = aov(survival ~ L, data=results_HC)
summary(salida.aov)

with(results_HC,pairwise.t.test(survival,S,p.adj="bonferroni"))
with(results_HC,pairwise.t.test(survival,M,p.adj="bonferroni"))
with(results_HC,pairwise.t.test(survival,L,p.adj="bonferroni"))
