#install.packages("ggplot2")
library("ggplot2")
library("readr")


#change this to the relative path to the blastResults folder in your own implementation of the project !!!
setwd("~/Desktop/Results/Raw_results_data_all_teams/blastResults")

dir.create("./../Rdata")



frequencyAnalysisPlots <- function(blastOutput){

csvFile <- read_table2(blastOutput,
                                col_names = FALSE,
                                comment = "#")




names(csvFile) <- c("Query_id",
                            "Subject_id",
                            "%_identity",
                            "alignment_length",
                            "mismatches",
                            "gap_openings",
                            "q.start",
                            "q.end",
                            "s.start",
                            "s.end",
                            "e-value",
                            "bit_score")

# I need to filter out matches between identical contigs to remove
# the biggest 100% coverage hits, or this comparison is meaningless.
# given that I alreayd KNOW two identical sequences
# will be... identical.

# gets rids of all fully identical of Seq A vs Seq A, etc...
csvFile2 <- csvFile[csvFile$Query_id != csvFile$Subject_id,]


filepath = paste("./../Rdata/", blastOutput, ".png", sep="")


png(filename= filepath)

par(mfrow=c(2,3),oma=c(0,0,2,0))

hist(csvFile2$alignment_length, main = "alignment_length")
abline(v=mean(csvFile2$alignment_length),col="blue")
hist(csvFile2$mismatches, main = "mismatches")
hist(csvFile2$gap_openings, main = "gap_openings")
hist(csvFile2$`e-value`, main = "e-value")
hist(csvFile2$bit_score, main = "bit_score")
hist(csvFile2$`%_identity`, main = "%_identity")
# creates the title of the big plot
title(blastOutput, outer=TRUE)

# closes off the png file.
dev.off()

# save the result of the summary function a textfile
sumFile <- paste("./../Rdata/", blastOutput, "_summary.txt")

sink(sumFile)
print(summary(csvFile2))
sink()


}
# Pattern matching to get only "autoblasts"
list <- list.files(path = getwd())  # list.files(path = getwd(), pattern= '(.*)_vs_\\1.csv$')

for(i in seq_along(list)){
  frequencyAnalysisPlots(list[i])
}
