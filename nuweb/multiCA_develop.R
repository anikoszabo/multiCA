library(devtools)
library(covr)
source('c:/RForge/Nuweb.R')
mc <- as.package("c:/RForge/multiCA")

use_git()
use_github()

nuweb(mc)
document(mc)
run_examples(mc)  # or dev_example("ran.CMData")
load_all(mc)

test(mc)
cov <- package_coverage(mc$path)
shine(cov)


check(mc, check_dir = "c:/RForge", check_version = TRUE, cran = TRUE, manual=TRUE)
install(mc)

#create data set
strk <- data.matrix(read.delim("z:/EOGeorge/MultiTrend/StrokeData.txt", row.names=1))
colnames(strk) <- gsub("X", "", colnames(strk))
stroke <- as.data.frame.table(strk)
names(stroke) <- c("Type", "Year", "Freq")
stroke$Year <- as.numeric(as.character(stroke$Year))

use_data(stroke, pkg=mc)
promptData(stroke)
