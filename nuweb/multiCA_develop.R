library(devtools)
library(covr)
source('nuweb/Nuweb.R')
mc <- as.package(".")

use_git()
use_github()
use_version("minor")
use_author(given = "Aniko", family="Szabo", role = c("aut","cre"), email = "aszabo@mcw.edu")

nuweb(mc)
document(mc)
run_examples(mc)  # or dev_example("ran.CMData")
load_all(mc)

test(mc)
cov <- package_coverage(mc$path)
shine(cov)


check(mc, cran = TRUE, remote = TRUE)
install(mc)

# releasing to CRAN

release_checks(mc)
spell_check()

check_win_release()
check_win_devel()
check_mac_release()

rhub::rhub_setup() # run once
rhub::rhub_doctor()
rhub::rhub_check()


urlchecker::url_check()

tools::dependsOnPkgs("multiCA")
revdepcheck::cran_revdeps("multiCA")
revdepcheck::revdep_check()

use_cran_comments()  #run first time

# Final step:
release()


#create data set
strk <- data.matrix(read.delim("z:/EOGeorge/MultiTrend/StrokeData.txt", row.names=1))
colnames(strk) <- gsub("X", "", colnames(strk))
stroke <- as.data.frame.table(strk)
names(stroke) <- c("Type", "Year", "Freq")
stroke$Year <- as.numeric(as.character(stroke$Year))

use_data(stroke, pkg=mc)
promptData(stroke)
