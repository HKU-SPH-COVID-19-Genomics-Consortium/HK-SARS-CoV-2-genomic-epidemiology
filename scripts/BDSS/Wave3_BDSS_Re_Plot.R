devtools::install_github("laduplessis/bdskytools", force=TRUE)
library(bdskytools)
library(lubridate)
library(ggplot2)

## Read in BEAST2 logfile for Wave-3 BDSS, discarding initial 10% of samples as burn-in
w3.fname <- "./Wave3_BDSS.log"
w3.lf <- readLogfile(w3.fname, burnin=0.1)

## Extract reproductiveNumber from logfile content
w3.Re_sky <- getSkylineSubset(w3.lf, "reproductiveNumber")

## Extract 95% HPD intervals of reproductiveNumber
w3.Re_hpd <- getMatrixHPD(w3.Re_sky)

## Set number of time-points at which reproductiveNumber was estimated in BDSS
w3.int_num <- 12

## Set time-gridpoints (n=27) from estimated treeHeight to time of most recent sample in Wave-3 dataset
w3.timegrid <- seq(0, 0.3054, length.out=27) # Estimated treeHeight as extracted from logfile using Tracer v1.7.1
w3.Re_gridded <- gridSkyline(w3.Re_sky, w3.lf$origin, w3.timegrid)
w3.Re_gridded_hpd <- getMatrixHPD(w3.Re_gridded)

## Invert time-gridpoints to plot temporal changes of reproductiveNumber from past to present
w3.recent <- 2020.8032786885246 # Date of most recent sample in Wave-3 dataset in decimal format as calculated using TempEst v1.5.3
w3.times <- w3.recent - w3.timegrid

# Transform function for time-axis labelling
x.date_transform <- function(x) {format(date_decimal(x), "%d/%b")}

## Plot temporal changes in reproductiveNumber
ggplot() +
  labs(x='Time (yr)', y='Effective reproductive number, Re') +
  scale_x_continuous(labels=x.date_transform, breaks=seq(w3.times[length(w3.times)], w3.times[1], by=0.019165), expand = c(0, 0)) + ## x-label every 1 week (0.019165 year)
  geom_ribbon(aes(x=w3.times, ymin=w3.Re_gridded_hpd[3,], ymax=w3.Re_gridded_hpd[1,]), fill="#DBD0DD") +
  geom_line(aes(x=w3.times, y=w3.Re_gridded_hpd[2,]), size=0.3, col="#654C6B") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  