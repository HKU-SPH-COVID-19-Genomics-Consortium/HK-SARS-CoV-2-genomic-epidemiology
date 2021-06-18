devtools::install_github("laduplessis/bdskytools", force=TRUE)
library(bdskytools)
library(lubridate)
library(ggplot2)

## Read in BEAST2 logfile for Wave-4 BDSS, discarding initial 10% of samples as burn-in
w4.fname <- "./Wave4_BDSS.log"
w4.lf <- readLogfile(w4.fname, burnin=0.1)

## Extract reproductiveNumber from logfile content
w4.Re_sky <- getSkylineSubset(w4.lf, "reproductiveNumber")

## Extract 95% HPD intervals of reproductiveNumber
w4.Re_hpd <- getMatrixHPD(w4.Re_sky)

## Set number of time-points at which reproductiveNumber was estimated in BDSS
w4.int_num <- 15

## Set time-gridpoints (n=27) from estimated treeHeight to time of most recent sample in Wave-4 dataset
w4.timegrid <- seq(0, 0.3855, length.out=35) # Estimated treeHeight as extracted from logfile using Tracer v1.7.1
w4.Re_gridded <- gridSkyline(w4.Re_sky, w4.lf$origin, w4.timegrid)
w4.Re_gridded_hpd <- getMatrixHPD(w4.Re_gridded)

## Invert time-gridpoints to plot temporal changes of reproductiveNumber from past to present
w4.recent <- 2021.0684931506848 # Date of most recent sample in Wave-4 dataset in decimal format as calculated using TempEst v1.5.3
w4.times <- w4.recent - w4.timegrid

# Transform function for time-axis labelling
x.date_transform <- function(x) {format(date_decimal(x), "%d/%b")}

## Plot temporal changes in reproductiveNumber
ggplot() +
  labs(x='Time (yr)', y='Effective reproductive number, Re') +
  scale_x_continuous(labels=x.date_transform, breaks=seq(w4.times[length(w4.times)], w4.times[1], by=0.019165), expand = c(0, 0)) + ## x-label every 1 week (0.019165 year)
  geom_ribbon(aes(x=w4.times, ymin=w4.Re_gridded_hpd[3,], ymax=w4.Re_gridded_hpd[1,]), fill="#FFECD5") +
  geom_line(aes(x=w4.times, y=w4.Re_gridded_hpd[2,]), size=0.3, col="#E27F10") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
