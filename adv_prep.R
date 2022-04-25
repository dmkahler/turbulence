## adv_prep.R ##

# This code pulls data from the Nortek Vector ADV and preps the endpoints for the profiling stops.  The 
# code will take both the .sen file with time and .dat file with pressure to determine the depth and 
# times.

## LOAD PACKAGES ##

library(readr)
library(dplyr)
library(ggplot2)
library(lubridate)
library(latex2exp)
library(gridExtra)
library(devtools)
install_github("LimpopoLab/hydrostats", force = TRUE)
library(hydrostats)

## READ DATA ##
fh <- "/Users/davidkahler/Documents/Hydrology_and_WRM/river_and_lake_mixing/ADV_data/109MON18" # filename header

# sen file
fn_sen <- paste(fh, "sen", sep = ".")
sen <- read_table(fn_sen, col_names = FALSE, col_types = "nnnnnnnnnnnnnnnn")
sen <- sen %>%
     rename(mon = X1, day = X2, yea = X3, hou = X4, mnt = X5, sec = X6, err = X7, sta = X8, bat = X9, ssp = X10, hed = X11, pit = X12, rol = X13, tmp = X14, a1 = X15, checksum = X16) %>%
     mutate(dt = ymd_hms(paste(yea,mon,day,hou,mnt,sec)))
# 1   Month                            (1-12)
# 2   Day                              (1-31)
# 3   Year
# 4   Hour                             (0-23)
# 5   Minute                           (0-59)
# 6   Second                           (0-59)
# 7   Error code
# 8   Status code
# 9   Battery voltage                  (V)
# 10   Soundspeed                       (m/s)
# 11   Heading                          (degrees)
# 12   Pitch                            (degrees)
# 13   Roll                             (degrees)
# 14   Temperature                      (degrees C)
# 15   Analog input
# 16   Checksum                         (1=failed)

# dat file
fn_dat <- paste(fh, "dat", sep = ".")
dat <- read_table(fn_dat, col_names = FALSE, col_types = "nndddnnnnnnnnnnnnn")
dat <- dat %>%
     rename(burst = X1, ensemble = X2, u = X3, v = X4, w = X5, amp1 = X6, amp2 = X7, amp3 = X8, snr1 = X9, snr2 = X10, snr3 = X11, corr1 = X12, corr2 = X13, corr3 = X14, p_dbar = X15, a1 = X16, a2 = X17, checksum = X18)
# 1   Burst counter
# 2   Ensemble counter                 (1-65536)
# 3   Velocity (Beam1|X|East)          (m/s)
# 4   Velocity (Beam2|Y|North)         (m/s)
# 5   Velocity (Beam3|Z|Up)            (m/s)
# 6   Amplitude (Beam1)                (counts)
# 7   Amplitude (Beam2)                (counts)
# 8   Amplitude (Beam3)                (counts)
# 9   SNR (Beam1)                      (dB)
# 10   SNR (Beam2)                      (dB)
# 11   SNR (Beam3)                      (dB)
# 12   Correlation (Beam1)              (%)
# 13   Correlation (Beam2)              (%)
# 14   Correlation (Beam3)              (%)
# 15   Pressure                         (dbar)        p_dbar
# 16   Analog input 1
# 17   Analog input 2
# 18   Checksum                         (1=failed)
sampling_rate = 64 # Hz, verify sampling rate in .hdr file under User setup

## DATA CHECK ##
records <- nrow(dat)
cover <- sampling_rate*nrow(sen) # manually compare this to records, should be close
error_checksum <- max(dat$checksum) # If there are any errors recorded, this will be = 1

## Check signal-to-noise ratio
snr1 <- ggplot(dat,aes(x=snr1)) +
     geom_histogram(breaks = (5*c(0:12)), color = "black", fill = "gray", na.rm = TRUE) +
     ylab("Beam 1") + 
     xlab(element_blank()) +
     theme(panel.background = element_rect(fill = "white", colour = "black")) + 
     theme(axis.text = element_text(face = "plain", size = 12))
snr2 <- ggplot(dat,aes(x=snr2)) +
     geom_histogram(breaks = (5*c(0:12)), color = "black", fill = "gray", na.rm = TRUE) +
     ylab("Beam 2") +
     xlab(element_blank()) +
     theme(panel.background = element_rect(fill = "white", colour = "black")) + 
     theme(axis.text = element_text(face = "plain", size = 12))
snr3 <- ggplot(dat,aes(x=snr3)) +
     geom_histogram(breaks = (5*c(0:12)), color = "black", fill = "gray", na.rm = TRUE) +
     ylab("Beam 3") +
     xlab("Signal-to-Noise Ratio (dB)") +
     theme(panel.background = element_rect(fill = "white", colour = "black")) + 
     theme(axis.text = element_text(face = "plain", size = 12))
snr_tot <- grid.arrange(snr1,snr2,snr3, nrow = 3) # SNR should be higher, but noise is common when out of water.  Compare improvement to next figure.

## CONVERT DATA ##
dat$dt <- (c(1:records))/sampling_rate # time is defined as seconds after start, 
dat$p_Pa <- 1e4*dat$p_dbar
p_atm <- mean(dat$p_Pa[1:(10*sampling_rate)]) # atmospheric pressure in Pa, assuming instrument started above water
under <- which(dat$p_Pa>(p_atm+5e2)) # finds all indexes where the instrument is deeper than 5cm
st <- 1000*(floor((under[5*64])/1000)) # takes the start 10 seconds into, then backs it out; this removes random single values
en <- 1000*ceiling((under[(length(under)-(5*64))])/1000)
dat2 <- dat[st:en,]

## Check signal-to-noise ratio
snr1 <- ggplot(dat2,aes(x=snr1)) +
     geom_histogram(breaks = (5*c(0:12)), color = "black", fill = "gray", na.rm = TRUE) +
     ylab("Beam 1") + 
     xlab(element_blank()) +
     theme(panel.background = element_rect(fill = "white", colour = "black")) + 
     theme(axis.text = element_text(face = "plain", size = 12))
snr2 <- ggplot(dat2,aes(x=snr2)) +
     geom_histogram(breaks = (5*c(0:12)), color = "black", fill = "gray", na.rm = TRUE) +
     ylab("Beam 2") +
     xlab(element_blank()) +
     theme(panel.background = element_rect(fill = "white", colour = "black")) + 
     theme(axis.text = element_text(face = "plain", size = 12))
snr3 <- ggplot(dat2,aes(x=snr3)) +
     geom_histogram(breaks = (5*c(0:12)), color = "black", fill = "gray", na.rm = TRUE) +
     ylab("Beam 3") +
     xlab("Signal-to-Noise Ratio (dB)") +
     theme(panel.background = element_rect(fill = "white", colour = "black")) + 
     theme(axis.text = element_text(face = "plain", size = 12))
snr_all <- grid.arrange(snr1,snr2,snr3, nrow = 3)
ggsave("snr_all.eps", snr_tot, device = "eps", dpi = 72)

## Big Profile
uall <- ggplot(dat2) +
     geom_line(aes(x=dt,y=u)) +
     ylab("u (m/s)") +
     xlab(element_blank()) +
     theme(panel.background = element_rect(fill = "white", colour = "black")) + 
     theme(axis.text = element_text(face = "plain", size = 12))
vall <- ggplot(dat2) +
     geom_line(aes(x=dt,y=v)) +
     ylab("v (m/s)") +
     xlab(element_blank()) +
     theme(panel.background = element_rect(fill = "white", colour = "black")) + 
     theme(axis.text = element_text(face = "plain", size = 12))
wall <- ggplot(dat2) +
     geom_line(aes(x=dt,y=w)) +
     ylab("w (m/s)") +
     xlab(element_blank()) +
     theme(panel.background = element_rect(fill = "white", colour = "black")) + 
     theme(axis.text = element_text(face = "plain", size = 12))
pall <- ggplot(dat2) +
     geom_line(aes(x=dt,y=(p_Pa/1000))) +
     ylab("Pressure (kPa)") +
     xlab("Time (s)") +
     theme(panel.background = element_rect(fill = "white", colour = "black")) + 
     theme(axis.text = element_text(face = "plain", size = 12))
guall <- ggplotGrob(uall)
gvall <- ggplotGrob(vall)
gwall <- ggplotGrob(wall)
gpall <- ggplotGrob(pall)
grid::grid.newpage()
uvwp <- grid::grid.draw(rbind(guall,gvall,gwall,gpall))
ggsave("uvwp.eps", uvwp, device = "eps", dpi = 72)

# average pressures
pressure <- array(NA, dim = c(dur,sampling_rate))
pres <- array(NA, dim = dur)
temp <- pres
dnst <- pres
dept <- pres
time <- pres
for (i in 1:dur) {
     for (j in 1:sampling_rate) {
          pressure[i,j] <- dat$p_Pa[((i-1)*sampling_rate+j)]
     }
     pres[i] <- mean(pressure[i,])  # pressure (Pa)
     temp[i] <- sen$tmp[i]          # temperature (C)
     dnst[i] <- waterrho(temp[i])   # water density (kg/m^3)
     time[i] <- as.numeric(sen$dt[i]) - start
}

num <- 0                               # set index for breaks
breaks <- array(0, dim = dur)          # preallocate a matrix to mark breaks
for (i in 12:dur) {
     ave_pres <- mean(pres[(i-11):(1-1)])
     if ((pres[i] - pres[i-1]) > (10 * ave_pres)) {
          num <- num + 1             # advance index
          breaks[i] <- 1             # record break
     }
}

loc <- data.frame(time,temp,pres,dnst) # time in seconds from start, temp in Celcius, pres in Pascals, dnst in N/m^3

ggplot(loc) +
     geom_line(aes(x=time, y=pres)) +
     theme(panel.background = element_rect(fill = "white", colour = "black")) +
     theme(aspect.ratio = 1) +
     theme(axis.text = element_text(face = "plain", size = 12))

## FIND ATMOSPHERIC PRESSURE

base <- mean(loc$pres[1:10]) # pressure for the first 10 seconds in Pa
b.sd <- stdev(loc$pres[1:10])



