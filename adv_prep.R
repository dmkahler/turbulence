## adv_prep.R ##

# This code pulls data from the Nortek Vector ADV and preps the endpoints for the profiling stops.  The 
# code will take both the .sen file with time and .dat file with pressure to determine the depth and 
# times.

## LOAD PACKAGES ##

library(readr)
library(dplyr)
library(tidyr)
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

## CONVERT DATA ##
dat$dt <- (c(1:records))/sampling_rate # time is defined as seconds after start, 
dat$p_Pa <- 1e4*dat$p_dbar
p_atm <- mean(dat$p_Pa[1:(10*sampling_rate)]) # atmospheric pressure in Pa, assuming instrument started above water
under <- which(dat$p_Pa>(p_atm+5e2)) # finds all indexes where the instrument is deeper than 5cm
st <- 1000*(floor((under[5*64])/1000)) # takes the start 10 seconds into, then backs it out; this removes random single values
en <- 1000*ceiling((under[(length(under)-(5*64))])/1000)
dat2 <- dat[st:en,]

## Check signal-to-noise ratio
# FIGURE 1
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
ggsave("snr_all.eps", snr_all, device = "eps", dpi = 72)
rm(snr1,snr2,snr3)

## Big Profile
# FIGURE 2
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
guall <- ggplotGrob(uall)
gvall <- ggplotGrob(vall)
gwall <- ggplotGrob(wall)
uvw <- grid::grid.draw(rbind(guall,gvall,gwall))
ggsave("uvw.eps", uvw, device = "eps", dpi = 72)
rm(uall,vall,wall,guall,gvall,gwall)

## Find depth changes
n <- nrow(dat2)
change <- 0.05*997*9.81 # establish the threshold for if the instrument is moving vertically = 5cm movement
breaks <- array(NA, dim = n)
j <- 1
last.break <- -1 * sampling_rate
for (i in sampling_rate:(n-sampling_rate)) {
     fst <- mean(dat2$p_Pa[(i-63):i], na.rm = TRUE)
     sec <- mean(dat2$p_Pa[i:(i+63)], na.rm = TRUE)
     # dpdt[i] <- abs(sec-fst) # was for development, may delete this line
     if ( (abs(sec-fst) > change) & (i > (last.break+(5*sampling_rate))) ) {
          breaks[j] <- i
          last.break <- i
          j <- j + 1
     }
}
breaks <- breaks[which(is.na(breaks)==FALSE)] # INDEX OF dat2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
breaks.time <- min(dat2$dt) + breaks/64

# FIGURE 3
snr_tab <- dat2 %>%
     select(dt,snr1,snr2,snr3) %>%
     pivot_longer(cols = c("snr1","snr2","snr3"), names_to = "beam", values_to = "snr")
snr <- ggplot(snr_tab) +
     geom_line(aes(x=dt, y=snr, color=beam)) +
     ylab("Signal-to-Noise (dB)") +
     theme(axis.title.x=element_blank(),
           axis.text.x=element_blank()) +
     theme(legend.position="top", legend.key = element_blank()) + 
     theme(panel.background = element_rect(fill = "white", colour = "black")) + 
     theme(axis.text = element_text(face = "plain", size = 12))
pall <- ggplot(dat2) +
     geom_line(aes(x=dt,y=(p_Pa/1000))) +
     geom_vline(xintercept = breaks.time, color = "blue") +
     ylab("Pressure (kPa)") +
     xlab("Time (s)") +
     theme(panel.background = element_rect(fill = "white", colour = "black")) + 
     theme(axis.text = element_text(face = "plain", size = 12))
gsnr <- ggplotGrob(snr)
gpall <- ggplotGrob(pall)
grid::grid.newpage()
snrp <- grid::grid.draw(rbind(gsnr,gpall))
ggsave("snrp.eps", snrp, device = "eps", dpi = 72)
rm(snr,pall,gsnr,gpall)

## IDENTIFY TIMES/INDEX
duration <- array(NA, dim = (length(breaks.time)))
start.index <- duration
start.time <- duration
end.index <- duration
end.time <- duration
elapsed <- duration
buffer <- 2 # time in seconds removed from the start and end of each
j <- 1
for (i in 2:length(breaks.time)) {
     duration[i] <- breaks.time[i] - breaks.time[i-1]
     if (duration[i] >= 59) {
          start.index[j] <- breaks[i-1] + (buffer * sampling_rate)
          start.time[j] <- breaks.time[i-1] + buffer
          end.index[j] <- breaks[i] - (buffer * sampling_rate)
          end.time[j] <- breaks.time[i] - buffer
          elapsed[j] <- duration[i] - (2 * buffer)
          j <- j + 1
     }
}
start.index <- start.index[which(is.na(start.index)==FALSE)]
start.time <- start.time[which(is.na(start.time)==FALSE)]
end.index <- end.index[which(is.na(end.index)==FALSE)]
end.time <- end.time[which(is.na(end.time)==FALSE)]
elapsed <- elapsed[which(is.na(elapsed)==FALSE)]
depths <- data.frame(start.index,end.index,start.time,end.time,elapsed)
rm(breaks,breaks.time,duration,start.index,start.time,end.index,end.time,elapsed)

## FIND TEMPERATURE AND DEPTH
for (i in 1:nrow(depths)) {
     srt <- ceiling(depths$start.time[i])
     fsh <- floor(depths$end.time[i])
     depths$temp[i] <- mean(sen$tmp[srt:fsh], na.rm = TRUE)
     depths$pres[i] <- mean(dat2$p_Pa[depths$start.index[i]:depths$end.index[i]], na.rm = TRUE) - p_atm # Gage pressure in Pa
}
depths <- depths[order(depths$pres),]
g <- 9.81 # m/s^2, acceleration due to gravity
rho <- waterrho(depths$temp[1], unit = "C") # returns kg/m^3
delta_P <- depths$pres[1]
depths$depth[1] <- delta_P / (g*rho)
for (i in 2:nrow(depths)) {
     rho <- waterrho(depths$temp[i], unit = "C") # returns kg/m^3
     delta_P <- depths$pres[i] - depths$pres[i-1]
     depths$depth[i] <- depths$depth[i-1] + (delta_P / (g*rho))
}

## ADD average velocities
for (i in 1:nrow(depths)) {
     depths$ubar[i] <- mean(dat2$u[depths$start.index[i]:depths$end.index[i]], na.rm = TRUE)
     depths$vbar[i] <- mean(dat2$v[depths$start.index[i]:depths$end.index[i]], na.rm = TRUE)
     depths$wbar[i] <- mean(dat2$w[depths$start.index[i]:depths$end.index[i]], na.rm = TRUE)
}
depths$height <- (max(depths$depth) + 0.15) - depths$depth
depthuvw <- depths %>%
     select(height,ubar,vbar,wbar) %>%
     mutate(u=-wbar) %>%
     rename(v=vbar,w=ubar) %>%
     pivot_longer(cols = c("u","v","w"), names_to = "Direction", values_to = "Velocity")
# TO KEEP dat2 UP-TO-DATE
dat2 <- dat2 %>%
     rename(save=u) %>% # saving to keep from overwriting
     mutate(u=-w) %>%
     select(-w) %>% # remove before replacing
     rename(w=save)
# Plotting
depth.uvw <- ggplot(depthuvw) +
     geom_point(aes(x=Velocity, y=height, color=Direction)) +
     ylab("Height from bottom (m)") +
     xlab("Velocity (m/s)") +
     theme(legend.position="right", legend.key = element_blank(), 
           legend.text = element_text(face = "plain", size = 12),
           legend.title = element_text(face = "plain", size = 12)) + 
     theme(panel.background = element_rect(fill = "white", colour = "black")) + 
     theme(aspect.ratio = 1) +
     theme(axis.text = element_text(face = "plain", size = 12)) +
     theme(axis.title = element_text(face = "plain", size = 12))
ggsave("depth_v_uvw.eps", depth.uvw, device = "eps", dpi = 72)


