# This code analyzes the data from the Nortek Vector ADV.
# taken from tke_reader.r, currently working on this code.

# to install the needed packages, run the following lines:
#install.packages("dplyr")
#install.packages("ggplot2")
#install.packages("lubridate")
library(dplyr)
library(ggplot2)
library(lubridate)
library(latex2exp)

setwd("c:/Users/duquesne/Documents/nortek/data") # lab laptop
#setwd("/Users/davidkahler/Documents/Hydrology_and_WRM/river_and_lake_mixing/ADV_data/") # David's computer
fh <- "903MON17" # filename header
fn_sen <- paste(fh, "sen", sep = ".")
sen <- read.table(fn_sen, header = FALSE, sep = "", dec = ".")
sen <- sen %>% rename(mon = V1, day = V2, yea = V3, hou = V4, mnt = V5, sec = V6, err = V7, sta = V8, bat = V9, ssp = V10, hed = V11, pit = V12, rol = V13, tmp = V14, a1 = V15, checksum = V16)
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
d <- 24*3600*as.numeric(as.Date(paste(sen$yea[1], sen$mon[1], sen$day[1], sep="-"), origin="1970-01-01")) # number of seconds that gives the day
h <- sen$sec[1]+60*sen$mnt[1]+3600*sen$hou[1] # time in seconds
starttime <- as_datetime(d + h) # lubridate datetime for the start of the data

fn_dat <- paste(fh, "dat", sep = ".")
dat <- read.table(fn_dat, header = FALSE, sep = "", dec = ".")
# or
# x <- file.choose()
# dat <- read.table(x, header = FALSE, sep = "", dec = ".")
dat <- dat %>% rename(burst = V1, ensemble = V2, w = V3, u = V4, v = V5, amp1 = V6, amp2 = V7, amp3 = V8, snr1 = V9, snr2 = V10, snr3 = V11, corr1 = V12, corr2 = V13, corr3 = V14, p_dbar = V15, a1 = V16, a2 = V17, checksum = V18)#v3 changed from u to w, v4 was changed from v to u, v5 was changed from w to v so it can match how we placed it in the water
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
#pck = read.table("MON103.pck", header = FALSE, sep = "", dec = ".")
#vhd = read.table("MON103.vhd", header = FALSE, sep = "", dec = ".")
sampling_rate = 64 # Hz, verify sampling rate in .hdr file under User setup

#u_new and v_new loop
x <- atan(0.18447139/0.3180326)
for(i in 1:nrow(dat)) {
     dat$u_new[i] <- dat$u[i] * cos(x) + dat$v[i] * sin(x)
     dat$v_new[i] <- -dat$u[i] * sin(x) + dat$v[i] * cos(x)
} 

# Data check
records <- nrow(dat)
seconds <- nrow(sen)
top_samples <- seconds*sampling_rate # this is the number of records that there would be if every second had all of the sampling rate's elements filled in, the value (records) should not exceed this.
numbers <- ((records<=top_samples)&&(records>(top_samples-(2*sampling_rate)))) # we allow twice the values to account for missing values at the first second and the last second.
error_checksum <- max(dat$checksum) # If there are any errors recorded, this will be = 1
par(mfrow = c(3,1), mar = c(4,4,1,1)) # mfrow=c(nrows, ncols), https://www.statmethods.net/advgraphs/layout.html
hist(dat$snr1, breaks = c(-100,0,5,10,15,20,25,30,35,40,45,50,55,60,100), xlim = c(0,60), ylab = "Beam 1", xlab = "", main = "")
hist(dat$snr2, breaks = c(-100,0,5,10,15,20,25,30,35,40,45,50,55,60,100), xlim = c(0,60), ylab = "Beam 2", xlab = "", main = "")
hist(dat$snr3, breaks = c(-100,0,5,10,15,20,25,30,35,40,45,50,55,60,100), xlim = c(0,60), ylab = "Beam 3", xlab = "Signal-to-Noise Ratio (dB)", main = "")

# bar is the averaging window for U_bar and to compute the deviations from the mean, easiest to express 
# as a multiple of the sampling rate, therefore measured in seconds: bar_s.  Can take non-integer values.
bar_s <- 15 # averaging window in seconds.  Should evaluate range from depth/mean representative velocity to entire period
bar <- round(bar_s * sampling_rate) # in indexed values [i], round() needed to ensure that it fits within the dataset
# gives decimal time for each data record
dat$time <- starttime + (c(0:(nrow(dat)-1)))/sampling_rate

## Examine data:
par(mfrow = c(1,1))
plot(dat$time,(1e4*dat$p_dbar), type = "l",ylab = "Pressure (Pa)", xlab = "Time (s)")
lines(c(min(dat$time),max(dat$time)),c(101325,101325)) # Places a line at what should be the surface of the water

## Figure out an estimate of pressure:
# temperature <- array(NA, dim = nrow(dat))
# for (i in 1:nrow(dat)) {
#       if (dat$p_dbar>45000) {
#             temperature[i] <- # WE NEED TO PULL TEMP DATA FROM SEN...
#       }
# }
atmos <- mean(1e4*dat$p_dbar[1:10]) # Pa, to subtract atmospheric pressure
dat$depth <- -(1e4*dat$p_dbar - atmos)/(9.81*997)
#xlim = c(ymd_hms("2021-06-29 15:09:00",ymd_hms("2021-06-29 15:11:00")))
par(mfrow = c(1,1), mar = c(4,4,1,1))
plot(dat$time,dat$depth, type = "l", ylim = c(-6,0), ylab = "Depth (m)", xlab = "Time")

par(mfrow = c(3,1), mar = c(4,4,1,1))
plot(hms::as_hms(dat$time),dat$u, ylim = c(-10, 10), type = "l",ylab = "u (m/s)", xlab = "")
plot(hms::as_hms(dat$time),dat$v, ylim = c(-10, 10), type = "l",ylab = "v (m/s)", xlab = "")
plot(hms::as_hms(dat$time),dat$w, ylim = c(-10, 10), type = "l",ylab = "w (m/s)", xlab = "Time (s, from midnight)")

# ANALYSIS BY AVERAGING WINDOW:
u <- array(NA, dim = c(ceiling(nrow(dat)/bar),bar))
v <- u
w <- v
time <- array(-9999, dim = c(ceiling(nrow(dat)/bar))) # will identify the start of every averaged window.
u_ave <- array(0, dim = c(nrow(u),2))
v_ave <- u_ave
w_ave <- u_ave
u_prime <- u
v_prime <- v
w_prime <- w
uu <- array(0, dim = c(nrow(u)))
vv <- uu
ww <- uu
uv <- uu
uw <- uu
vw <- uu
for (i in 1:(nrow(u))) {
     time[i] <- dat$time[(bar*(i-1)) + 1] # this is the start time of the averaging window at time step, i
     for (j in 1:bar) { # this will cycle over each averaging window
          dat_index <- (bar*(i-1)) + j
          if (is.na(dat$u[dat_index])==FALSE) {
               if (dat$checksum[dat_index]>0) {
                    print(paste0("error at ", dat$ensemble[dat_index]))
               }
               u[i,j] <- dat$u[dat_index] # organize data
               u_ave[i,2] = u_ave[i,2] + 1
               v[i,j] <- dat$v[dat_index]
               v_ave[i,2] = v_ave[i,2] + 1
               w[i,j] <- dat$w[dat_index]
               w_ave[i,2] = w_ave[i,2] + 1
          }
     }
     u_ave[i,1] <- mean(u[i,], na.rm = TRUE)
     v_ave[i,1] <- mean(v[i,], na.rm = TRUE)
     w_ave[i,1] <- mean(w[i,], na.rm = TRUE)
     for (j in 1:bar) {
          u_prime[i,j] <- u[i,j] - u_ave[i,1]
          v_prime[i,j] <- v[i,j] - v_ave[i,1]
          w_prime[i,j] <- w[i,j] - w_ave[i,1]
     }
     uu[i] <- mean((u_prime[i,]^2))
     vv[i] <- mean((v_prime[i,]^2))
     ww[i] <- mean((w_prime[i,]^2))
     uv[i] <- mean((u_prime[i,]*v_prime[i,]))
     uw[i] <- mean((u_prime[i,]*w_prime[i,]))
     vw[i] <- mean((v_prime[i,]*w_prime[i,]))
}

par(mfrow = c(3,1), mar = c(4,4,1,1))
# REM: x-axis is datetime range
plot(hms::as_hms(as_datetime(as.numeric(time))),u_ave[,1], ylim = c(-1, 1), type = "l",ylab = "u (m/s)", xlab = "")
plot(hms::as_hms(as_datetime(as.numeric(time))),v_ave[,1], ylim = c(-1, 1), type = "l",ylab = "v (m/s)", xlab = "")
plot(hms::as_hms(as_datetime(as.numeric(time))),w_ave[,1], ylim = c(-1, 1), type = "l",ylab = "w (m/s)", xlab = "Time (s)")

plot(hms::as_hms(as_datetime(as.numeric(time))),uu, ylim = c(-1, 1), type = "l",ylab = "uu (m/s)", xlab = "")
plot(hms::as_hms(as_datetime(as.numeric(time))),vv, ylim = c(-1, 1), type = "l",ylab = "vv (m/s)", xlab = "")
plot(hms::as_hms(as_datetime(as.numeric(time))),ww, ylim = c(-1, 1), type = "l",ylab = "ww (m/s)", xlab = "Time (s)")

plot(hms::as_hms(as_datetime(as.numeric(time))),uv, ylim = c(-1, 1), type = "l",ylab = "uv (m/s)", xlab = "")
plot(hms::as_hms(as_datetime(as.numeric(time))),vw, ylim = c(-1, 1), type = "l",ylab = "vw (m/s)", xlab = "")
plot(hms::as_hms(as_datetime(as.numeric(time))),vw, ylim = c(-1, 1), type = "l",ylab = "vw (m/s)", xlab = "Time (s)")

# to zoom in on an area of interest, find indices:
start <- as.numeric(ymd_hms("2021-09-03 19:18:00")) # Enter start time here as "YYYY-MM-DD HH:MM:SS" in 24-hour time
end <- as.numeric(ymd_hms("2021-09-03 19:36:00")) # End time, same format
diff_s <- abs(time-start) # finds the difference between the time entries and start time (in case we don't hit it exactly)
diff_e <- abs(time-end)
m_s <- 10*bar_s # allocate variable for minimum finding.  larger than what will be found
m_e <- m_s
s <- NA
e <- NA
for (i in 1:(length(time))) {
     if (diff_s[i] < m_s) {
          m_s <- diff_s[i]
          s <- i # index of start
     }
     if (diff_e[i] <= m_e) {
          m_e <- diff_e[i]
          e <- i # index of end
     }
}

par(mfrow = c(3,1), mar = c(4,4,1,1))
plot(hms::as_hms(as_datetime(as.numeric(time[s:e]))), uu[s:e], ylim = c(0, 0.3), type = "l", ylab = "uu", xlab = "")
plot(hms::as_hms(as_datetime(as.numeric(time[s:e]))), vv[s:e], ylim = c(0, 0.4), type = "l", ylab = "vv", xlab = "")
plot(hms::as_hms(as_datetime(as.numeric(time[s:e]))), ww[s:e], ylim = c(0, 0.6), type = "l", ylab = "ww", xlab = "Time (s)")

par(mfrow = c(3,1), mar = c(4,4,1,1))
plot(uv, ylim = c(-0.3, 0.3), type = "l", ylab = "uv", xlab = "")
plot(uw, ylim = c(-0.3, 0.3), type = "l", ylab = "uw", xlab = "")
plot(vw, ylim = c(-0.3, 0.3), type = "l", ylab = "vw", xlab = "Time (s)")

# ANALYSIS AT EACH DEPTH
# Check 16:57 to 16:59
start <- as.numeric(ymd_hms("2021-09-03 19:10:00")) # Enter start time here as "YYYY-MM-DD HH:MM:SS" in 24-hour time
end <- as.numeric(ymd_hms("2021-09-03 19:35:30")) # End time, same format
for (i in 1:nrow(dat)) {
     if (as.numeric(dat$time[i]) < start) {
          s <- i
     }
     if (as.numeric(dat$time[i]) < end) {
          e <- i
     }
}
s <- s + 1
e <- e + 1
u_ave <- mean(dat$u[s:e])
v_ave <- mean(dat$v[s:e])
w_ave <- mean(dat$w[s:e])
ui <- dat$u[s:e] - u_ave
vi <- dat$v[s:e] - v_ave
wi <- dat$w[s:e] - w_ave
uiuj <- array(NA, dim = c(3,3))
uiuj[1,1] <- mean(ui^2)
uiuj[1,2] <- mean(ui*vi)
uiuj[1,3] <- mean(ui*wi)
uiuj[2,2] <- mean(vi^2)
uiuj[2,3] <- mean(vi*wi)
uiuj[3,3] <- mean(wi^2)
U <- sqrt((u_ave^2) + (v_ave^2) + (w_ave^2)) #average velocity in all directions for time window
avgdepth <- mean(dat$depth[s:e]) #average depth for 2 minute time span

#create dataframes for ggplots
df <- data.frame(ui, vi, wi)

#plots at each depth
ggplot(df, aes(x= wi, y= vi)) +
     geom_point() +
     xlim(-3, 3) +
     ylim(-6, 6) +
     xlab(TeX('$w_i$')) +
     ylab(TeX('$v_i$'))+
     theme(panel.background = element_rect(fill = "white", colour = "black")) +
     theme(aspect.ratio = 1)+
     theme(axis.text = element_text(face = "plain", size = 12))

# depth and dissolved oxygen for Monongahela and Allegheny Rivers

#629Mon
#cable crossing
file <- file.choose()
cable_crossing_DO <- read.csv(file)
cable_crossing_DO <- rename(cable_crossing_DO, time = ï..time) # This appears to be needed for Windows... 

ggplot(cable_crossing_DO, aes(x = DO.mg.l. , y= Depth.m.))+
     geom_point() +
     labs(x = "DO (mg/l)", y = "Depth (m)") +
     ggtitle("DO profile in Monongahela River") +
     theme(panel.background = element_rect(fill = "white", colour = "black")) +
     theme(aspect.ratio = 1) +
     theme(axis.text = element_text(face = "plain", size = 12))

#birmingham bridge
file <- file.choose()
birmingham_bridge_DO <- read.csv(file)
birmingham_bridge_DO <- rename(birmingham_bridge_DO, time =ï..time)

ggplot(birmingham_bridge_DO, aes(x = DO.mg.l. , y= Depth.m.))+
     geom_point() +
     labs(x = "DO (mg/l)", y = "Depth (m)") +
     ggtitle("DO profile in Monongahela River") +
     theme(panel.background = element_rect(fill = "white", colour = "black")) +
     theme(aspect.ratio = 1) +
     theme(axis.text = element_text(face = "plain", size = 12))

#528AR
#need to pull depth from 528AR13 before this one is finished 
file <- file.choose()
DO528AR13 <- read.csv(file)
DO528AR13 <- rename(DO528AR13, time = ï..time) # This appears to be needed for Windows... 

ggplot(DO528AR13, aes(x = DO.mg.l. , y= Depth.m.))+
     geom_point() +
     labs(x = "DO (mg/l)", y = "Depth (m)") +
     ggtitle("DO profile in Allegheny River") +
     theme(panel.background = element_rect(fill = "white", colour = "black")) +
     theme(aspect.ratio = 1) +
     theme(axis.text = element_text(face = "plain", size = 12))

#528AR214
file <- file.choose()
DO528AR214 <- read.csv(file)
DO528AR214 <- rename(DO528AR214, time = ï..time) # This appears to be needed for Windows... 

ggplot(DO528AR214, aes(x = DO.mg.l. , y= Depth.m.))+
     geom_point() +
     labs(x = "DO (mg/l)", y = "Depth (m)") +
     ggtitle("DO profile in Allegheny River") +
     theme(panel.background = element_rect(fill = "white", colour = "black")) +
     theme(aspect.ratio = 1) +
     theme(axis.text = element_text(face = "plain", size = 12))


# Spectra
par(mfrow = c(3,1), mar = c(4,4,2,2))
hist(uv[start:stop,1], breaks = c(-10000,-1.5,-0.1,-0.09,-0.08,-0.07,-0.06,-0.05,-0.04,-0.03,-0.02,-0.01,0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,1.5,10000), xlim = c(-1.5,1.5), ylab = "u'v'", xlab = "", main = "")
hist(uw[start:stop,1], breaks = c(-10000,-1.5,-0.1,-0.09,-0.08,-0.07,-0.06,-0.05,-0.04,-0.03,-0.02,-0.01,0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,1.5,10000), xlim = c(-1.5,1.5), ylab = "u'w'", xlab = "", main = "")
hist(vw[start:stop,1], breaks = c(-10000,-1.5,-0.1,-0.09,-0.08,-0.07,-0.06,-0.05,-0.04,-0.03,-0.02,-0.01,0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,1.5,10000), xlim = c(-1.5,1.5), ylab = "v'w'", xlab = "", main = "")