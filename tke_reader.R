# This code analyzes the data from the Nortek Vector ADV.

library(dplyr)

setwd("/Users/davidkahler/Documents/Hydrology_and_WRM/river_and_lake_mixing/ADV data/")
pck = read.table("TESTDE02.pck", header = FALSE, sep = "", dec = ".")
sen = read.table("TESTDE02.sen", header = FALSE, sep = "", dec = ".")
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
vhd = read.table("TESTDE02.vhd", header = FALSE, sep = "", dec = ".")
dat = read.table("TESTDE02.dat", header = FALSE, sep = "", dec = ".")
dat <- dat %>% rename(burst = V1, ensemble = V2, u = V3, v = V4, w = V5, amp1 = V6, amp2 = V7, amp3 = V8, snr1 = V9, snr2 = V10, snr3 = V11, corr1 = V12, corr2 = V13, corr3 = V14, p_dbar = V15, a1 = V16, a2 = V17, checksum = V18)
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
# 15   Pressure                         (dbar)
# 16   Analog input 1
# 17   Analog input 2
# 18   Checksum                         (1=failed)
sampling_rate = 64 # Hz, verify sampling rate in .hdr file under User setup

# Data check
max(dat$checksum)
hist(dat$snr1)
hist(dat$snr2)
hist(dat$snr3)

datetime <- array(-9999, dim = c(nrow(sen),2)) # col 1: days (day 1 is 01 Jan 2020), col 2: seconds of the day
u <- array(-9999, dim = c(nrow(sen),sampling_rate))
v <- array(-9999, dim = c(nrow(sen),sampling_rate))
w <- array(-9999, dim = c(nrow(sen),sampling_rate))
u_ave <- array(0, dim = c(nrow(sen),2))
v_ave <- array(0, dim = c(nrow(sen),2))
w_ave <- array(0, dim = c(nrow(sen),2))
u_prime <- u
v_prime <- v
w_prime <- w
for (i in 1:nrow(sen)) {
      datetime[i,1] <- ( as.numeric(as.Date(paste(sen$yea[i], sen$mon[i], sen$day[i], sep = " "), format = "%Y %m %d")) - as.numeric(as.Date("2019 12 31", format = "%Y %m %d")) )
      datetime[i,2] <- sen$sec[i] + (sen$mnt[i]-1)*60 + (sen$hou[1]-1)*3600
      for (j in 1:sampling_rate) {
            dat_index <- (sampling_rate*(i-1)) + j
            if (is.na(dat$u[dat_index])==FALSE) {
                  if (dat$checksum[dat_index]>0) {
                        print(paste0("error at ", dat$ensemble[dat_index]))
                  }
                  u[i,j] <- dat$u[dat_index]
                  u_ave[i,2] = u_ave[i,2] + 1
                  v[i,j] <- dat$v[dat_index]
                  v_ave[i,2] = v_ave[i,2] + 1
                  w[i,j] <- dat$w[dat_index]
                  w_ave[i,2] = w_ave[i,2] + 1
            }
      }
      u_ave[i,1] <- mean(u[i,])
      v_ave[i,1] <- mean(v[i,])
      w_ave[i,1] <- mean(w[i,])
      for (j in 1:sampling_rate) {
            if (u_ave[i,2]==sampling_rate) {
                  u_prime[i,j] <- u[i,j] - u_ave[i,1]
                  v_prime[i,j] <- v[i,j] - v_ave[i,1]
                  w_prime[i,j] <- w[i,j] - w_ave[i,1]
            }
      }
}

plot(u_ave[,1], ylim = c(-0.5, 0.5), type = "l")
plot(v_ave[,1], ylim = c(-0.5, 0.5), type = "l")
plot(w_ave[,1], ylim = c(-0.5, 0.5), type = "l")
