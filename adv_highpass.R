# Filter ADV data to remove oscillations due to instrument motion

library(dplyr)
library(ggplot2)
library(doParallel)

# Run adv_prep.R first.  Uses dat2 and depths, specifically

registerDoParallel(detectCores())
rst <- foreach (i = 1:nrow(depths), .combine = rbind) %dopar% {
     # Determine u', v', w'
     h <- depths$height[i]
     d <- depths$depth[i]
     dt <- dat2$dt[depths$start.index[i]:depths$end.index[i]]
     u <- dat2$u[depths$start.index[i]:depths$end.index[i]]
     ubar <- mean(u)
     uprime <- u-ubar
     v <- dat2$v[depths$start.index[i]:depths$end.index[i]]
     vbar <- mean(v)
     vprime <- v-vbar
     w <- dat2$w[depths$start.index[i]:depths$end.index[i]]
     wbar <- mean(w)
     wprime <- w-wbar
     uvec <- data.frame(dt,u,ubar,uprime,v,vbar,vprime,w,wbar,wprime)
     rm(dt,u,ubar,uprime,v,vbar,vprime,w,wbar,wprime)
     
     # Plot u', v', and w'
     # up <- ggplot(uvec) +
     #      geom_line(aes(x=dt,y=uprime)) +
     #      ylab("u' (m/s)") +
     #      xlab(element_blank()) +
     #      theme(panel.background = element_rect(fill = "white", colour = "black")) + 
     #      theme(axis.text = element_text(face = "plain", size = 12))
     # vp <- ggplot(uvec) +
     #      geom_line(aes(x=dt,y=vprime)) +
     #      ylab("v' (m/s)") +
     #      xlab(element_blank()) +
     #      theme(panel.background = element_rect(fill = "white", colour = "black")) + 
     #      theme(axis.text = element_text(face = "plain", size = 12))
     # wp <- ggplot(uvec) +
     #      geom_line(aes(x=dt,y=wprime)) +
     #      ylab("w' (m/s)") +
     #      xlab("Time (s)") +
     #      theme(panel.background = element_rect(fill = "white", colour = "black")) + 
     #      theme(axis.text = element_text(face = "plain", size = 12))
     # gup <- ggplotGrob(up)
     # gvp <- ggplotGrob(vp)
     # gwp <- ggplotGrob(wp)
     # setEPS() # https://www.geeksforgeeks.org/export-plot-to-eps-file-in-r/
     # postscript(paste0("uvw_prime.i",as.character(i),".h",as.character(round(h, digits = 4)),".eps"))
     # grid::grid.newpage()
     # grid::grid.draw(rbind(gup,gvp,gwp))
     # dev.off()
     
     # Filter - just filtering u' and w'
     N <- nrow(uvec)
     uvec$index <- c(0:(N-1))
     uvec$uprime.fft <- fft(uvec$uprime)
     uvec$wprime.fft <- fft(uvec$wprime)
     u.filter.dist <- 100 # found by visual inspection
     w.filter.dist <- 150
     
     # show filter region
     uprime.fft <- ggplot(uvec) +
          geom_line(aes(x=index, y=Mod(uprime.fft))) +
          geom_vline(xintercept = u.filter.dist, color = "red") +
          geom_vline(xintercept = (N-u.filter.dist-1), color = "red") +
          xlab(element_blank()) +
          ylab("u' frequency domain") +
          theme(panel.background = element_rect(fill = "white", colour = "black")) +
          theme(axis.text = element_text(face = "plain", size = 12)) +
          theme(axis.title = element_text(face = "plain", size = 12))
     wprime.fft <- ggplot(uvec) +
          geom_line(aes(x=index, y=Mod(wprime.fft))) +
          geom_vline(xintercept = w.filter.dist, color = "red") +
          geom_vline(xintercept = (N-w.filter.dist-1), color = "red") +
          xlab(element_blank()) +
          ylab("w' frequency domain") +
          theme(panel.background = element_rect(fill = "white", colour = "black")) + 
          theme(axis.text = element_text(face = "plain", size = 12)) +
          theme(axis.title = element_text(face = "plain", size = 12))
     uprime.fft.g <- ggplotGrob(uprime.fft)
     wprime.fft.g <- ggplotGrob(wprime.fft)
     setEPS() # https://www.geeksforgeeks.org/export-plot-to-eps-file-in-r/
     postscript(paste0("uw_fft.i",as.character(i),".h",as.character(round(h, digits = 4)),".eps"))
     grid::grid.newpage()
     grid::grid.draw(rbind(uprime.fft.g,wprime.fft.g))
     dev.off()
     
     # Impose filter
     u.filter <- rep(1, N) # creates a vector of ones
     w.filter <- rep(1, N)
     u.filter[1:u.filter.dist] = 0 # changes the start and end i.filter values to zero
     u.filter[(N-u.filter.dist):N] = 0 # end values
     w.filter[1:w.filter.dist] = 0 # changes the start and end i.filter values to zero
     w.filter[(N-w.filter.dist):N] = 0 # end values
     uvec$u.filtered = u.filter * uvec$uprime.fft # reconstructs
     uvec$w.filtered = w.filter * uvec$wprime.fft # reconstructs
     uvec$u.ifft <- Re(fft(uvec$u.filtered, inverse=TRUE) / N)
     uvec$w.ifft <- Re(fft(uvec$w.filtered, inverse=TRUE) / N)
     
     # plot time series
     up <- ggplot(uvec) +
          geom_line(aes(x=dt,y=uprime)) +
          geom_line(aes(x=dt,y=u.ifft),color="blue") +
          ylab("u' (m/s)") +
          xlab(element_blank()) +
          theme(panel.background = element_rect(fill = "white", colour = "black")) + 
          theme(axis.text = element_text(face = "plain", size = 12))
     wp <- ggplot(uvec) +
          geom_line(aes(x=dt,y=wprime)) +
          geom_line(aes(x=dt,y=w.ifft),color="blue") +
          ylab("w' (m/s)") +
          xlab("Time (s)") +
          theme(panel.background = element_rect(fill = "white", colour = "black")) + 
          theme(axis.text = element_text(face = "plain", size = 12))
     gup <- ggplotGrob(up)
     gwp <- ggplotGrob(wp)
     setEPS() # https://www.geeksforgeeks.org/export-plot-to-eps-file-in-r/
     postscript(paste0("uw_prime.ts.i",as.character(i),".h",as.character(round(h, digits = 4)),".eps"))
     grid::grid.newpage()
     grid::grid.draw(rbind(gup,gwp))
     dev.off()
     
     # Calculate u'w' averaged over the ensemble 
     uvec$upwp <- uvec$uprime*uvec$wprime
     upwp_bar <- mean(uvec$upwp, na.rm = TRUE)
     uvec$upwp_filtered <- uvec$u.ifft*uvec$w.ifft
     upwp_filtered <- mean(uvec$upwp_filtered, na.rm = TRUE)
     prime <- data.frame(
          x = -0.20,
          y = -0.45,
          label = paste0("u'w'=",upwp_bar)
     )
     filtered <- data.frame(
          x = -0.20,
          y = -0.50,
          label = paste0("u'w'=",upwp_filtered," (filtered)")
     )
     
     RS <- ggplot(uvec) +
          geom_point(aes(x=uprime,y=wprime)) +
          geom_point(aes(x=u.ifft,y=w.ifft),linetype="dashed",color="blue") +
          xlim(c(-0.5,0.5)) +
          xlab("u'") +
          ylim(c(-0.5,0.5)) +
          ylab("w'") +
          geom_text(data=prime, aes(x=x, y=y, label=label)) +
          geom_text(data=filtered, aes(x=x, y=y, label=label), color="blue") +
          theme(panel.background = element_rect(fill = "white", colour = "black")) + 
          theme(aspect.ratio = 1) +
          theme(axis.text = element_text(face = "plain", size = 12)) +
          theme(axis.title = element_text(face = "plain", size = 12))
     ggsave(paste0("uw_prime.i",as.character(i),".h",as.character(round(h, digits = 4)),".eps"), RS, device = "eps")
     
     print(c(i,upwp_bar,upwp_filtered)) # outputs to parallel output, rst, for Reynolds stress tensor
}

rst <- rst[order(rst[,1]),]
depths$upwp <- rst[,2]
depths$upwp_filtered <- rst[,3]






