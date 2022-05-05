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
     
     # Plot
     up <- ggplot(uvec) +
          geom_line(aes(x=dt,y=uprime)) +
          ylab("u' (m/s)") +
          xlab(element_blank()) +
          theme(panel.background = element_rect(fill = "white", colour = "black")) + 
          theme(axis.text = element_text(face = "plain", size = 12))
     vp <- ggplot(uvec) +
          geom_line(aes(x=dt,y=vprime)) +
          ylab("v' (m/s)") +
          xlab(element_blank()) +
          theme(panel.background = element_rect(fill = "white", colour = "black")) + 
          theme(axis.text = element_text(face = "plain", size = 12))
     wp <- ggplot(uvec) +
          geom_line(aes(x=dt,y=wprime)) +
          ylab("w' (m/s)") +
          xlab("Time (s)") +
          theme(panel.background = element_rect(fill = "white", colour = "black")) + 
          theme(axis.text = element_text(face = "plain", size = 12))
     gup <- ggplotGrob(up)
     gvp <- ggplotGrob(vp)
     gwp <- ggplotGrob(wp)
     setEPS() # https://www.geeksforgeeks.org/export-plot-to-eps-file-in-r/
     postscript(paste0("uvw_prime.i",as.character(i),".h",as.character(round(h, digits = 4)),".eps"))
     grid::grid.newpage()
     grid::grid.draw(rbind(gup,gvp,gwp))
     #ggsave(paste0("uvw_prime.i",as.character(i),".h",as.character(round(h, digits = 4))),".eps", uvwprimes, device = "eps", dpi = 72)
     dev.off()
     
     uvec$upwp_bar <- uvec$uprime*uvec$wprime
     upwp_bar <- mean(uprime*wprime, na.rm = TRUE)
     annotation <- data.frame(
          x = -0.2,
          y = -0.5,
          label = paste0("u'w'=",upwp_bar)
     )
     #depths$upwp_bar[i] <- upwp_bar # removed due to parallel, replacing at end of parallel loop
     RS <- ggplot(uvec) +
          geom_point(aes(x=uprime,y=wprime)) +
          xlim(c(-0.5,0.5)) +
          xlab("u'") +
          ylim(c(-0.5,0.5)) +
          ylab("w'") +
          geom_text(data=annotation, aes( x=x, y=y, label=label), 
                    color="blue") +
          theme(panel.background = element_rect(fill = "white", colour = "black")) + 
          theme(aspect.ratio = 1) +
          theme(axis.text = element_text(face = "plain", size = 12)) +
          theme(axis.title = element_text(face = "plain", size = 12))
     ggsave(paste0("uw_prime.i",as.character(i),".h",as.character(round(h, digits = 4)),".eps"), RS, device = "eps")
     
     print(c(i,upwp_bar)) # outputs to parallel output, rst, for Reynolds stress tensor
}

rst <- rst[order(rst[,1]),]
depths$upwp <- rst[,2]







