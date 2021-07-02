###########################################################################
###########################################################################
###########################################################################
##                      Simulaciones a 303K                              ##
###########################################################################
###########################################################################
###########################################################################

###################################
## Plot comparativo pocket  303K ##
###################################

library('tidyverse')

#Info replicas WT
r_WT_1_303 = read.table('303K/WT/1/mdpout_descriptors.txt', h=T)
r_WT_2_303 = read.table('303K/WT/2/mdpout_descriptors.txt', h=T)
r_WT_2_303['snapshot'] = r_WT_2_303['snapshot'] + 10001
r_WT_303 <- rbind(r_WT_1_303, r_WT_2_303)
r_WT_303['Time'] = seq(from = 0, to = 500.025, by = 0.025)
ylim =c(0, 1000)
plot(smooth.spline(r_WT_303[, 'Time'], r_WT_303[, 'pock_volume'], df=40), col='blue', lwd=3, ylim=ylim, ty='l', xlab='Time [ns]', ylab='volume (Å3)')
par(new=T)
#Info replicas S214H
r_S214H_1_303 = read.table('303K/S214H/1/mdpout_descriptors.txt', h=T)
r_S214H_2_303 = read.table('303K/S214H/2/mdpout_descriptors.txt', h=T)
r_S214H_2_303['snapshot'] = r_S214H_2_303['snapshot'] + 10001
r_S214H_303 <- rbind(r_S214H_1_303, r_S214H_2_303)
r_S214H_303['Time'] = seq(from = 0, to = 500.025, by = 0.025)
plot(smooth.spline(r_S214H_303[, 'Time'], r_S214H_303[, 'pock_volume'], df=40), col='yellow', lwd=3, ylim=ylim, ty='l', xlab='Time [ns]', ylab='volume (Å3)')
par(new=T)
#Info repolicas TSPETase
r_TSPETase_1_303 = read.table('303K/TSPETase/1/mdpout_descriptors.txt', h=T)
r_TSPETase_2_303 = read.table('303K/TSPETase/2/mdpout_descriptors.txt', h=T)
r_TSPETase_2_303['snapshot'] = r_TSPETase_2_303['snapshot'] + 10001
r_TSPETase_303 <- rbind(r_TSPETase_1_303, r_TSPETase_2_303)
r_TSPETase_303['Time'] = seq(from = 0, to = 500.025, by = 0.025)
plot(smooth.spline(r_TSPETase_303[, 'Time'], r_TSPETase_303[, 'pock_volume'], df=40), col='green', lwd=3, ylim=ylim, ty='l', xlab='Time [ns]', ylab='volume (Å3)')
par(new=T)
#Info replicas DuraPETase
r_DuraPETase_1_303 = read.table('303K/DuraPETase/1/mdpout_descriptors.txt', h=T)
r_DuraPETase_2_303 = read.table('303K/DuraPETase/2/mdpout_descriptors.txt', h=T)
r_DuraPETase_2_303['snapshot'] = r_DuraPETase_2_303['snapshot'] + 10001
r_DuraPETase_303 <- rbind(r_DuraPETase_1_303, r_DuraPETase_2_303)
r_DuraPETase_303['Time'] = seq(from = 0, to = 500.025, by = 0.025)
plot(smooth.spline(r_DuraPETase_303[, 'Time'], r_DuraPETase_303[, 'pock_volume'], df=40), col='red', lwd=3, ylim=ylim, ty='l', xlab='Time [ns]', ylab='volume (Å3)')
par(new=T)
abline(v = 250, lwd=3, lty = 2)
## Add Legend to the plot
legend(0, 1000, legend=c('WT', 'S214H', 'TSPETase', 'DuraPETase'), col=c('blue', 'yellow', 'green', 'red'),
       lty=1:1, cex=0.8)
## Add a title
title(main = 'Pocket volume variation during simulations at 303K', cex.main=1.0, col='black')

###################################
## Plot comparativo pocket  315K ##
###################################

library('tidyverse')

#Info replicas WT
r_WT_1_315 = read.table('315K/WT/1/mdpout_descriptors.txt', h=T)
r_WT_2_315 = read.table('315K/WT/2/mdpout_descriptors.txt', h=T)
r_WT_2_315['snapshot'] = r_WT_2_315['snapshot'] + 10001
r_WT_315 <- rbind(r_WT_1_315, r_WT_2_315)
r_WT_315['Time'] = seq(from = 0, to = 500.025, by = 0.025)
ylim =c(0, 1000)
plot(smooth.spline(r_WT_315[, 'Time'], r_WT_315[, 'pock_volume'], df=40), col='blue', lwd=3, ylim=ylim, ty='l', xlab='Time [ns]', ylab='volume (Å3)')
par(new=T)
#Info replicas S214H
r_S214H_1_315 = read.table('315K/S214H/1/mdpout_descriptors.txt', h=T)
r_S214H_2_315 = read.table('315K/S214H/2/mdpout_descriptors.txt', h=T)
r_S214H_2_315['snapshot'] = r_S214H_2_315['snapshot'] + 10001
r_S214H_315 <- rbind(r_S214H_1_315, r_S214H_2_315)
r_S214H_315['Time'] = seq(from = 0, to = 500.025, by = 0.025)
plot(smooth.spline(r_S214H_315[, 'Time'], r_S214H_315[, 'pock_volume'], df=40), col='yellow', lwd=3, ylim=ylim, ty='l', xlab='Time [ns]', ylab='volume (Å3)')
par(new=T)
#Info repolicas TSPETase
r_TSPETase_1_315 = read.table('315K/TSPETase/1/mdpout_descriptors.txt', h=T)
r_TSPETase_2_315 = read.table('315K/TSPETase/2/mdpout_descriptors.txt', h=T)
r_TSPETase_2_315['snapshot'] = r_TSPETase_2_303['snapshot'] + 10001
r_TSPETase_315 <- rbind(r_TSPETase_1_315, r_TSPETase_2_315)
r_TSPETase_315['Time'] = seq(from = 0, to = 500.025, by = 0.025)
plot(smooth.spline(r_TSPETase_315[, 'Time'], r_TSPETase_315[, 'pock_volume'], df=40), col='green', lwd=3, ylim=ylim, ty='l', xlab='Time [ns]', ylab='volume (Å3)')
par(new=T)
#Info replicas DuraPETase
r_DuraPETase_1_315 = read.table('315K/DuraPETase/1/mdpout_descriptors.txt', h=T)
r_DuraPETase_2_315 = read.table('315K/DuraPETase/2/mdpout_descriptors.txt', h=T)
r_DuraPETase_2_303['snapshot'] = r_DuraPETase_2_315['snapshot'] + 10001
r_DuraPETase_315 <- rbind(r_DuraPETase_1_315, r_DuraPETase_2_315)
r_DuraPETase_315['Time'] = seq(from = 0, to = 500.025, by = 0.025)
plot(smooth.spline(r_DuraPETase_315[, 'Time'], r_DuraPETase_315[, 'pock_volume'], df=40), col='red', lwd=3, ylim=ylim, ty='l', xlab='Time [ns]', ylab='volume (Å3)')
par(new=T)
abline(v = 250, lwd=3, lty = 2)
## Add Legend to the plot
legend(0, 1000, legend=c('WT', 'S214H', 'TSPETase', 'DuraPETase'), col=c('blue', 'yellow', 'green', 'red'),
       lty=1:1, cex=0.8)
## Add a title
title(main = 'Pocket volume variation during simulations at 315K', cex.main=1.0, col='black')


##Extract information of DuraPETase
which(grepl(11.88, r_DuraPETase_315$pock_volume))
max(r_DuraPETase_315[r_DuraPETase_315$pock_volume>0,"pock_volume"])
which(grepl(1812.84, r_DuraPETase_315$pock_volume))
