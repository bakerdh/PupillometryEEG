library(signal)

outputplot <- 0

alltargets <- array(0,dim=c(2,2,60))
allmasks <- alltargets
allwaveforms <- array(0,dim=c(2,2,60,1200))
allspectra <- allwaveforms
timeseq <- seq(1/120,10,length.out=120*10)
targetindex <- (2*10)+1
maskindex <- (1.6*10)+1

for (block in 1:2){
infofile <- read.csv(paste('~/Google Drive/Current work/Research/Melatonin/DHBcrf/Block',block,'/info.csv',sep=''))
pupilstarttime <- as.numeric(as.character(infofile[4,2]))
localstarttime <- as.numeric(as.character(infofile[5,2]))

if (block==1){psychoutput <- read.csv('~/Google Drive/Current work/Research/Melatonin/DHBcrf/Block1/DHB_CRFstudy_2019_Aug_23_1626.csv')}
if (block==2){psychoutput <- read.csv('~/Google Drive/Current work/Research/Melatonin/DHBcrf/Block2/DHB_CRFstudy_2019_Aug_23_1643.csv')}

trialtimes <- psychoutput$trialonset - pupilstarttime
condorder <- psychoutput$condition

pdata <- read.csv(paste('~/Google Drive/Current work/Research/Melatonin/DHBcrf/Block',block,'/pupil_positions.csv',sep=''))
pdata2 <- pdata[,c(1,3,4,14)]
pdata2[,1] <- pdata2[,1] - localstarttime
pdata2 <- pdata2[which(pdata2[,3]>0.2),]

for (eye in 1:2){
eyedata <- pdata2[which(pdata2[,2]==(2-eye)),]
for (cond in 1:60){
a <- which(eyedata[,1]>trialtimes[cond]+2)
b <- which(eyedata[,1]<(trialtimes[cond]+12))
i <- intersect(a,b)
trial <- eyedata[i,]
trial[,1] <- trial[,1] - trial[1,1]
if (nrow(trial)>3){
resampled <- interp1(trial[,1],trial[,4],timeseq,method='linear',extrap=TRUE)
fspec <- (fft(resampled)/length(resampled))
alltargets[block,eye,condorder[cond]] <- fspec[targetindex]
allmasks[block,eye,condorder[cond]] <- fspec[maskindex]
allwaveforms[block,eye,condorder[cond],] <- resampled - mean(resampled)
allspectra[block,eye,condorder[cond],] <- fspec
}}}
}

contrastsdB <- 20*log10(c(6,12,24,48,96))

colvect <- c('red','blue','darkgreen','grey','purple','orange')
plotlims <- c(12,40,0,0.04)  # define the x and y limits of the plot (minx,maxx,miny,maxy)
ticklocsx <- seq(12,40,6)    # locations of tick marks on x axis
ticklocsy <- seq(0,0.04,0.01)    # locations of tick marks on y axis
ticklabelsx <-ticklocsx        # set labels for x ticks
ticklabelsy <- ticklocsy    # set labels for y ticks

if(outputplot==2){pdf("CRF1.pdf", bg="transparent", height = 5.5, width = 5.5)}

par(pty="s")  # make axis square
plot(x=NULL,y=NULL,axes=FALSE, ann=FALSE, xlim=plotlims[1:2], ylim=plotlims[3:4])   # create an empty axis of the correct dimensions
axis(1, at=ticklocsx, tck=0.01, lab=F, lwd=2)     # plot tick marks (no labels)
axis(2, at=ticklocsy, tck=0.01, lab=F, lwd=2)
axis(3, at=ticklocsx, tck=0.01, lab=F, lwd=2)
axis(4, at=ticklocsy, tck=0.01, lab=F, lwd=2)
mtext(text = ticklabelsx, side = 1, at=ticklocsx)     # add the tick labels
mtext(text = ticklabelsy, side = 2, at=ticklocsy, line=0.2, las=1)  # the 'line' command moves away from the axis, the 'las' command rotates to vertical
box(lwd=2)      # draw a box around the graph
title(xlab="Target contrast (dB)", col.lab=rgb(0,0,0), line=1.2, cex.lab=1.5)    # titles for axes
title(ylab="Amplitude (mm) at 2Hz", col.lab=rgb(0,0,0), line=1.5, cex.lab=1.5)

meanamps <- apply(alltargets,3,mean)
for (cond in 1:3){
  temp <- meanamps[((10*(cond-1))+1):(10*cond)]
  toplot <- abs(temp[1:5] + temp[6:10])/2
  lines(contrastsdB,toplot, col=colvect[cond], lwd=3, cex=0.5)     # draw a line connecting the points
  points(contrastsdB,toplot, pch = 21, col='black', bg=colvect[cond], cex=1.6, lwd=3)   # draw the data points themselves
}  
legend(12,0.04,c('Mon','Bin','Dich'),pch=21,pt.bg=colvect[1:3],pt.lwd=3,pt.cex=1.6,box.lwd=2)
if(outputplot>0){dev.off()}  # this line goes after you've finished plotting (to output the example below, move it to the bottom of the script)


if(outputplot==2){pdf("CRF2.pdf", bg="transparent", height = 5.5, width = 5.5)}

par(pty="s")  # make axis square
plot(x=NULL,y=NULL,axes=FALSE, ann=FALSE, xlim=plotlims[1:2], ylim=plotlims[3:4])   # create an empty axis of the correct dimensions
axis(1, at=ticklocsx, tck=0.01, lab=F, lwd=2)     # plot tick marks (no labels)
axis(2, at=ticklocsy, tck=0.01, lab=F, lwd=2)
axis(3, at=ticklocsx, tck=0.01, lab=F, lwd=2)
axis(4, at=ticklocsy, tck=0.01, lab=F, lwd=2)
mtext(text = ticklabelsx, side = 1, at=ticklocsx)     # add the tick labels
mtext(text = ticklabelsy, side = 2, at=ticklocsy, line=0.2, las=1)  # the 'line' command moves away from the axis, the 'las' command rotates to vertical
box(lwd=2)      # draw a box around the graph
title(xlab="Target contrast (dB)", col.lab=rgb(0,0,0), line=1.2, cex.lab=1.5)    # titles for axes
title(ylab="Amplitude (mm) at 2Hz", col.lab=rgb(0,0,0), line=1.5, cex.lab=1.5)

meanamps <- apply(alltargets,3,mean)
condlist <- c(1,5,6)
for (cond in 1:3){
  temp <- meanamps[((10*(condlist[cond]-1))+1):(10*condlist[cond])]
  toplot <- abs(temp[1:5] + temp[6:10])/2
  lines(contrastsdB,toplot, col=colvect[condlist[cond]], lwd=3, cex=0.5)     # draw a line connecting the points
  points(contrastsdB,toplot, pch = 21, col='black', bg=colvect[condlist[cond]], cex=1.6, lwd=3)   # draw the data points themselves
}  
legend(12,0.04,c('Mon','Bin X','Dich X'),pch=21,pt.bg=colvect[condlist],pt.lwd=3,pt.cex=1.6,box.lwd=2)
if(outputplot>0){dev.off()}  # this line goes after you've finished plotting (to output the example below, move it to the bottom of the script)



if(outputplot==2){pdf("CRF3.pdf", bg="transparent", height = 5.5, width = 5.5)}

par(pty="s")  # make axis square
plot(x=NULL,y=NULL,axes=FALSE, ann=FALSE, xlim=plotlims[1:2], ylim=plotlims[3:4])   # create an empty axis of the correct dimensions
axis(1, at=ticklocsx, tck=0.01, lab=F, lwd=2)     # plot tick marks (no labels)
axis(2, at=ticklocsy, tck=0.01, lab=F, lwd=2)
axis(3, at=ticklocsx, tck=0.01, lab=F, lwd=2)
axis(4, at=ticklocsy, tck=0.01, lab=F, lwd=2)
mtext(text = ticklabelsx, side = 1, at=ticklocsx)     # add the tick labels
mtext(text = ticklabelsy, side = 2, at=ticklocsy, line=0.2, las=1)  # the 'line' command moves away from the axis, the 'las' command rotates to vertical
box(lwd=2)      # draw a box around the graph
title(xlab="Target contrast (dB)", col.lab=rgb(0,0,0), line=1.2, cex.lab=1.5)    # titles for axes
title(ylab="Amplitude (mm) at 1.6Hz", col.lab=rgb(0,0,0), line=1.5, cex.lab=1.5)

meanamps <- apply(allmasks,3,mean)
for (cond in 4:6){
  temp <- meanamps[((10*(cond-1))+1):(10*cond)]
  toplot <- abs(temp[1:5] + temp[6:10])/2
lines(contrastsdB,toplot, col=colvect[cond], lwd=3, cex=0.5)     # draw a line connecting the points
points(contrastsdB,toplot, pch = 21, col='black', bg=colvect[cond], cex=1.6, lwd=3)   # draw the data points themselves
}
legend(12,0.04,c('Mon X','Bin X','Dich X'),pch=21,pt.bg=colvect[4:6],pt.lwd=3,pt.cex=1.6,box.lwd=2)
if(outputplot>0){dev.off()}  # this line goes after you've finished plotting (to output the example below, move it to the bottom of the script)



if(outputplot==2){pdf("consensual1.pdf", bg="transparent", height = 5.5, width = 5.5)}

par(pty="s")  # make axis square
plot(x=NULL,y=NULL,axes=FALSE, ann=FALSE, xlim=plotlims[1:2], ylim=plotlims[3:4])   # create an empty axis of the correct dimensions
axis(1, at=ticklocsx, tck=0.01, lab=F, lwd=2)     # plot tick marks (no labels)
axis(2, at=ticklocsy, tck=0.01, lab=F, lwd=2)
axis(3, at=ticklocsx, tck=0.01, lab=F, lwd=2)
axis(4, at=ticklocsy, tck=0.01, lab=F, lwd=2)
mtext(text = ticklabelsx, side = 1, at=ticklocsx)     # add the tick labels
mtext(text = ticklabelsy, side = 2, at=ticklocsy, line=0.2, las=1)  # the 'line' command moves away from the axis, the 'las' command rotates to vertical
box(lwd=2)      # draw a box around the graph
title(xlab="Target contrast (dB)", col.lab=rgb(0,0,0), line=1.2, cex.lab=1.5)    # titles for axes
title(ylab="Amplitude (mm)", col.lab=rgb(0,0,0), line=1.5, cex.lab=1.5)

meanamps <- apply(alltargets,c(2,3),mean)

  toplot <- abs(meanamps[1,1:5] + meanamps[2,6:10])/2
  lines(contrastsdB,toplot, col='pink', lwd=3, cex=0.5)     # draw a line connecting the points
  points(contrastsdB,toplot, pch = 21, col='black', bg='pink', cex=1.6, lwd=3)   # draw the data points themselves

  toplot <- abs(meanamps[2,1:5] + meanamps[1,6:10])/2
  lines(contrastsdB,toplot, col='black', lwd=3, cex=0.5)     # draw a line connecting the points
  points(contrastsdB,toplot, pch = 21, col='black', bg='black', cex=1.6, lwd=3)   # draw the data points themselves
  
  meanamps <- apply(allmasks,c(2,3),mean)
  toplot <- abs(meanamps[1,31:35] + meanamps[2,36:40])/2
  lines(contrastsdB,toplot, col='pink', lwd=3, cex=0.5)     # draw a line connecting the points
  points(contrastsdB,toplot, pch = 22, col='black', bg='pink', cex=1.6, lwd=3)   # draw the data points themselves
  
  toplot <- abs(meanamps[2,31:35] + meanamps[1,36:40])/2
  lines(contrastsdB,toplot, col='black', lwd=3, cex=0.5)     # draw a line connecting the points
  points(contrastsdB,toplot, pch = 22, col='black', bg='black', cex=1.6, lwd=3)   # draw the data points themselves
  
legend(12,0.04,c('Ipsilateral 2Hz','Consensual 2Hz','Ipsilateral 1.6Hz','Consensual 1.6Hz'),pch=c(21,21,22,22),pt.bg=c('pink','black'),pt.lwd=3,pt.cex=1.6,box.lwd=2)
if(outputplot>0){dev.off()}  # this line goes after you've finished plotting (to output the example below, move it to the bottom of the script)



if(outputplot==2){pdf("consensual2.pdf", bg="transparent", height = 5.5, width = 5.5)}

par(pty="s")  # make axis square
plot(x=NULL,y=NULL,axes=FALSE, ann=FALSE, xlim=plotlims[1:2], ylim=plotlims[3:4])   # create an empty axis of the correct dimensions
axis(1, at=ticklocsx, tck=0.01, lab=F, lwd=2)     # plot tick marks (no labels)
axis(2, at=ticklocsy, tck=0.01, lab=F, lwd=2)
axis(3, at=ticklocsx, tck=0.01, lab=F, lwd=2)
axis(4, at=ticklocsy, tck=0.01, lab=F, lwd=2)
mtext(text = ticklabelsx, side = 1, at=ticklocsx)     # add the tick labels
mtext(text = ticklabelsy, side = 2, at=ticklocsy, line=0.2, las=1)  # the 'line' command moves away from the axis, the 'las' command rotates to vertical
box(lwd=2)      # draw a box around the graph
title(xlab="Target contrast (dB)", col.lab=rgb(0,0,0), line=1.2, cex.lab=1.5)    # titles for axes
title(ylab="Amplitude (mm)", col.lab=rgb(0,0,0), line=1.5, cex.lab=1.5)

meanamps <- apply(alltargets,c(2,3),mean)

toplot <- abs(meanamps[1,51:55] + meanamps[2,56:60])/2
lines(contrastsdB,toplot, col='pink', lwd=3, cex=0.5)     # draw a line connecting the points
points(contrastsdB,toplot, pch = 21, col='black', bg='pink', cex=1.6, lwd=3)   # draw the data points themselves

toplot <- abs(meanamps[2,51:55] + meanamps[1,56:60])/2
lines(contrastsdB,toplot, col='black', lwd=3, cex=0.5)     # draw a line connecting the points
points(contrastsdB,toplot, pch = 21, col='black', bg='black', cex=1.6, lwd=3)   # draw the data points themselves

meanamps <- apply(allmasks,c(2,3),mean)
toplot <- abs(meanamps[1,51:55] + meanamps[2,56:60])/2
lines(contrastsdB,toplot, col='pink', lwd=3, cex=0.5)     # draw a line connecting the points
points(contrastsdB,toplot, pch = 22, col='black', bg='pink', cex=1.6, lwd=3)   # draw the data points themselves

toplot <- abs(meanamps[2,51:55] + meanamps[1,56:60])/2
lines(contrastsdB,toplot, col='black', lwd=3, cex=0.5)     # draw a line connecting the points
points(contrastsdB,toplot, pch = 22, col='black', bg='black', cex=1.6, lwd=3)   # draw the data points themselves

legend(28,0.04,c('Ipsilateral 2Hz','Consensual 2Hz','Ipsilateral 1.6Hz','Consensual 1.6Hz'),pch=c(21,21,22,22),pt.bg=c('pink','black'),pt.lwd=3,pt.cex=1.6,box.lwd=2)
if(outputplot>0){dev.off()}  # this line goes after you've finished plotting (to output the example below, move it to the bottom of the script)





plotlims <- c(0,10,0,6)  # define the x and y limits of the plot (minx,maxx,miny,maxy)
ticklocsx <- seq(0,10,2)    # locations of tick marks on x axis
ticklocsy <- seq(0,6,1)    # locations of tick marks on y axis
ticklabelsx <-ticklocsx        # set labels for x ticks
ticklabelsy <- ticklocsy    # set labels for y ticks

for (block in 1:2){
for (eye in 1:2){
  if(outputplot==2){pdf(paste("Block",block,"Eye",eye,".pdf",sep=''), bg="transparent", height = 5.5, width = 5.5)}
  
plot(x=NULL,y=NULL,axes=FALSE, ann=FALSE, xlim=plotlims[1:2], ylim=plotlims[3:4])   # create an empty axis of the correct dimensions
axis(1, at=ticklocsx, tck=0.01, lab=F, lwd=2)     # plot tick marks (no labels)
axis(2, at=ticklocsy, tck=0.01, lab=F, lwd=2)
axis(3, at=ticklocsx, tck=0.01, lab=F, lwd=2)
axis(4, at=ticklocsy, tck=0.01, lab=F, lwd=2)
mtext(text = ticklabelsx, side = 1, at=ticklocsx)     # add the tick labels
mtext(text = ticklabelsy, side = 2, at=ticklocsy, line=0.2, las=1)  # the 'line' command moves away from the axis, the 'las' command rotates to vertical
box(lwd=2)      # draw a box around the graph
title(xlab="Time (s)", col.lab=rgb(0,0,0), line=1.2, cex.lab=1.5)    # titles for axes
title(ylab="Pupil diameter (mm)", col.lab=rgb(0,0,0), line=1.5, cex.lab=1.5)

for (cond in 1:60){lines(timeseq,(cond/10)+allwaveforms[block,eye,cond,],lwd=2)}
if(outputplot>0){dev.off()}  # this line goes after you've finished plotting (to output the example below, move it to the bottom of the script)

}}
  
  
  
frequencies <- (1:100)/10
plotlims <- c(0,10,0,6)  # define the x and y limits of the plot (minx,maxx,miny,maxy)
ticklocsx <- seq(0,10,1)    # locations of tick marks on x axis
ticklocsy <- seq(0,6,1)    # locations of tick marks on y axis
ticklabelsx <-ticklocsx        # set labels for x ticks
ticklabelsy <- ticklocsy    # set labels for y ticks
meanspectra <- abs(apply(allspectra,c(3,4),mean))

# for (block in 1:2){
#   for (eye in 1:2){
if(outputplot==2){pdf("meanspectra.pdf", bg="transparent", height = 5.5, width = 5.5)}

    plot(x=NULL,y=NULL,axes=FALSE, ann=FALSE, xlim=plotlims[1:2], ylim=plotlims[3:4])   # create an empty axis of the correct dimensions
    axis(1, at=ticklocsx, tck=0.01, lab=F, lwd=2)     # plot tick marks (no labels)
    axis(2, at=ticklocsy, tck=0.01, lab=F, lwd=2)
    axis(3, at=ticklocsx, tck=0.01, lab=F, lwd=2)
    axis(4, at=ticklocsy, tck=0.01, lab=F, lwd=2)
    mtext(text = ticklabelsx, side = 1, at=ticklocsx)     # add the tick labels
    mtext(text = ticklabelsy, side = 2, at=ticklocsy, line=0.2, las=1)  # the 'line' command moves away from the axis, the 'las' command rotates to vertical
    box(lwd=2)      # draw a box around the graph
    title(xlab="Frequency (Hz)", col.lab=rgb(0,0,0), line=1.2, cex.lab=1.5)    # titles for axes
    title(ylab="Pupil diameter (mm)", col.lab=rgb(0,0,0), line=1.5, cex.lab=1.5)
    
    for (cond in 1:60){
      toplot <- meanspectra[cond,2:101]*10
      toplot <- toplot - mean(toplot)
      lines(frequencies,(cond/10)+toplot,lwd=2)}
  # }}
lines(c(0.4,0.4),c(0,6),col='red')
lines(c(3.6,3.6),c(0,6),col='red')
lines(c(1.6,1.6),c(0,6),col='green')
lines(c(2,2),c(0,6),col='blue')
lines(c(3.2,3.2),c(0,6),col='green')
lines(c(4,4),c(0,6),col='blue')
if(outputplot>0){dev.off()}  # this line goes after you've finished plotting (to output the example below, move it to the bottom of the script)

meanSNR <- meanspectra
for (cond in 1:60){
for (f in 2:1199){
  meanSNR[cond,f] <- meanspectra[cond,f]/mean(meanspectra[cond,f+c(-1,1)])
}}

plotlims <- c(0,10,0,300)  # define the x and y limits of the plot (minx,maxx,miny,maxy)
ticklocsx <- seq(0,10,1)    # locations of tick marks on x axis
ticklocsy <- seq(0,300,50)    # locations of tick marks on y axis
ticklabelsx <-ticklocsx        # set labels for x ticks
ticklabelsy <- ticklocsy    # set labels for y ticks

if(outputplot==2){pdf("meanSNR.pdf", bg="transparent", height = 5.5, width = 5.5)}

plot(x=NULL,y=NULL,axes=FALSE, ann=FALSE, xlim=plotlims[1:2], ylim=plotlims[3:4])   # create an empty axis of the correct dimensions
axis(1, at=ticklocsx, tck=0.01, lab=F, lwd=2)     # plot tick marks (no labels)
axis(2, at=ticklocsy, tck=0.01, lab=F, lwd=2)
axis(3, at=ticklocsx, tck=0.01, lab=F, lwd=2)
axis(4, at=ticklocsy, tck=0.01, lab=F, lwd=2)
mtext(text = ticklabelsx, side = 1, at=ticklocsx)     # add the tick labels
mtext(text = ticklabelsy, side = 2, at=ticklocsy, line=0.2, las=1)  # the 'line' command moves away from the axis, the 'las' command rotates to vertical
box(lwd=2)      # draw a box around the graph
title(xlab="Frequency (Hz)", col.lab=rgb(0,0,0), line=1.2, cex.lab=1.5)    # titles for axes
title(ylab="SNR", col.lab=rgb(0,0,0), line=1.5, cex.lab=1.5)

for (cond in 1:60){
  toplot <- meanSNR[cond,2:101]
  toplot <- toplot - mean(toplot)
  lines(frequencies,(cond*5)+toplot,lwd=2)}
# }}
lines(c(0.4,0.4),c(0,300),col='red')
lines(c(3.6,3.6),c(0,300),col='red')
lines(c(1.6,1.6),c(0,300),col='green')
lines(c(2,2),c(0,300),col='blue')
lines(c(3.2,3.2),c(0,300),col='green')
lines(c(4,4),c(0,300),col='blue')

if(outputplot>0){dev.off()}  # this line goes after you've finished plotting (to output the example below, move it to the bottom of the script)


# targetSNRs <- meanSNR[,21]
# maskSNRs <- meanSNR[,17]
# 
# colvect <- c('red','blue','darkgreen','grey','purple','orange')
# plotlims <- c(12,40,0,40)  # define the x and y limits of the plot (minx,maxx,miny,maxy)
# ticklocsx <- seq(12,40,6)    # locations of tick marks on x axis
# ticklocsy <- seq(0,40,10)    # locations of tick marks on y axis
# ticklabelsx <-ticklocsx        # set labels for x ticks
# ticklabelsy <- ticklocsy    # set labels for y ticks
# 
# par(pty="s")  # make axis square
# plot(x=NULL,y=NULL,axes=FALSE, ann=FALSE, xlim=plotlims[1:2], ylim=plotlims[3:4])   # create an empty axis of the correct dimensions
# axis(1, at=ticklocsx, tck=0.01, lab=F, lwd=2)     # plot tick marks (no labels)
# axis(2, at=ticklocsy, tck=0.01, lab=F, lwd=2)
# axis(3, at=ticklocsx, tck=0.01, lab=F, lwd=2)
# axis(4, at=ticklocsy, tck=0.01, lab=F, lwd=2)
# mtext(text = ticklabelsx, side = 1, at=ticklocsx)     # add the tick labels
# mtext(text = ticklabelsy, side = 2, at=ticklocsy, line=0.2, las=1)  # the 'line' command moves away from the axis, the 'las' command rotates to vertical
# box(lwd=2)      # draw a box around the graph
# title(xlab="Target contrast (dB)", col.lab=rgb(0,0,0), line=1.2, cex.lab=1.5)    # titles for axes
# title(ylab="SNR", col.lab=rgb(0,0,0), line=1.5, cex.lab=1.5)
# 
# for (cond in 1:3){
#   temp <- targetSNRs[((10*(cond-1))+1):(10*cond)]
#   toplot <- abs(temp[1:5] + temp[6:10])/2
#   lines(contrastsdB,toplot, col=colvect[cond], lwd=3, cex=0.5)     # draw a line connecting the points
#   points(contrastsdB,toplot, pch = 21, col='black', bg=colvect[cond], cex=1.6, lwd=3)   # draw the data points themselves
# }  
# legend(12,40,c('Mon','Bin','Dich'),pch=21,pt.bg=colvect[1:3],pt.lwd=3,pt.cex=1.6,box.lwd=2)
# 
# 
# par(pty="s")  # make axis square
# plot(x=NULL,y=NULL,axes=FALSE, ann=FALSE, xlim=plotlims[1:2], ylim=plotlims[3:4])   # create an empty axis of the correct dimensions
# axis(1, at=ticklocsx, tck=0.01, lab=F, lwd=2)     # plot tick marks (no labels)
# axis(2, at=ticklocsy, tck=0.01, lab=F, lwd=2)
# axis(3, at=ticklocsx, tck=0.01, lab=F, lwd=2)
# axis(4, at=ticklocsy, tck=0.01, lab=F, lwd=2)
# mtext(text = ticklabelsx, side = 1, at=ticklocsx)     # add the tick labels
# mtext(text = ticklabelsy, side = 2, at=ticklocsy, line=0.2, las=1)  # the 'line' command moves away from the axis, the 'las' command rotates to vertical
# box(lwd=2)      # draw a box around the graph
# title(xlab="Target contrast (dB)", col.lab=rgb(0,0,0), line=1.2, cex.lab=1.5)    # titles for axes
# title(ylab="SNR", col.lab=rgb(0,0,0), line=1.5, cex.lab=1.5)
# 
# for (cond in 4:6){
#   temp <- targetSNRs[((10*(cond-1))+1):(10*cond)]
#   toplot <- abs(temp[1:5] + temp[6:10])/2
#   lines(contrastsdB,toplot, col=colvect[cond], lwd=3, cex=0.5)     # draw a line connecting the points
#   points(contrastsdB,toplot, pch = 21, col='black', bg=colvect[cond], cex=1.6, lwd=3)   # draw the data points themselves
# }  
# legend(12,40,c('Mon X','Bin X','Dich X'),pch=21,pt.bg=colvect[4:6],pt.lwd=3,pt.cex=1.6,box.lwd=2)
# 
# 
# par(pty="s")  # make axis square
# plot(x=NULL,y=NULL,axes=FALSE, ann=FALSE, xlim=plotlims[1:2], ylim=plotlims[3:4])   # create an empty axis of the correct dimensions
# axis(1, at=ticklocsx, tck=0.01, lab=F, lwd=2)     # plot tick marks (no labels)
# axis(2, at=ticklocsy, tck=0.01, lab=F, lwd=2)
# axis(3, at=ticklocsx, tck=0.01, lab=F, lwd=2)
# axis(4, at=ticklocsy, tck=0.01, lab=F, lwd=2)
# mtext(text = ticklabelsx, side = 1, at=ticklocsx)     # add the tick labels
# mtext(text = ticklabelsy, side = 2, at=ticklocsy, line=0.2, las=1)  # the 'line' command moves away from the axis, the 'las' command rotates to vertical
# box(lwd=2)      # draw a box around the graph
# title(xlab="Target contrast (dB)", col.lab=rgb(0,0,0), line=1.2, cex.lab=1.5)    # titles for axes
# title(ylab="SNR", col.lab=rgb(0,0,0), line=1.5, cex.lab=1.5)
# 
# for (cond in 4:6){
#   temp <- maskSNRs[((10*(cond-1))+1):(10*cond)]
#   toplot <- abs(temp[1:5] + temp[6:10])/2
#   lines(contrastsdB,toplot, col=colvect[cond], lwd=3, cex=0.5)     # draw a line connecting the points
#   points(contrastsdB,toplot, pch = 21, col='black', bg=colvect[cond], cex=1.6, lwd=3)   # draw the data points themselves
# }  
# legend(12,40,c('Mon X','Bin X','Dich X'),pch=21,pt.bg=colvect[4:6],pt.lwd=3,pt.cex=1.6,box.lwd=2)
# 
