# WARNING: Each participant analysed will download around 1.3GB of data locally!
# this script downloads and analyses pupillometry and EEG data 

participant <- 'P302'
outputplot <- 0


localdir <- '~/Documents/local/'
if (!file.exists(localdir)){dir.create(localdir)}   # create a local directory to store data and outputs
EEGdir <- '~/Documents/local/EEG/'
if (!file.exists(EEGdir)){dir.create(EEGdir)}   # create a local directory to store EEG data
PPdir <- '~/Documents/local/Pupil/'
if (!file.exists(PPdir)){dir.create(PPdir)}   # create a local directory to store pupil data
Pydir <- '~/Documents/local/Psychopy/'
if (!file.exists(Pydir)){dir.create(Pydir)}   # create a local directory to store Psychopy output
figdir <- '~/Documents/local/Figures/'
if (!file.exists(figdir)){dir.create(figdir)}   # create a local directory to store figures
datadir <- '~/Documents/local/Data/'
if (!file.exists(datadir)){dir.create(datadir)}   # create a local directory to store processed data


packagelist <- c('signal','remotes','tictoc')
missingpackages <- packagelist[!packagelist %in% installed.packages()[,1]]
if (length(missingpackages)>0){install.packages(missingpackages)}
if (packagelist[!'osfr' %in% installed.packages()[,1]]){remotes::install_github("centerforopenscience/osfr")}
packagelist <- c(packagelist,'osfr')
toinstall <- packagelist[which(!packagelist %in% (.packages()))]
invisible(lapply(toinstall,library,character.only=TRUE))

tic()

osf_auth(token = '3dEYuhZNmwWbG3xuhRIiVRo0T2oniOkEP3Ip8i1LPG4PNjeTRln54eGNAG7oyTT9xozWwJ')
osfproject <- osf_retrieve_node("tbema")
componentlist <- osf_ls_nodes(osfproject)
EEGID <- match('Raw EEG data',as.character(unlist(componentlist[,1])))
EEGtoken <- componentlist[EEGID,2]
EEGfiles <- osf_ls_files(EEGtoken)
PPID <- match('Raw pupillometry data',as.character(unlist(componentlist[,1])))
PPtoken <- componentlist[PPID,2]
PPfiles <- osf_ls_files(PPtoken)
PyID <- match('Psychopy logs',as.character(unlist(componentlist[,1])))
Pytoken <- componentlist[PyID,2]
Pyfiles <- osf_ls_files(Pytoken)

for (n in 1:nrow(EEGfiles)){
  if (pmatch(participant,as.character(EEGfiles[n,1]),nomatch=0)){
    if (!file.exists(paste(EEGdir,as.character(EEGfiles[n,1]),sep=''))){
      osf_download(EEGfiles[n,],paste(EEGdir,as.character(EEGfiles[n,1]),sep=''))
    }}
  if (pmatch('header',as.character(EEGfiles[n,1]),nomatch=0)){
    if (!file.exists(paste(EEGdir,'headerfile.csv',sep=''))){
      osf_download(EEGfiles[n,],paste(EEGdir,'headerfile.csv',sep=''))
    }}  
}
for (n in 1:nrow(PPfiles)){
  if (pmatch(participant,as.character(PPfiles[n,1]),nomatch=0)){
    if (!file.exists(paste(PPdir,as.character(PPfiles[n,1]),sep=''))){
  osf_download(PPfiles[n,],paste(PPdir,as.character(PPfiles[n,1]),sep=''))
}}}
for (n in 1:nrow(Pyfiles)){
  if (pmatch(participant,as.character(Pyfiles[n,1]),nomatch=0)){
    if (!file.exists(paste(Pydir,as.character(Pyfiles[n,1]),sep=''))){
      osf_download(Pyfiles[n,],paste(Pydir,as.character(Pyfiles[n,1]),sep=''))
}}}

toc()


SDthresh <- 3
legaltriggers <- 1

pupiltargets <- array(0,dim=c(3,2,60))
pupilmasks <- pupiltargets
pupilwaveforms <- array(0,dim=c(3,2,60,1200))
pupilspectra <- pupilwaveforms
EEGtargets <- array(0,dim=c(3,64,60))
EEGmasks <- EEGtargets
EEGwaveforms <- array(0,dim=c(3,64,60,10000))
EEGspectra <- EEGwaveforms
timeseq <- seq(1/120,10,length.out=120*10)
EEGtimes <- seq(1/1000,10,1/1000)
targetindex <- (2*10)+1
maskindex <- (1.6*10)+1
showEEG <- 0

psychopyfiles <- dir(path=Pydir,pattern=participant)
# EEGfiles <- dir(path=paste(datadirectory,'EEG/',participant,'/',sep=''),pattern='*.csv.gz')

for (block in 1:3){
  infofile <- read.csv(paste(PPdir,participant,'_S',block,'_info.csv',sep=''))
  pupilstarttime <- as.numeric(as.character(infofile[4,2]))
  localstarttime <- as.numeric(as.character(infofile[5,2]))
  
  psychoutput <- read.csv(paste(Pydir,psychopyfiles[block],sep=''))

  trialtimes <- psychoutput$trialonset - pupilstarttime
  condorder <- psychoutput$condition
  
  if (file.exists(paste(EEGdir,participant,'_S',block,'_EEG.gz',sep=''))){
    showEEG <- 1
  EEGdata <- read.csv(paste(EEGdir,participant,'_S',block,'_EEG.gz',sep=''))
  electrodes <- colnames(EEGdata)
  
    triggertimes <- NULL
    counter <- 0
    lasttrigger <- -10000
    for (n in 1:nrow(EEGdata)){
      if(EEGdata$Trigger[n] %in% legaltriggers){
        if(n>(lasttrigger+10000)){
        counter <- counter + 1
        triggertimes[counter] <- n
        lasttrigger <- n
      }}
    }
  
    for (cond in 1:60){
      for (ch in 1:64){
        trial <- EEGdata[(triggertimes[cond]+(2001:12000)),ch+2]
        fspec <- (fft(trial)/length(trial))
        EEGtargets[block,ch,condorder[cond]] <- fspec[targetindex]
        EEGmasks[block,ch,condorder[cond]] <- fspec[maskindex]
        EEGwaveforms[block,ch,condorder[cond],] <- trial - mean(trial)
        EEGspectra[block,ch,condorder[cond],] <- fspec
      }}
  }
  
  pdata <- read.csv(paste(PPdir,participant,'_S',block,'_pupil_positions.csv',sep=''))
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
        pupiltargets[block,eye,condorder[cond]] <- fspec[targetindex]
        pupilmasks[block,eye,condorder[cond]] <- fspec[maskindex]
        pupilwaveforms[block,eye,condorder[cond],] <- resampled - mean(resampled)
        pupilspectra[block,eye,condorder[cond],] <- fspec
      }}}
}



cleanmeansP <- matrix(0,nrow=6,ncol=5)
cleanmasksP <- matrix(0,nrow=6,ncol=5)
cleanmeansE <- array(0,dim=c(64,6,5))
cleanmasksE <- array(0,dim=c(64,6,5))
for (cond in 1:6){
  startindex <- (10*(cond-1))
  for (level in 1:5){
    temp <- c(pupiltargets[,,startindex+level],pupiltargets[,,startindex+level+5])
    compmean <- mean(temp)
    absdiffs <- abs(temp - compmean)
    threshcut <- SDthresh*sd(absdiffs)
    cleanmeansP[cond,level] <- abs(mean(temp[which(absdiffs<threshcut)]))

    temp <- c(pupilmasks[,,startindex+level],pupiltargets[,,startindex+level+5])
    compmean <- mean(temp)
    absdiffs <- abs(temp - compmean)
    threshcut <- SDthresh*sd(absdiffs)
    cleanmasksP[cond,level] <- abs(mean(temp[which(absdiffs<threshcut)]))

    if (showEEG){
    for (ch in 1:64){
    temp <- c(EEGtargets[,ch,startindex+level],EEGtargets[,ch,startindex+level+5])
    compmean <- mean(temp)
    absdiffs <- abs(temp - compmean)
    threshcut <- SDthresh*sd(absdiffs)
    cleanmeansE[ch,cond,level] <- abs(mean(temp[which(absdiffs<threshcut)]))
    
    temp <- c(EEGmasks[,ch,startindex+level],EEGmasks[,ch,startindex+level+5])
    compmean <- mean(temp)
    absdiffs <- abs(temp - compmean)
    threshcut <- SDthresh*sd(absdiffs)
    cleanmasksE[ch,cond,level] <- abs(mean(temp[which(absdiffs<threshcut)]))
    }
    }
  }
}

contrastsdB <- 20*log10(c(6,12,24,48,96))
colvect <- c('red','blue','darkgreen','grey','purple','orange')
plotlims <- c(12,40,0,1)  # define the x and y limits of the plot (minx,maxx,miny,maxy)
ticklocsx <- seq(12,40,6)    # locations of tick marks on x axis
ticklocsy <- seq(0,1,0.2)    # locations of tick marks on y axis
ticklabelsx <-ticklocsx        # set labels for x ticks
ticklabelsy <- ticklocsy    # set labels for y ticks

if (showEEG){
targetelectrodes <- c('POz','Oz','O1','O2')
electrodeindices <- match(targetelectrodes,electrodes)-2

targetstoplot <- apply(cleanmeansE[electrodeindices,,],c(2,3),mean,na.rm=TRUE)
maskstoplot <- apply(cleanmasksE[electrodeindices,,],c(2,3),mean,na.rm=TRUE)

if(outputplot==2){pdf(paste(figdir,"CRF1e.pdf",sep=''), bg="transparent", height = 5.5, width = 5.5)}

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
title(ylab="Amplitude (µV) at 2Hz", col.lab=rgb(0,0,0), line=1.5, cex.lab=1.5)

for (cond in 1:3){
  lines(contrastsdB,targetstoplot[cond,], col=colvect[cond], lwd=3, cex=0.5)     # draw a line connecting the points
  points(contrastsdB,targetstoplot[cond,], pch = 21, col='black', bg=colvect[cond], cex=1.6, lwd=3)   # draw the data points themselves
}  
legend(12,1,c('Mon','Bin','Dich'),pch=21,pt.bg=colvect[1:3],pt.lwd=3,pt.cex=1.6,box.lwd=2)
if(outputplot>0){dev.off()}  # this line goes after you've finished plotting (to output the example below, move it to the bottom of the script)


if(outputplot==2){pdf(paste(figdir,"CRF2e.pdf",sep=''), bg="transparent", height = 5.5, width = 5.5)}

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
title(ylab="Amplitude (µV) at 2Hz", col.lab=rgb(0,0,0), line=1.5, cex.lab=1.5)

condlist <- c(1,5,6)
for (cond in 1:3){
  lines(contrastsdB,targetstoplot[condlist[cond],], col=colvect[condlist[cond]], lwd=3, cex=0.5)     # draw a line connecting the points
  points(contrastsdB,targetstoplot[condlist[cond],], pch = 21, col='black', bg=colvect[condlist[cond]], cex=1.6, lwd=3)   # draw the data points themselves
}  
legend(12,1,c('Mon','Bin X','Dich X'),pch=21,pt.bg=colvect[condlist],pt.lwd=3,pt.cex=1.6,box.lwd=2)
if(outputplot>0){dev.off()}  # this line goes after you've finished plotting (to output the example below, move it to the bottom of the script)



if(outputplot==2){pdf(paste(figdir,"CRF3e.pdf",sep=''), bg="transparent", height = 5.5, width = 5.5)}

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
title(ylab="Amplitude (µV) at 1.6Hz", col.lab=rgb(0,0,0), line=1.5, cex.lab=1.5)

for (cond in 4:6){
  lines(contrastsdB,maskstoplot[cond,], col=colvect[cond], lwd=3, cex=0.5)     # draw a line connecting the points
  points(contrastsdB,maskstoplot[cond,], pch = 21, col='black', bg=colvect[cond], cex=1.6, lwd=3)   # draw the data points themselves
}
legend(12,1,c('Mon X','Bin X','Dich X'),pch=21,pt.bg=colvect[4:6],pt.lwd=3,pt.cex=1.6,box.lwd=2)
if(outputplot>0){dev.off()}  # this line goes after you've finished plotting (to output the example below, move it to the bottom of the script)



hdata <- read.csv(paste(EEGdir,'headerfile.csv',sep=''),header=TRUE)

xpos <- 1:64
ypos <- 1:64
montageE <- toupper(as.character(hdata$Electrode))
for (ch in 1:64){
  i <- match(toupper(electrodes[ch+2]),montageE)
  xpos[ch] <- hdata$X_position[i]
  ypos[ch] <- hdata$Y_position[i]
}


v4Interp <- function(df, xo, yo, rmax = .75, gridRes = 67) {
  ## Create a function to perform Matlab's v4 interpolation.
  ## Takes as input a data-frame with columns x, y, and z (x co-ordinates, y co-ordinates, and amplitude)
  ## and variables xo and yo, the co-ordinates which will be use to create a grid for interpolation
  xo <- matrix(rep(xo,length(yo)),nrow = length(xo),ncol = length(yo))
  yo <- t(matrix(rep(yo,length(xo)),nrow = length(yo),ncol = length(xo)))
  xy <- df$x + df$y*sqrt(as.complex(-1))
  d <- matrix(rep(xy,length(xy)),nrow = length(xy), ncol = length(xy))
  d <- abs(d - t(d))
  diag(d) <- 1
  g <- (d^2) * (log(d)-1)   # Green's function.
  diag(g) <- 0
  weights <- qr.solve(g,df$z)
  xy <- t(xy)
  outmat <- matrix(nrow = gridRes,ncol = gridRes)
  for (i in 1:gridRes){
    for (j in 1:gridRes) {
      test4 <- abs((xo[i,j] + sqrt(as.complex(-1))*yo[i,j]) - xy)
      g <- (test4^2) * (log(test4)-1)
      outmat[i,j] <- g %*% weights}}
  outDf <- data.frame(x = xo[,1],outmat)
  names(outDf)[1:length(yo[1,])+1] <- yo[1,]
  return(outDf)}

## Create data frame to be used for interpolation - the function needs columns labelled x, y, and z
rmax <- 0.55   #specify a maximum boundary for the grid
gridRes <- 100 #specify the interpolation grid resolution

toplot <- cleanmeansE[,2,5]
toplot[which(is.na(toplot))] <- 0
testDat<- data.frame(x = xpos, y = -ypos, z = toplot)

#Create the interpolation grid
xo <- seq(min(-rmax, testDat$x), max(rmax, testDat$x), length = gridRes)
yo <- seq(max(rmax, testDat$y), min(-rmax, testDat$y), length = gridRes)

interpV4 <- v4Interp(testDat, xo, yo, rmax, gridRes)

zo2 <- as.matrix(interpV4[,2:ncol(interpV4)])

xo2 <- matrix(rep(xo,length(yo)),nrow = length(xo),ncol = length(yo))
yo2 <- t(matrix(rep(yo,length(xo)),nrow = length(yo),ncol = length(xo)))
outsidecircle <- sqrt(xo2^2 + yo2^2) > 0.51
zo2[outsidecircle] <- 0

ramp2 <- colorRamp(c("black","darkred","red","yellow","white"))  # create a ramp from one colour to another
colmatrix2 <- rgb(ramp2(seq(0, 1, length = 101)), max = 255)

# this section produces a scalp plot for the 5Hz component
if(outputplot==1){pdf(paste(figdir,"head1.pdf",sep=''), bg="transparent", height = 5.5, width = 5.5)}
if(outputplot==2){tiff(paste(figdir,"head1.tiff",sep=''), height = 600, width = 600, units="px", bg="white")}

plotlims <- c(-rmax,rmax,-rmax,rmax)  # define the x and y limits of the plot (minx,maxx,miny,maxy)
par(pty="s")  # make axis square
plot(x=NULL,y=NULL,ann=FALSE, axes=FALSE, xlim=plotlims[1:2], ylim=plotlims[3:4])   # create an empty axis of the correct dimensions
image(xo,xo,zo2,col=colmatrix2,add=TRUE,useRaster=TRUE)
maskx <- c(hdata$OutlineX[1:51]*2.2,hdata$OutlineX[51:1])
masky <- c(hdata$OutlineY[1:51]*2.2,hdata$OutlineY[51:1])
polygon(maskx,masky,border=NA,col="white")
maskx <- c(hdata$OutlineX[51:101]*2.2,hdata$OutlineX[101:51])
masky <- c(hdata$OutlineY[51:101]*2.2,hdata$OutlineY[101:51])
polygon(maskx,masky,border=NA,col="white")

blackelectrodes <- match(toupper(targetelectrodes),toupper(as.character(electrodes[3:66])))
points(xpos[blackelectrodes],ypos[blackelectrodes],pch=16,col="grey",cex=2)

lines(hdata$OutlineX,hdata$OutlineY,col="black",lwd=2)
lines(hdata$NoseX,hdata$NoseY,col="black",lwd=2)
lines(hdata$LearX,hdata$LearY,col="black",lwd=2)
lines(hdata$RearX,hdata$RearY,col="black",lwd=2)

if(outputplot>0){dev.off()}  # this line goes after you've finished plotting (to output the example below, move it to the bottom of the script)
}





plotlims <- c(12,40,0,0.04)  # define the x and y limits of the plot (minx,maxx,miny,maxy)
ticklocsx <- seq(12,40,6)    # locations of tick marks on x axis
ticklocsy <- seq(0,0.04,0.01)    # locations of tick marks on y axis
ticklabelsx <-ticklocsx        # set labels for x ticks
ticklabelsy <- ticklocsy    # set labels for y ticks

if(outputplot==2){pdf(paste(figdir,"CRF1p.pdf",sep=''), bg="transparent", height = 5.5, width = 5.5)}

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

for (cond in 1:3){
  lines(contrastsdB,cleanmeansP[cond,], col=colvect[cond], lwd=3, cex=0.5)     # draw a line connecting the points
  points(contrastsdB,cleanmeansP[cond,], pch = 21, col='black', bg=colvect[cond], cex=1.6, lwd=3)   # draw the data points themselves
}  
legend(12,0.04,c('Mon','Bin','Dich'),pch=21,pt.bg=colvect[1:3],pt.lwd=3,pt.cex=1.6,box.lwd=2)
if(outputplot>0){dev.off()}  # this line goes after you've finished plotting (to output the example below, move it to the bottom of the script)


if(outputplot==2){pdf(paste(figdir,"CRF2p.pdf",sep=''), bg="transparent", height = 5.5, width = 5.5)}

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

condlist <- c(1,5,6)
for (cond in 1:3){
  lines(contrastsdB,cleanmeansP[condlist[cond],], col=colvect[condlist[cond]], lwd=3, cex=0.5)     # draw a line connecting the points
  points(contrastsdB,cleanmeansP[condlist[cond],], pch = 21, col='black', bg=colvect[condlist[cond]], cex=1.6, lwd=3)   # draw the data points themselves
}  
legend(12,0.04,c('Mon','Bin X','Dich X'),pch=21,pt.bg=colvect[condlist],pt.lwd=3,pt.cex=1.6,box.lwd=2)
if(outputplot>0){dev.off()}  # this line goes after you've finished plotting (to output the example below, move it to the bottom of the script)



if(outputplot==2){pdf(paste(figdir,"CRF3p.pdf",sep=''), bg="transparent", height = 5.5, width = 5.5)}

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

for (cond in 4:6){
  lines(contrastsdB,cleanmasksP[cond,], col=colvect[cond], lwd=3, cex=0.5)     # draw a line connecting the points
  points(contrastsdB,cleanmasksP[cond,], pch = 21, col='black', bg=colvect[cond], cex=1.6, lwd=3)   # draw the data points themselves
}
legend(12,0.04,c('Mon X','Bin X','Dich X'),pch=21,pt.bg=colvect[4:6],pt.lwd=3,pt.cex=1.6,box.lwd=2)
if(outputplot>0){dev.off()}  # this line goes after you've finished plotting (to output the example below, move it to the bottom of the script)



# if(outputplot==2){pdf("consensual1.pdf", bg="transparent", height = 5.5, width = 5.5)}
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
# title(ylab="Amplitude (mm)", col.lab=rgb(0,0,0), line=1.5, cex.lab=1.5)
# 
# meanamps <- apply(pupiltargets,c(2,3),mean)
# 
# toplot <- abs(meanamps[1,1:5] + meanamps[2,6:10])/2
# lines(contrastsdB,toplot, col='pink', lwd=3, cex=0.5)     # draw a line connecting the points
# points(contrastsdB,toplot, pch = 21, col='black', bg='pink', cex=1.6, lwd=3)   # draw the data points themselves
# 
# toplot <- abs(meanamps[2,1:5] + meanamps[1,6:10])/2
# lines(contrastsdB,toplot, col='black', lwd=3, cex=0.5)     # draw a line connecting the points
# points(contrastsdB,toplot, pch = 21, col='black', bg='black', cex=1.6, lwd=3)   # draw the data points themselves
# 
# meanamps <- apply(pupilmasks,c(2,3),mean)
# toplot <- abs(meanamps[1,31:35] + meanamps[2,36:40])/2
# lines(contrastsdB,toplot, col='pink', lwd=3, cex=0.5)     # draw a line connecting the points
# points(contrastsdB,toplot, pch = 22, col='black', bg='pink', cex=1.6, lwd=3)   # draw the data points themselves
# 
# toplot <- abs(meanamps[2,31:35] + meanamps[1,36:40])/2
# lines(contrastsdB,toplot, col='black', lwd=3, cex=0.5)     # draw a line connecting the points
# points(contrastsdB,toplot, pch = 22, col='black', bg='black', cex=1.6, lwd=3)   # draw the data points themselves
# 
# legend(12,0.04,c('Ipsilateral 2Hz','Consensual 2Hz','Ipsilateral 1.6Hz','Consensual 1.6Hz'),pch=c(21,21,22,22),pt.bg=c('pink','black'),pt.lwd=3,pt.cex=1.6,box.lwd=2)
# if(outputplot>0){dev.off()}  # this line goes after you've finished plotting (to output the example below, move it to the bottom of the script)
# 
# 
# 
# if(outputplot==2){pdf("consensual2.pdf", bg="transparent", height = 5.5, width = 5.5)}
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
# title(ylab="Amplitude (mm)", col.lab=rgb(0,0,0), line=1.5, cex.lab=1.5)
# 
# meanamps <- apply(pupiltargets,c(2,3),mean)
# 
# toplot <- abs(meanamps[1,51:55] + meanamps[2,56:60])/2
# lines(contrastsdB,toplot, col='pink', lwd=3, cex=0.5)     # draw a line connecting the points
# points(contrastsdB,toplot, pch = 21, col='black', bg='pink', cex=1.6, lwd=3)   # draw the data points themselves
# 
# toplot <- abs(meanamps[2,51:55] + meanamps[1,56:60])/2
# lines(contrastsdB,toplot, col='black', lwd=3, cex=0.5)     # draw a line connecting the points
# points(contrastsdB,toplot, pch = 21, col='black', bg='black', cex=1.6, lwd=3)   # draw the data points themselves
# 
# meanamps <- apply(pupilmasks,c(2,3),mean)
# toplot <- abs(meanamps[1,51:55] + meanamps[2,56:60])/2
# lines(contrastsdB,toplot, col='pink', lwd=3, cex=0.5)     # draw a line connecting the points
# points(contrastsdB,toplot, pch = 22, col='black', bg='pink', cex=1.6, lwd=3)   # draw the data points themselves
# 
# toplot <- abs(meanamps[2,51:55] + meanamps[1,56:60])/2
# lines(contrastsdB,toplot, col='black', lwd=3, cex=0.5)     # draw a line connecting the points
# points(contrastsdB,toplot, pch = 22, col='black', bg='black', cex=1.6, lwd=3)   # draw the data points themselves
# 
# legend(28,0.04,c('Ipsilateral 2Hz','Consensual 2Hz','Ipsilateral 1.6Hz','Consensual 1.6Hz'),pch=c(21,21,22,22),pt.bg=c('pink','black'),pt.lwd=3,pt.cex=1.6,box.lwd=2)
# if(outputplot>0){dev.off()}  # this line goes after you've finished plotting (to output the example below, move it to the bottom of the script)
# 
# 
# 
# 
# 
# plotlims <- c(0,10,0,6)  # define the x and y limits of the plot (minx,maxx,miny,maxy)
# ticklocsx <- seq(0,10,2)    # locations of tick marks on x axis
# ticklocsy <- seq(0,6,1)    # locations of tick marks on y axis
# ticklabelsx <-ticklocsx        # set labels for x ticks
# ticklabelsy <- ticklocsy    # set labels for y ticks
# 
# for (block in 1:2){
#   for (eye in 1:2){
#     if(outputplot==2){pdf(paste("Block",block,"Eye",eye,".pdf",sep=''), bg="transparent", height = 5.5, width = 5.5)}
#     
#     plot(x=NULL,y=NULL,axes=FALSE, ann=FALSE, xlim=plotlims[1:2], ylim=plotlims[3:4])   # create an empty axis of the correct dimensions
#     axis(1, at=ticklocsx, tck=0.01, lab=F, lwd=2)     # plot tick marks (no labels)
#     axis(2, at=ticklocsy, tck=0.01, lab=F, lwd=2)
#     axis(3, at=ticklocsx, tck=0.01, lab=F, lwd=2)
#     axis(4, at=ticklocsy, tck=0.01, lab=F, lwd=2)
#     mtext(text = ticklabelsx, side = 1, at=ticklocsx)     # add the tick labels
#     mtext(text = ticklabelsy, side = 2, at=ticklocsy, line=0.2, las=1)  # the 'line' command moves away from the axis, the 'las' command rotates to vertical
#     box(lwd=2)      # draw a box around the graph
#     title(xlab="Time (s)", col.lab=rgb(0,0,0), line=1.2, cex.lab=1.5)    # titles for axes
#     title(ylab="Pupil diameter (mm)", col.lab=rgb(0,0,0), line=1.5, cex.lab=1.5)
#     
#     for (cond in 1:60){lines(timeseq,(cond/10)+pupilwaveforms[block,eye,cond,],lwd=2)}
#     if(outputplot>0){dev.off()}  # this line goes after you've finished plotting (to output the example below, move it to the bottom of the script)
#     
#   }}
# 
# 
# 
# frequencies <- (1:100)/10
# plotlims <- c(0,10,0,6)  # define the x and y limits of the plot (minx,maxx,miny,maxy)
# ticklocsx <- seq(0,10,1)    # locations of tick marks on x axis
# ticklocsy <- seq(0,6,1)    # locations of tick marks on y axis
# ticklabelsx <-ticklocsx        # set labels for x ticks
# ticklabelsy <- ticklocsy    # set labels for y ticks
# meanspectra <- abs(apply(pupilspectra,c(3,4),mean))
# 
# # for (block in 1:2){
# #   for (eye in 1:2){
# if(outputplot==2){pdf("meanspectra.pdf", bg="transparent", height = 5.5, width = 5.5)}
# 
# plot(x=NULL,y=NULL,axes=FALSE, ann=FALSE, xlim=plotlims[1:2], ylim=plotlims[3:4])   # create an empty axis of the correct dimensions
# axis(1, at=ticklocsx, tck=0.01, lab=F, lwd=2)     # plot tick marks (no labels)
# axis(2, at=ticklocsy, tck=0.01, lab=F, lwd=2)
# axis(3, at=ticklocsx, tck=0.01, lab=F, lwd=2)
# axis(4, at=ticklocsy, tck=0.01, lab=F, lwd=2)
# mtext(text = ticklabelsx, side = 1, at=ticklocsx)     # add the tick labels
# mtext(text = ticklabelsy, side = 2, at=ticklocsy, line=0.2, las=1)  # the 'line' command moves away from the axis, the 'las' command rotates to vertical
# box(lwd=2)      # draw a box around the graph
# title(xlab="Frequency (Hz)", col.lab=rgb(0,0,0), line=1.2, cex.lab=1.5)    # titles for axes
# title(ylab="Pupil diameter (mm)", col.lab=rgb(0,0,0), line=1.5, cex.lab=1.5)
# 
# for (cond in 1:60){
#   toplot <- meanspectra[cond,2:101]*10
#   toplot <- toplot - mean(toplot)
#   lines(frequencies,(cond/10)+toplot,lwd=2)}
# # }}
# lines(c(0.4,0.4),c(0,6),col='red')
# lines(c(3.6,3.6),c(0,6),col='red')
# lines(c(1.6,1.6),c(0,6),col='green')
# lines(c(2,2),c(0,6),col='blue')
# lines(c(3.2,3.2),c(0,6),col='green')
# lines(c(4,4),c(0,6),col='blue')
# if(outputplot>0){dev.off()}  # this line goes after you've finished plotting (to output the example below, move it to the bottom of the script)
# 
# meanSNR <- meanspectra
# for (cond in 1:60){
#   for (f in 2:1199){
#     meanSNR[cond,f] <- meanspectra[cond,f]/mean(meanspectra[cond,f+c(-1,1)])
#   }}
# 
# plotlims <- c(0,10,0,300)  # define the x and y limits of the plot (minx,maxx,miny,maxy)
# ticklocsx <- seq(0,10,1)    # locations of tick marks on x axis
# ticklocsy <- seq(0,300,50)    # locations of tick marks on y axis
# ticklabelsx <-ticklocsx        # set labels for x ticks
# ticklabelsy <- ticklocsy    # set labels for y ticks
# 
# if(outputplot==2){pdf("meanSNR.pdf", bg="transparent", height = 5.5, width = 5.5)}
# 
# plot(x=NULL,y=NULL,axes=FALSE, ann=FALSE, xlim=plotlims[1:2], ylim=plotlims[3:4])   # create an empty axis of the correct dimensions
# axis(1, at=ticklocsx, tck=0.01, lab=F, lwd=2)     # plot tick marks (no labels)
# axis(2, at=ticklocsy, tck=0.01, lab=F, lwd=2)
# axis(3, at=ticklocsx, tck=0.01, lab=F, lwd=2)
# axis(4, at=ticklocsy, tck=0.01, lab=F, lwd=2)
# mtext(text = ticklabelsx, side = 1, at=ticklocsx)     # add the tick labels
# mtext(text = ticklabelsy, side = 2, at=ticklocsy, line=0.2, las=1)  # the 'line' command moves away from the axis, the 'las' command rotates to vertical
# box(lwd=2)      # draw a box around the graph
# title(xlab="Frequency (Hz)", col.lab=rgb(0,0,0), line=1.2, cex.lab=1.5)    # titles for axes
# title(ylab="SNR", col.lab=rgb(0,0,0), line=1.5, cex.lab=1.5)
# 
# for (cond in 1:60){
#   toplot <- meanSNR[cond,2:101]
#   toplot <- toplot - mean(toplot)
#   lines(frequencies,(cond*5)+toplot,lwd=2)}
# # }}
# lines(c(0.4,0.4),c(0,300),col='red')
# lines(c(3.6,3.6),c(0,300),col='red')
# lines(c(1.6,1.6),c(0,300),col='green')
# lines(c(2,2),c(0,300),col='blue')
# lines(c(3.2,3.2),c(0,300),col='green')
# lines(c(4,4),c(0,300),col='blue')
# 
# if(outputplot>0){dev.off()}  # this line goes after you've finished plotting (to output the example below, move it to the bottom of the script)
