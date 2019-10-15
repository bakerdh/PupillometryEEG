

outputplot <- 2    # 0 means plot to plot window, 1 means separate pdfs, 2 means combined pdf (via ps files)


localdir <- 'local/'    # all files are stored in the project directory /local/ which git is told to ignore
if (!file.exists(localdir)){dir.create(localdir)}   # create a local directory to store data and outputs
figdir <- 'local/Figures/'
if (!file.exists(figdir)){dir.create(figdir)}   # create a local directory to store figures
datadir <- 'local/Data/'
if (!file.exists(datadir)){dir.create(datadir)}   # create a local directory to store processed data
EEGdir <- 'local/EEG/'
if (!file.exists(EEGdir)){dir.create(EEGdir)}   # create a local directory to store EEG data

packagelist <- c('tictoc','grImport','tiff')
missingpackages <- packagelist[!packagelist %in% installed.packages()[,1]]
if (length(missingpackages)>0){install.packages(missingpackages)}
toinstall <- packagelist[which(!packagelist %in% (.packages()))]
invisible(lapply(toinstall,library,character.only=TRUE))

tic()

d <- dir(datadir,pattern='summary.RData')

tempmasksP <- array(0,dim=c(length(d),6,5))
tempmeansP <- tempmasksP
tempmasksE <- array(0,dim=c(length(d),64,6,5))
tempmeansE <- tempmasksE
tempmeansE2 <- tempmeansE
tempmasksE2 <- tempmasksE
tempspectraP <- array(0,dim=c(length(d),60,300))
tempspectraE <- array(0,dim=c(length(d),64,60,300))

for (s in 1:length(d)){
load(paste(datadir,d[s],sep=''))
  
  tempmasksP[s,,] <- cleanmasksP
  tempmeansP[s,,] <- cleanmeansP
  tempspectraP[s,,] <- meanspectraP
  tempspectraE[s,,,] <- meanspectraE
  tempmeansE[s,,,] <- cleanmeansE
  tempmasksE[s,,,] <- cleanmasksE
  tempmeansE2[s,,,] <- cleanmeansE2
  tempmasksE2[s,,,] <- cleanmasksE2
}
cleanmasksP <- apply(tempmasksP,c(2,3),mean)
cleanmeansP <- apply(tempmeansP,c(2,3),mean)
meanspectraP <- apply(tempspectraP,c(2,3),mean)
meanspectraE <- apply(tempspectraE,c(2,3,4),mean)
cleanmeansE <- apply(tempmeansE,c(2,3,4),mean)
cleanmasksE <- apply(tempmasksE,c(2,3,4),mean)
cleanmeansE2 <- apply(tempmeansE2,c(2,3,4),mean)
cleanmasksE2 <- apply(tempmasksE2,c(2,3,4),mean)

cleanmasksPse <- apply(tempmasksP,c(2,3),sd)/sqrt(s)
cleanmeansPse <- apply(tempmeansP,c(2,3),sd)/sqrt(s)
cleanmeansEse <- apply(tempmeansE,c(2,3,4),sd)/sqrt(s)
cleanmasksEse <- apply(tempmasksE,c(2,3,4),sd)/sqrt(s)
cleanmeansE2se <- apply(tempmeansE2,c(2,3,4),sd)/sqrt(s)
cleanmasksE2se <- apply(tempmasksE2,c(2,3,4),sd)/sqrt(s)


hdata <- read.csv(paste(EEGdir,'headerfile.csv',sep=''),header=TRUE)

contrastsdB <- 20*log10(c(6,12,24,48,96))
colvect <- c('red','blue','darkgreen','grey','purple','orange')
plotlims <- c(12,40,0,2)  # define the x and y limits of the plot (minx,maxx,miny,maxy)
ticklocsx <- seq(12,40,6)    # locations of tick marks on x axis
ticklocsy <- seq(0,2,0.5)    # locations of tick marks on y axis
ticklabelsx <-ticklocsx        # set labels for x ticks
ticklabelsy <- ticklocsy    # set labels for y ticks


  targetelectrodes <- c('POz','Oz','O1','O2')
  electrodeindices <- match(targetelectrodes,electrodes)-2
  
  targetstoplot <- apply(cleanmeansE[electrodeindices,,],c(2,3),mean,na.rm=TRUE)
  maskstoplot <- apply(cleanmasksE[electrodeindices,,],c(2,3),mean,na.rm=TRUE)
  
  if(outputplot==1){pdf(paste(figdir,participant,"_CRF1e.pdf",sep=''), bg="transparent", height = 5.5, width = 5.5)}
  if(outputplot==2){postscript(paste("CRF1e.ps",sep=''), horizontal = FALSE, onefile = FALSE, paper = "special", height = 5.5, width = 5.5)}
  
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
  
  # arrows(contrastsdB,targetstoplot[cond,],x1=contrastsdB, y1=targetstoplot[cond,]-SEdata, length=0.015, angle=90, lwd=2, col='black')  # add lower error bar
  # arrows(contrastsdB,targetstoplot[cond,],x1=contrastsdB, y1=targetstoplot[cond,]+SEdata, length=0.015, angle=90, lwd=2, col='black')  # add upper error bar
  
  for (cond in 1:3){
    lines(contrastsdB,targetstoplot[cond,], col=colvect[cond], lwd=3, cex=0.5)     # draw a line connecting the points
    points(contrastsdB,targetstoplot[cond,], pch = 21, col='black', bg=colvect[cond], cex=1.6, lwd=3)   # draw the data points themselves
  }  
  legend(12,2,c('Mon','Bin','Dich'),pch=21,pt.bg=colvect[1:3],pt.lwd=3,pt.cex=1.6,box.lwd=2)
  if(outputplot>0){dev.off()}  # this line goes after you've finished plotting (to output the example below, move it to the bottom of the script)
  
  
  if(outputplot==1){pdf(paste(figdir,participant,"_CRF2e.pdf",sep=''), bg="transparent", height = 5.5, width = 5.5)}
  if(outputplot==2){postscript(paste("CRF2e.ps",sep=''), horizontal = FALSE, onefile = FALSE, paper = "special", height = 5.5, width = 5.5)}
  
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
  legend(12,2,c('Mon','Bin X','Dich X'),pch=21,pt.bg=colvect[condlist],pt.lwd=3,pt.cex=1.6,box.lwd=2)
  if(outputplot>0){dev.off()}  # this line goes after you've finished plotting (to output the example below, move it to the bottom of the script)
  
  
  if(outputplot==1){pdf(paste(figdir,participant,"_CRF3e.pdf",sep=''), bg="transparent", height = 5.5, width = 5.5)}
  if(outputplot==2){postscript(paste("CRF3e.ps",sep=''), horizontal = FALSE, onefile = FALSE, paper = "special", height = 5.5, width = 5.5)}
  
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
  legend(12,2,c('Mon X','Bin X','Dich X'),pch=21,pt.bg=colvect[4:6],pt.lwd=3,pt.cex=1.6,box.lwd=2)
  if(outputplot>0){dev.off()}  # this line goes after you've finished plotting (to output the example below, move it to the bottom of the script)
  
  
  
  plotlims <- c(12,40,0,1)  # define the x and y limits of the plot (minx,maxx,miny,maxy)
  ticklocsx <- seq(12,40,6)    # locations of tick marks on x axis
  ticklocsy <- seq(0,1,0.2)    # locations of tick marks on y axis
  ticklabelsx <-ticklocsx        # set labels for x ticks
  ticklabelsy <- ticklocsy    # set labels for y ticks
  
  targetstoplot <- apply(cleanmeansE2[electrodeindices,,],c(2,3),mean,na.rm=TRUE)
  maskstoplot <- apply(cleanmasksE2[electrodeindices,,],c(2,3),mean,na.rm=TRUE)
  
  if(outputplot==1){pdf(paste(figdir,participant,"_CRF1e2.pdf",sep=''), bg="transparent", height = 5.5, width = 5.5)}
  if(outputplot==2){postscript(paste("CRF1e2.ps",sep=''), horizontal = FALSE, onefile = FALSE, paper = "special", height = 5.5, width = 5.5)}
  
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
  title(ylab="Amplitude (µV) at 4Hz", col.lab=rgb(0,0,0), line=1.5, cex.lab=1.5)
  
  for (cond in 1:3){
    lines(contrastsdB,targetstoplot[cond,], col=colvect[cond], lwd=3, cex=0.5)     # draw a line connecting the points
    points(contrastsdB,targetstoplot[cond,], pch = 21, col='black', bg=colvect[cond], cex=1.6, lwd=3)   # draw the data points themselves
  }  
  legend(12,1,c('Mon','Bin','Dich'),pch=21,pt.bg=colvect[1:3],pt.lwd=3,pt.cex=1.6,box.lwd=2)
  if(outputplot>0){dev.off()}  # this line goes after you've finished plotting (to output the example below, move it to the bottom of the script)
  
  
  if(outputplot==1){pdf(paste(figdir,participant,"_CRF2e2.pdf",sep=''), bg="transparent", height = 5.5, width = 5.5)}
  if(outputplot==2){postscript(paste("CRF2e2.ps",sep=''), horizontal = FALSE, onefile = FALSE, paper = "special", height = 5.5, width = 5.5)}
  
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
  title(ylab="Amplitude (µV) at 4Hz", col.lab=rgb(0,0,0), line=1.5, cex.lab=1.5)
  
  condlist <- c(1,5,6)
  for (cond in 1:3){
    lines(contrastsdB,targetstoplot[condlist[cond],], col=colvect[condlist[cond]], lwd=3, cex=0.5)     # draw a line connecting the points
    points(contrastsdB,targetstoplot[condlist[cond],], pch = 21, col='black', bg=colvect[condlist[cond]], cex=1.6, lwd=3)   # draw the data points themselves
  }  
  legend(12,1,c('Mon','Bin X','Dich X'),pch=21,pt.bg=colvect[condlist],pt.lwd=3,pt.cex=1.6,box.lwd=2)
  if(outputplot>0){dev.off()}  # this line goes after you've finished plotting (to output the example below, move it to the bottom of the script)
  
  
  if(outputplot==1){pdf(paste(figdir,participant,"_CRF3e2.pdf",sep=''), bg="transparent", height = 5.5, width = 5.5)}
  if(outputplot==2){postscript(paste("CRF3e2.ps",sep=''), horizontal = FALSE, onefile = FALSE, paper = "special", height = 5.5, width = 5.5)}
  
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
  title(ylab="Amplitude (µV) at 3.2Hz", col.lab=rgb(0,0,0), line=1.5, cex.lab=1.5)
  
  for (cond in 4:6){
    lines(contrastsdB,maskstoplot[cond,], col=colvect[cond], lwd=3, cex=0.5)     # draw a line connecting the points
    points(contrastsdB,maskstoplot[cond,], pch = 21, col='black', bg=colvect[cond], cex=1.6, lwd=3)   # draw the data points themselves
  }
  legend(12,1,c('Mon X','Bin X','Dich X'),pch=21,pt.bg=colvect[4:6],pt.lwd=3,pt.cex=1.6,box.lwd=2)
  if(outputplot>0){dev.off()}  # this line goes after you've finished plotting (to output the example below, move it to the bottom of the script)
  
  
  
  
  
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
  if(outputplot==1){pdf(paste(figdir,participant,"_head1.pdf",sep=''), bg="transparent", height = 5.5, width = 5.5)}
  if(outputplot==2){tiff(paste("head1.tiff",sep=''), height = 600, width = 600, units="px", bg="white")}
  
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
  
  
  toplot <- cleanmeansE2[,2,5]
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
  
  # this section produces a scalp plot for the 5Hz component
  if(outputplot==1){pdf(paste(figdir,participant,"_head2.pdf",sep=''), bg="transparent", height = 5.5, width = 5.5)}
  if(outputplot==2){tiff(paste("head2.tiff",sep=''), height = 600, width = 600, units="px", bg="white")}
  
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






plotlims <- c(12,40,0,0.04)  # define the x and y limits of the plot (minx,maxx,miny,maxy)
ticklocsx <- seq(12,40,6)    # locations of tick marks on x axis
ticklocsy <- seq(0,0.04,0.01)    # locations of tick marks on y axis
ticklabelsx <-ticklocsx        # set labels for x ticks
ticklabelsy <- ticklocsy    # set labels for y ticks

if(outputplot==1){pdf(paste(figdir,participant,"_CRF1p.pdf",sep=''), bg="transparent", height = 5.5, width = 5.5)}
if(outputplot==2){postscript(paste("CRF1p.ps",sep=''), horizontal = FALSE, onefile = FALSE, paper = "special", height = 5.5, width = 5.5)}

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

if(outputplot==1){pdf(paste(figdir,participant,"_CRF2p.pdf",sep=''), bg="transparent", height = 5.5, width = 5.5)}
if(outputplot==2){postscript(paste("CRF2p.ps",sep=''), horizontal = FALSE, onefile = FALSE, paper = "special", height = 5.5, width = 5.5)}

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


if(outputplot==1){pdf(paste(figdir,participant,"_CRF3p.pdf",sep=''), bg="transparent", height = 5.5, width = 5.5)}
if(outputplot==2){postscript(paste("CRF3p.ps",sep=''), horizontal = FALSE, onefile = FALSE, paper = "special", height = 5.5, width = 5.5)}

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



if(outputplot==1){pdf(paste(figdir,participant,"_SpecP.pdf",sep=''), bg="transparent", height = 5.5, width = 8)}
if(outputplot==2){postscript(paste("SpecP.ps",sep=''), horizontal = FALSE, onefile = FALSE, paper = "special", height = 5.5, width = 8)}

frequencies <- (0:299)/10

plotlims <- c(1,5,0,0.1)  # define the x and y limits of the plot (minx,maxx,miny,maxy)
ticklocsx <- seq(1,5,1)    # locations of tick marks on x axis
ticklocsy <- seq(0,0.1,0.02)    # locations of tick marks on y axis
ticklabelsx <-ticklocsx        # set labels for x ticks
ticklabelsy <- ticklocsy    # set labels for y ticks

plot(x=NULL,y=NULL,axes=FALSE, ann=FALSE, xlim=plotlims[1:2], ylim=plotlims[3:4])   # create an empty axis of the correct dimensions
axis(1, at=ticklocsx, tck=0.01, lab=F, lwd=2)     # plot tick marks (no labels)
axis(2, at=ticklocsy, tck=0.01, lab=F, lwd=2)
mtext(text = ticklabelsx, side = 1, at=ticklocsx)     # add the tick labels
mtext(text = ticklabelsy, side = 2, at=ticklocsy, line=0.2, las=1)  # the 'line' command moves away from the axis, the 'las' command rotates to vertical
title(xlab="Frequency (Hz)", col.lab=rgb(0,0,0), line=1.2, cex.lab=1.5)    # titles for axes
title(ylab="Amplitude (mm)", col.lab=rgb(0,0,0), line=1.5, cex.lab=1.5)

ftoplot <- abs((meanspectraP[15,] + meanspectraP[20,])/2)
lines(frequencies[11:51],ftoplot[11:51], col='black', lwd=3, cex=0.5)     # draw a line connecting the points

if(outputplot>0){dev.off()}  # this line goes after you've finished plotting (to output the example below, move it to the bottom of the script)



if(outputplot==1){pdf(paste(figdir,participant,"_SpecE.pdf",sep=''), bg="transparent", height = 5.5, width = 8)}
if(outputplot==2){postscript(paste("SpecE.ps",sep=''), horizontal = FALSE, onefile = FALSE, paper = "special", height = 5.5, width = 8)}

frequencies <- (0:299)/10

plotlims <- c(1,5,0,2)  # define the x and y limits of the plot (minx,maxx,miny,maxy)
ticklocsx <- seq(1,5,1)    # locations of tick marks on x axis
ticklocsy <- seq(0,2,0.5)    # locations of tick marks on y axis
ticklabelsx <-ticklocsx        # set labels for x ticks
ticklabelsy <- ticklocsy    # set labels for y ticks

plot(x=NULL,y=NULL,axes=FALSE, ann=FALSE, xlim=plotlims[1:2], ylim=plotlims[3:4])   # create an empty axis of the correct dimensions
axis(1, at=ticklocsx, tck=0.01, lab=F, lwd=2)     # plot tick marks (no labels)
axis(2, at=ticklocsy, tck=0.01, lab=F, lwd=2)
mtext(text = ticklabelsx, side = 1, at=ticklocsx)     # add the tick labels
mtext(text = ticklabelsy, side = 2, at=ticklocsy, line=0.2, las=1)  # the 'line' command moves away from the axis, the 'las' command rotates to vertical
title(xlab="Frequency (Hz)", col.lab=rgb(0,0,0), line=1.2, cex.lab=1.5)    # titles for axes
title(ylab="Amplitude (µV)", col.lab=rgb(0,0,0), line=1.5, cex.lab=1.5)

ftoplot <- abs(colMeans((meanspectraE[28:31,15,] + meanspectraE[28:31,20,])/2))
lines(frequencies[11:51],ftoplot[11:51], col='black', lwd=3, cex=0.5)     # draw a line connecting the points

if(outputplot>0){dev.off()}  # this line goes after you've finished plotting (to output the example below, move it to the bottom of the script)



# if we're outputting as postscript, it's because we're merging to create the full figure
if(outputplot==2){
  
  PostScriptTrace(paste('CRF1p.ps',sep=''))
  p1 <- readPicture(paste('CRF1p.ps.xml',sep=''))
  PostScriptTrace(paste('CRF2p.ps',sep=''))
  p2 <- readPicture(paste('CRF2p.ps.xml',sep=''))
  PostScriptTrace(paste('CRF3p.ps',sep=''))
  p3 <- readPicture(paste('CRF3p.ps.xml',sep=''))
  PostScriptTrace(paste('SpecP.ps',sep=''))
  p4 <- readPicture(paste('SpecP.ps.xml',sep=''))
  
  PostScriptTrace(paste('CRF1e.ps',sep=''))
  e1 <- readPicture(paste('CRF1e.ps.xml',sep=''))
  PostScriptTrace(paste('CRF2e.ps',sep=''))
  e2 <- readPicture(paste('CRF2e.ps.xml',sep=''))
  PostScriptTrace(paste('CRF3e.ps',sep=''))
  e3 <- readPicture(paste('CRF3e.ps.xml',sep=''))
  PostScriptTrace(paste('CRF1e2.ps',sep=''))
  e4 <- readPicture(paste('CRF1e2.ps.xml',sep=''))
  PostScriptTrace(paste('CRF2e2.ps',sep=''))
  e5 <- readPicture(paste('CRF2e2.ps.xml',sep=''))
  PostScriptTrace(paste('CRF3e2.ps',sep=''))
  e6 <- readPicture(paste('CRF3e2.ps.xml',sep=''))
  PostScriptTrace(paste('SpecE.ps',sep=''))
  e7 <- readPicture(paste('SpecE.ps.xml',sep=''))
  
  h1 <- readTIFF('head1.tiff')
  h2 <- readTIFF('head2.tiff')
  
  pdf(paste(figdir,"MeanSummary.pdf",sep=''), bg="transparent", height = 10, width = 15)
  par(mar=c(0.1,0.1,0.1,0.1))
  plot(x=NULL,y=NULL,axes=FALSE, ann=FALSE, xlim=c(0,1), ylim=c(0,1))   # create an empty axis of the correct dimensions
  
  aspratio <- 10/15  # this is the aspect ratio of the output pdf
  imwidth <- 0.24
  xstart <- 0.7
  ystart <- 0.4
  rasterImage(h1,xstart,ystart,xstart+imwidth*aspratio,ystart+imwidth) # insert the head plot first so the white border doesn't overlap the other graphs
  xstart <- 0.85
  ystart <- 0.4
  rasterImage(h2,xstart,ystart,xstart+imwidth*aspratio,ystart+imwidth) # insert the head plot first so the white border doesn't overlap the other graphs
  
  grid.picture(p1,x=0.12,y=0.8,width=0.2,height=1)
  grid.picture(p2,x=0.32,y=0.8,width=0.2,height=1)
  grid.picture(p3,x=0.52,y=0.8,width=0.2,height=1)
  grid.picture(p4,x=0.8,y=0.8,width=0.35,height=1)
  
  grid.picture(e1,x=0.12,y=0.5,width=0.2,height=1)
  grid.picture(e2,x=0.32,y=0.5,width=0.2,height=1)
  grid.picture(e3,x=0.52,y=0.5,width=0.2,height=1)
  
  grid.picture(e4,x=0.12,y=0.2,width=0.2,height=1)
  grid.picture(e5,x=0.32,y=0.2,width=0.2,height=1)
  grid.picture(e6,x=0.52,y=0.2,width=0.2,height=1)
  grid.picture(e7,x=0.8,y=0.2,width=0.35,height=1)
  
  dev.off()
  file.remove(c('CRF1p.ps','CRF2p.ps','CRF3p.ps','CRF1e.ps','CRF2e.ps','CRF3e.ps','CRF1e2.ps','CRF2e2.ps','CRF3e2.ps','SpecE.ps','SpecP.ps','head1.tiff','head2.tiff'))
  file.remove(c('CRF1p.ps.xml','CRF2p.ps.xml','CRF3p.ps.xml','CRF1e.ps.xml','CRF2e.ps.xml','CRF3e.ps.xml','CRF1e2.ps.xml','CRF2e2.ps.xml','CRF3e2.ps.xml','SpecE.ps.xml','SpecP.ps.xml'))
  file.remove(c('captureCRF1p.ps','captureCRF2p.ps','captureCRF3p.ps','captureCRF1e.ps','captureCRF2e.ps','captureCRF3e.ps','captureCRF1e2.ps','captureCRF2e2.ps','captureCRF3e2.ps','captureSpecE.ps','captureSpecP.ps'))
}

toc()