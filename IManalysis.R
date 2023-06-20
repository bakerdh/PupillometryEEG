
# check which packages are installed, install the missing ones, and activate
packagelist <- c('knitr','remotes','tictoc','R.matlab','bookdown','grImport','png','tiff','pals','ez','gtools','signal','boot','quickpsy','rstan','coda','parallel','utils','osfr') # list of CRAN packages
missingpackages <- packagelist[!packagelist %in% installed.packages()[,1]]
if (length(missingpackages)>0){install.packages(missingpackages)}
if (!'FourierStats' %in% installed.packages()[,1]){remotes::install_github("bakerdh/FourierStats")}
packagelist <- c(packagelist,'FourierStats')
toinstall <- packagelist[which(!packagelist %in% (.packages()))]
invisible(lapply(toinstall,library,character.only=TRUE))

cleanup2D <- function(data){
  # helper function to clean up complex data by removing any data points more than 
  # SDthresh standard deviations from the mean (in Mahalanobis distance units)
  cartdata <- data.frame(Re(data),Im(data))
  output <- 0
  if (sum(cartdata)!=0){
    D <- sqrt(stats::mahalanobis(cartdata, colMeans(cartdata), cov(cartdata)))
    i <- which(D<SDthresh) 
    output <- mean(temp[i])}
  return(output)}

SDthresh <<- 3
nbootstraps <<- 1000

colpal <- c('#FE5000','#8783CF','#228B22','#808080','#6d008b','#8B8000')
targetelectrodes <- c('Oz','POz','O1','O2')

localdir <- 'local/'    # all files are stored in the project directory /local/ which git is told to ignore
if (!file.exists(localdir)){dir.create(localdir)}   # create a local directory to store data and outputs
rawdir <- 'local/rawdata/'
if (!file.exists(rawdir)){dir.create(rawdir)}   # create a local directory to store raw data
figdir <- 'Figures/'
if (!file.exists(figdir)){dir.create(figdir)}   # create a local directory to store figures
datadir <- 'local/spectra/'
if (!file.exists(datadir)){dir.create(datadir)}   # create a local directory to store processed data
  
  legaltriggers <- 1
  
  timeseq <- seq(1/120,10,length.out=120*10)
  timeseq2 <- seq(1/120,14,length.out=120*14)
  EEGtimes <- seq(1/1000,10,1/1000)
  targetindex <- (2*10)+1
  maskindex <- (1.6*10)+1
  targetindex2 <- (2*2*10)+1
  maskindex2 <- (2*1.6*10)+1
  
  for (participant in 1:30){
    
    # only download and process data if this has not been done already for this participant
    if (!file.exists(paste0(datadir,'P',100+participant,'spectra.RData'))){   
      
      # if there is no directory for this participant's data
      if (!file.exists(paste0(rawdir,'P',100+participant))){
        dir.create(paste0(rawdir,'P',100+participant))
        
        d <- dir(paste0(rawdir,'P',100+participant))
        
        if (length(d)==0){
          
          # if the tar file of this participant's data doesn't exist
          if (!file.exists(paste0(rawdir,'P',100+participant,'.tar'))){
            # download it from OSF
            if (!exists('osffiles')){
              osfnode <- 'x8u4v'
              osfproject <- osf_retrieve_node(osfnode)
              osffiles <- osf_ls_files(osfproject,n_max=300)}
            fid <- which(osffiles$name==paste0('P',100+participant,'.tar'))
            osf_download(osffiles[fid,],rawdir,progress=TRUE)
          }
          
          # then unzip the tar file
          untar(paste0(rawdir,'P',100+participant,'.tar'),exdir=paste0(rawdir,'P',100+participant))
          
          # and delete the tar file to save storage space
          file.remove(paste0(rawdir,'P',100+participant,'.tar'))
        }
      }
      
      d <- dir(paste0(rawdir,'P',100+participant),full.names=TRUE)
      
      pupiltargets <- array(0,dim=c(3,2,60))
      pupilmasks <- pupiltargets
      pupilwaveforms <- array(0,dim=c(3,2,60,120*14))
      pupilspectra <- array(0,dim=c(3,2,60,300))
      EEGtargets <- array(0,dim=c(3,64,60))
      EEGmasks <- EEGtargets
      EEGtargets2 <- EEGtargets
      EEGmasks2 <- EEGtargets
      EEGwaveforms <- array(0,dim=c(3,64,60,1000*14))
      EEGspectra <- array(0,dim=c(3,64,60,300))
      
      psychopyfiles <- dir(path=paste0(rawdir,'P',100+participant),pattern=paste0('P',100+participant,'_CRFstudy'), full.names = TRUE)
      infofiles <- dir(path=paste0(rawdir,'P',100+participant),pattern='*_info.csv', full.names = TRUE)
      EEGfiles <- dir(path=paste0(rawdir,'P',100+participant),pattern='*_EEG.csv.gz', full.names = TRUE)
      pupilfiles <- dir(path=paste0(rawdir,'P',100+participant),pattern='*_pupil_positions.csv', full.names = TRUE)
      
      for (block in 1:length(EEGfiles)){
        
        infofile <- read.csv(infofiles[block])
        pupilstarttime <- as.numeric(as.character(infofile[4,2]))
        localstarttime <- as.numeric(as.character(infofile[5,2]))
        
        psychoutput <- read.csv(psychopyfiles[block])
        
        trialtimes <- psychoutput$trialonset - pupilstarttime
        condorder <- psychoutput$condition
        
        EEGdata <- read.csv(EEGfiles[block])
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
            EEGtargets2[block,ch,condorder[cond]] <- fspec[targetindex2]
            EEGmasks2[block,ch,condorder[cond]] <- fspec[maskindex2]
            EEGspectra[block,ch,condorder[cond],1:300] <- fspec[1:300]
            
            trial <- EEGdata[(triggertimes[cond]+(-999:13000)),ch+2]
            EEGwaveforms[block,ch,condorder[cond],] <- trial - mean(trial[1:1000])
          }}
        
        
        pdata <- read.csv(pupilfiles[block])
        pdata2 <- pdata[,c(1,3,4,14)]
        pdata2[,1] <- pdata2[,1] - localstarttime
        pdata2 <- pdata2[which(pdata2[,3]>0.0),]
        
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
              pupilspectra[block,eye,condorder[cond],1:300] <- fspec[1:300]
            }
            
            a <- which(eyedata[,1]>(trialtimes[cond]-1))
            b <- which(eyedata[,1]<(trialtimes[cond]+13))
            i <- intersect(a,b)
            trial <- eyedata[i,]
            trial[,1] <- trial[,1] - trial[1,1]
            if (nrow(trial)>3){
              resampled <- interp1(trial[,1],trial[,4],timeseq2,method='linear',extrap=TRUE)
              pupilwaveforms[block,eye,condorder[cond],] <- resampled - mean(resampled[1:120])
            }
            
          }}
        
      }
      
      cleanmeansP <- matrix(0,nrow=6,ncol=5)
      cleanmasksP <- matrix(0,nrow=6,ncol=5)
      cleanmeansE <- array(0,dim=c(64,6,5))
      cleanmasksE <- array(0,dim=c(64,6,5))
      cleanmeansE2 <- array(0,dim=c(64,6,5))
      cleanmasksE2 <- array(0,dim=c(64,6,5))
      
      for (cond in 1:6){
        startindex <- (10*(cond-1))
        for (level in 1:5){
          temp <- c(pupiltargets[,,startindex+level],pupiltargets[,,startindex+level+5])
          cleanmeansP[cond,level] <- cleanup2D(temp)
        }
        
        temp <- c(pupilmasks[,,startindex+level],pupilmasks[,,startindex+level+5])
        cleanmasksP[cond,level] <- cleanup2D(temp)
        
        for (ch in 1:64){
          temp <- c(EEGtargets[,ch,startindex+level],EEGtargets[,ch,startindex+level+5])
          cleanmeansE[ch,cond,level] <- cleanup2D(temp)
          
          temp <- c(EEGmasks[,ch,startindex+level],EEGmasks[,ch,startindex+level+5])
          cleanmasksE[ch,cond,level] <- cleanup2D(temp)
          
          temp <- c(EEGtargets2[,ch,startindex+level],EEGtargets2[,ch,startindex+level+5])
          cleanmeansE2[ch,cond,level] <- cleanup2D(temp)
          
          temp <- c(EEGmasks2[,ch,startindex+level],EEGmasks2[,ch,startindex+level+5])
          cleanmasksE2[ch,cond,level] <- cleanup2D(temp) 
        }
      }
      
      electrodeindices <- match(targetelectrodes,electrodes)-2
      
      meanwavesP <- apply(pupilwaveforms[,,c(15,20),],4,mean,na.rm=TRUE)
      meanwavesE <- apply(EEGwaveforms[,electrodeindices,c(15,20),],4,mean,na.rm=TRUE)
      meanspectraP <- apply(pupilspectra[,,c(15,20),],4,mean,na.rm=TRUE)
      meanspectraE <- apply(EEGspectra[,electrodeindices,c(15,20),],4,mean,na.rm=TRUE)
      save(file=paste(datadir,'P',100+participant,'spectra.RData',sep=''),list=c('pupilspectra','EEGspectra','electrodeindices'))
      
    }
    
    
  }
  
  meanspectraP <- array(0,dim=c(30,6,5,300))
  meanspectraE <- array(0,dim=c(30,6,5,300))

  for (participant in 1:30){  
    
    load(paste(datadir,'P',100+participant,'spectra.RData',sep=''))
    
    for (lev in 1:5){
  meanspectraP[participant,1,lev,] <- apply(pupilspectra[,,c(lev,lev+5),],4,mean,na.rm=TRUE)
  meanspectraE[participant,1,lev,] <- apply(EEGspectra[,electrodeindices,c(lev,lev+5),],4,mean,na.rm=TRUE)
  
  meanspectraP[participant,2,lev,] <- apply(pupilspectra[,,c(lev,lev+5)+10,],4,mean,na.rm=TRUE)
  meanspectraE[participant,2,lev,] <- apply(EEGspectra[,electrodeindices,c(lev,lev+5)+10,],4,mean,na.rm=TRUE)
  
  meanspectraP[participant,3,lev,] <- apply(pupilspectra[,,c(lev,lev+5)+20,],4,mean,na.rm=TRUE)
  meanspectraE[participant,3,lev,] <- apply(EEGspectra[,electrodeindices,c(lev,lev+5)+20,],4,mean,na.rm=TRUE)
  
  meanspectraP[participant,4,lev,] <- apply(pupilspectra[,,c(lev,lev+5)+30,],4,mean,na.rm=TRUE)
  meanspectraE[participant,4,lev,] <- apply(EEGspectra[,electrodeindices,c(lev,lev+5)+30,],4,mean,na.rm=TRUE)

  meanspectraP[participant,5,lev,] <- apply(pupilspectra[,,c(lev,lev+5)+40,],4,mean,na.rm=TRUE)
  meanspectraE[participant,5,lev,] <- apply(EEGspectra[,electrodeindices,c(lev,lev+5)+40,],4,mean,na.rm=TRUE)

  meanspectraP[participant,6,lev,] <- apply(pupilspectra[,,c(lev,lev+5)+50,],4,mean,na.rm=TRUE)
  meanspectraE[participant,6,lev,] <- apply(EEGspectra[,electrodeindices,c(lev,lev+5)+50,],4,mean,na.rm=TRUE)
  
    }
  }
  
groupspecP <- abs(apply(meanspectraP,2:4,mean))
groupspecE <- abs(apply(meanspectraE,2:4,mean))

frequencies <- (0:299)/10

plotlims <- c(0,6,0,0.04) 
ticklocsx <- seq(0,6,1)    # locations of tick marks on x axis
ticklocsy <- seq(0,0.04,0.01)    # locations of tick marks on y axis
ticklabelsx <-ticklocsx        # set labels for x ticks
ticklabelsy <- ticklocsy    # set labels for y ticks

plot(x=NULL,y=NULL,axes=FALSE, ann=FALSE, xlim=plotlims[1:2], ylim=plotlims[3:4])   # create an empty axis of the correct dimensions
axis(1, at=ticklocsx, tck=0.01, lab=F, lwd=2)     # plot tick marks (no labels)
axis(2, at=ticklocsy, tck=0.01, lab=F, lwd=2)
mtext(text = ticklabelsx, side = 1, at=ticklocsx)     # add the tick labels
mtext(text = ticklabelsy, side = 2, at=ticklocsy, line=0.2, las=1)  
title(xlab="Frequency (Hz)", col.lab=rgb(0,0,0), line=1.2, cex.lab=1.5)  
title(ylab="Amplitude (mm)", col.lab=rgb(0,0,0), line=1.8, cex.lab=1.5)

# lines(frequencies[1:61],abs(groupspecP[1,5,1:61]), lwd=3, cex=0.5) 
# 
# lines(frequencies[1:61],abs(groupspecP[4,5,1:61]), lwd=3, cex=0.5) 

lines(frequencies[1:61],abs(groupspecP[6,4,1:61]), lwd=3, cex=0.5) 

lines(c(0.4,0.4),c(0,0.04),lty=2,col='blue')




plotlims <- c(0,6,0,0.5) 
ticklocsx <- seq(0,6,1)    # locations of tick marks on x axis
ticklocsy <- seq(0,0.5,0.1)    # locations of tick marks on y axis
ticklabelsx <-ticklocsx        # set labels for x ticks
ticklabelsy <- ticklocsy    # set labels for y ticks

plot(x=NULL,y=NULL,axes=FALSE, ann=FALSE, xlim=plotlims[1:2], ylim=plotlims[3:4])   # create an empty axis of the correct dimensions
axis(1, at=ticklocsx, tck=0.01, lab=F, lwd=2)     # plot tick marks (no labels)
axis(2, at=ticklocsy, tck=0.01, lab=F, lwd=2)
mtext(text = ticklabelsx, side = 1, at=ticklocsx)     # add the tick labels
mtext(text = ticklabelsy, side = 2, at=ticklocsy, line=0.2, las=1)  
title(xlab="Frequency (Hz)", col.lab=rgb(0,0,0), line=1.2, cex.lab=1.5)  
title(ylab="Amplitude (uV)", col.lab=rgb(0,0,0), line=1.8, cex.lab=1.5)

lines(frequencies[1:61],abs(groupspecE[5,4,1:61]), lwd=3, cex=0.5) 

lines(c(3.6,3.6),c(0,0.45),lty=2,col='blue')

