# script to reverse engineer trigger timing using coincidence of blinks between frontal electrodes and eye tracking data
# DHB 12/10/19

# first section is proof of concept on a participant without missing triggers
# align eye tracking and EEG data based on blinks, then check if re-aligned triggers coincide

participant <- 'P304'
block <- 1

packagelist <- c('signal','tictoc','R.utils')
missingpackages <- packagelist[!packagelist %in% installed.packages()[,1]]
if (length(missingpackages)>0){install.packages(missingpackages)}
toinstall <- packagelist[which(!packagelist %in% (.packages()))]
invisible(lapply(toinstall,library,character.only=TRUE))


localdir <- 'local/'
EEGdir <- '/local/EEG/'
PPdir <- 'local/Pupil/'
Pydir <- 'local/Psychopy/'
figdir <- 'local/Figures/'
datadir <- 'local/Data/'


infofile <- read.csv(paste(PPdir,participant,'_S',block,'_info.csv',sep=''))
pupilstarttime <- as.numeric(as.character(infofile[4,2]))
localstarttime <- as.numeric(as.character(infofile[5,2]))

psychopyfiles <- dir(path=Pydir,pattern=participant)
psychoutput <- read.csv(paste(Pydir,psychopyfiles[block],sep=''))

trialtimes <- psychoutput$trialonset - pupilstarttime
condorder <- psychoutput$condition

EEGdata <- read.csv(paste(EEGdir,participant,'_S',block,'_EEG.gz',sep=''))
electrodes <- colnames(EEGdata)

legaltriggers <- 1
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

pdata <- read.csv(paste(PPdir,participant,'_S',block,'_pupil_positions.csv',sep=''))
pdata2 <- pdata[,c(1,3,4,14)]
pdata2[,1] <- pdata2[,1] - localstarttime
eyedataL <- pdata2[which(pdata2[,2]==(1)),]
eyedataR <- pdata2[which(pdata2[,2]==(0)),]

timeseq <- seq(0,max(pdata2[,1]),1/120)
resampledL <- interp1(eyedataL[,1],eyedataL[,3],timeseq,method='linear',extrap=TRUE)
resampledR <- interp1(eyedataR[,1],eyedataR[,3],timeseq,method='linear',extrap=TRUE)
blinksignal <- 1-round(resampledL*resampledR)

eegtimes <- EEGdata$Time/1000
eegsignal1 <- EEGdata$Fp1
eegsignal1 <- eegsignal1 - min(eegsignal1)
eegsignal1 <- eegsignal1/max(eegsignal1)
eegsignal2 <- EEGdata$Fp2
eegsignal2 <- eegsignal2 - min(eegsignal2)
eegsignal2 <- eegsignal2/max(eegsignal2)

timediff <- (triggertimes[1]/1000) - trialtimes[1]
# eegtimes <- eegtimes - timediff
triggertimesadjusted <- (triggertimes/1000) - timediff

EEGresamp1 <- interp1(eegtimes,eegsignal1,timeseq,method='linear',extrap=TRUE)
EEGresamp2 <- interp1(eegtimes,eegsignal2,timeseq,method='linear',extrap=TRUE)
meanEEG <- (EEGresamp1 + EEGresamp2)/2

xcor <- ccf(meanEEG,blinksignal, lag.max = 2000)
xcor$lag <- xcor$lag[2000:4001]
xcor$acf <- xcor$acf[2000:4001]

firstderiv <- diff(xcor$acf)
plot(xcor$lag[2:length(xcor$lag)],firstderiv,type='l')
meanlag <- ((xcor$lag[which(firstderiv==max(firstderiv))]) + (xcor$lag[which(firstderiv==min(firstderiv))]))/2
lines(c(meanlag,meanlag),c(-1,1),col='blue')
lagval <- meanlag/120
# lagval <- (xcor$lag[which(firstderiv==max(firstderiv))])/120
eegtimesadjusted <- eegtimes - lagval
triggertimesadjusted2 <- (triggertimes/1000) - lagval

EEGresamp1 <- interp1(eegtimesadjusted,eegsignal1,timeseq,method='linear',extrap=TRUE)
EEGresamp2 <- interp1(eegtimesadjusted,eegsignal2,timeseq,method='linear',extrap=TRUE)
meanEEG <- (EEGresamp1 + EEGresamp2)/2

timestoplot <- 12000:24000
plot(timeseq[timestoplot],blinksignal[timestoplot],type='l',col='black')
for (n in 1:length(trialtimes)){lines(trialtimes[c(n,n)],c(0,1),col='red')}
for (n in 1:length(triggertimes)){lines(triggertimesadjusted[c(n,n)],c(0,1),col='blue',lty=2)}
for (n in 1:length(triggertimes)){lines(triggertimesadjusted2[c(n,n)],c(0,1),col='orange',lty=2)}
lines(timeseq[timestoplot],0.4+meanEEG[timestoplot],col='green')

triggertimesadjusted-triggertimesadjusted2





participant <- 'P302'
sourceEEG <- '~/Desktop/Pupildata/EEG/'

for (block in 1:3){

infofile <- read.csv(paste(PPdir,participant,'_S',block,'_info.csv',sep=''))
pupilstarttime <- as.numeric(as.character(infofile[4,2]))
localstarttime <- as.numeric(as.character(infofile[5,2]))

psychopyfiles <- dir(path=Pydir,pattern=participant)
psychoutput <- read.csv(paste(Pydir,psychopyfiles[block],sep=''))

trialtimes <- psychoutput$trialonset - pupilstarttime
condorder <- psychoutput$condition

EEGdata <- read.csv(paste(sourceEEG,participant,'_S',block,'_EEG.gz',sep=''))
electrodes <- colnames(EEGdata)

pdata <- read.csv(paste(PPdir,participant,'_S',block,'_pupil_positions.csv',sep=''))
pdata2 <- pdata[,c(1,3,4,14)]
pdata2[,1] <- pdata2[,1] - localstarttime
eyedataL <- pdata2[which(pdata2[,2]==(1)),]
eyedataR <- pdata2[which(pdata2[,2]==(0)),]

timeseq <- seq(0,max(pdata2[,1]),1/120)
resampledL <- interp1(eyedataL[,1],eyedataL[,3],timeseq,method='linear',extrap=TRUE)
resampledR <- interp1(eyedataR[,1],eyedataR[,3],timeseq,method='linear',extrap=TRUE)
blinksignal <- 1-round(resampledL*resampledR)

eegtimes <- EEGdata$Time/1000
eegsignal1 <- EEGdata$Fp1
eegsignal1 <- eegsignal1 - min(eegsignal1)
eegsignal1 <- eegsignal1/max(eegsignal1)
eegsignal2 <- EEGdata$Fp2
eegsignal2 <- eegsignal2 - min(eegsignal2)
eegsignal2 <- eegsignal2/max(eegsignal2)

# timediff <- (triggertimes[1]/1000) - trialtimes[1]
# # eegtimes <- eegtimes - timediff
# triggertimesadjusted <- (triggertimes/1000) - timediff

EEGresamp1 <- interp1(eegtimes,eegsignal1,timeseq,method='linear',extrap=TRUE)
EEGresamp2 <- interp1(eegtimes,eegsignal2,timeseq,method='linear',extrap=TRUE)
meanEEG <- (EEGresamp1 + EEGresamp2)/2

xcor <- ccf(meanEEG,blinksignal, lag.max = 2000)
xcor$lag <- xcor$lag[2000:4001]
xcor$acf <- xcor$acf[2000:4001]

firstderiv <- diff(xcor$acf)
plot(xcor$lag[2:length(xcor$lag)],firstderiv,type='l')
meanlag <- ((xcor$lag[which(firstderiv==max(firstderiv[500:1000]))]) + (xcor$lag[which(firstderiv==min(firstderiv[500:1000]))]))/2
lines(c(meanlag,meanlag),c(-1,1),col='blue')
lagval <- meanlag/120
# lagval <- (xcor$lag[which(firstderiv==max(firstderiv))])/120
eegtimesadjusted <- eegtimes - lagval
# triggertimesadjusted2 <- (triggertimes/1000) - lagval

EEGresamp1 <- interp1(eegtimesadjusted,eegsignal1,timeseq,method='linear',extrap=TRUE)
EEGresamp2 <- interp1(eegtimesadjusted,eegsignal2,timeseq,method='linear',extrap=TRUE)
meanEEG <- (EEGresamp1 + EEGresamp2)/2

timestoplot <- 1:(100*120)
plot(timeseq[timestoplot],blinksignal[timestoplot],type='l',col='black')
for (n in 1:length(trialtimes)){lines(trialtimes[c(n,n)],c(0,1),col='red')}
# for (n in 1:length(triggertimes)){lines(triggertimesadjusted[c(n,n)],c(0,1),col='blue',lty=2)}
# for (n in 1:length(triggertimes)){lines(triggertimesadjusted2[c(n,n)],c(0,1),col='orange',lty=2)}
lines(timeseq[timestoplot],0.4+meanEEG[timestoplot],col='green')

EEGtrigs <- round((trialtimes + lagval)*1000)

EEGdata$Trigger[EEGtrigs] <- 1

write.csv(EEGdata,file=paste(EEGdir,participant,'_S',block,'_EEG.csv',sep=''),row.names=FALSE)
gzip(paste(EEGdir,participant,'_S',block,'_EEG.csv',sep=''),destname=paste(EEGdir,participant,'_S',block,'_EEG.gz',sep=''),overwrite=TRUE,remove=TRUE)

}






