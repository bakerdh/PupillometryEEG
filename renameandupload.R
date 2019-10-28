# script to collate the various files and upload to OSF
# checks which files already exist and which need to be uploaded
# important to run the matlab script to convert EEG data first

participant <- 'P310'
sourcedir <- '~/Desktop/Pupildata/'

localdir <- 'local/'    # all files are stored in the project directory /local/ which git is told to ignore
if (!file.exists(localdir)){dir.create(localdir)}   # create a local directory to store data and outputs
EEGdir <- 'local/EEG/'
if (!file.exists(EEGdir)){dir.create(EEGdir)}   # create a local directory to store EEG data
PPdir <- 'local/Pupil/'
if (!file.exists(PPdir)){dir.create(PPdir)}   # create a local directory to store pupil data
Pydir <- 'local/Psychopy/'
if (!file.exists(Pydir)){dir.create(Pydir)}   # create a local directory to store Psychopy output
figdir <- 'local/Figures/'
if (!file.exists(figdir)){dir.create(figdir)}   # create a local directory to store figures
datadir <- 'local/Data/'
if (!file.exists(datadir)){dir.create(datadir)}   # create a local directory to store processed data

packagelist <- c('remotes','tictoc','utils')
missingpackages <- packagelist[!packagelist %in% installed.packages()[,1]]
if (length(missingpackages)>0){install.packages(missingpackages)}
if (packagelist[!'osfr' %in% installed.packages()[,1]]){remotes::install_github("centerforopenscience/osfr")}
packagelist <- c(packagelist,'osfr')
toinstall <- packagelist[which(!packagelist %in% (.packages()))]
invisible(lapply(toinstall,library,character.only=TRUE))

tic()

my_wd <- getwd() # save your current working directory path

# move EEG files to the correct place and rename
d <- dir(paste(sourcedir,'EEG/',participant,sep=''),pattern='*.csv.gz')
if (length(d)){
for (n in 1:length(d)){
  file.copy(paste(sourcedir,'EEG/',participant,'/',d[n],sep=''),paste(sourcedir,'EEG/',participant,'_S',n,'_EEG.gz',sep=''))
  file.remove(paste(sourcedir,'EEG/',participant,'/',d[n],sep=''))
}}
d <- dir(paste(sourcedir,'EEG/',participant,sep=''),full.names=TRUE)
if (!file.exists(paste(sourcedir,'EEG/',participant,'.zip',sep=''))){zip(paste(sourcedir,'EEG/',participant,'.zip',sep=''),files=d)}   # zip the EEG directory

# copy and rename info and pupil position files
for (n in 1:3){
  if (!file.exists(paste(sourcedir,'Eyetracking/',participant,'/',participant,'_S',n,'_info.csv',sep=''))){
  file.copy(paste(sourcedir,'Eyetracking/',participant,'/00',n-1,'/info.csv',sep=''),paste(sourcedir,'Eyetracking/',participant,'/',participant,'_S',n,'_info.csv',sep=''))}
  if (!file.exists(paste(sourcedir,'Eyetracking/',participant,'/',participant,'_S',n,'_pupil_positions.csv',sep=''))){
  file.copy(paste(sourcedir,'Eyetracking/',participant,'/00',n-1,'/exports/000/pupil_positions.csv',sep=''),paste(sourcedir,'Eyetracking/',participant,'/',participant,'_S',n,'_pupil_positions.csv',sep=''))}
}

# zip up the eyetracking directories
for (n in 1:3){
  if (!file.exists(paste(sourcedir,'Eyetracking/',participant,'/',participant,'_S',n,'.zip',sep=''))){
  destpath <- paste(sourcedir,'Eyetracking/',participant,'/00',n-1,sep='')
  setwd(destpath)
  zip(paste(sourcedir,'Eyetracking/',participant,'/00',n-1,'.zip',sep=''),files=list.files(destpath))
  file.rename(paste(sourcedir,'Eyetracking/',participant,'/00',n-1,'.zip',sep=''),paste(sourcedir,'Eyetracking/',participant,'/',participant,'_S',n,'.zip',sep=''))
}}

setwd(my_wd) # reset working directory path

osf_auth(token = '3dEYuhZNmwWbG3xuhRIiVRo0T2oniOkEP3Ip8i1LPG4PNjeTRln54eGNAG7oyTT9xozWwJ')
osfproject <- osf_retrieve_node("tbema")
componentlist <- osf_ls_nodes(osfproject)
EEGID <- match('Raw EEG data',as.character(unlist(componentlist[,1])))
EEGtoken <- componentlist[EEGID,2]
EEGfiles <- osf_ls_files(EEGtoken,n_max=300)
PPID <- match('Raw pupillometry data',as.character(unlist(componentlist[,1])))
PPtoken <- componentlist[PPID,2]
PPfiles <- osf_ls_files(PPtoken,n_max=300)
PyID <- match('Psychopy logs',as.character(unlist(componentlist[,1])))
Pytoken <- componentlist[PyID,2]
Pyfiles <- osf_ls_files(Pytoken,n_max=300)

componentlist2 <- osf_ls_nodes(EEGtoken)
EEGID2 <- match('ASAfiles',as.character(unlist(componentlist2[,1])))
EEGtoken2 <- componentlist2[EEGID2,2]
EEGfiles2 <- osf_ls_files(EEGtoken2,n_max=300)

componentlist3 <- osf_ls_nodes(PPtoken)
PPID2 <- match('Videofiles',as.character(unlist(componentlist3[,1])))
PPtoken2 <- componentlist3[PPID2,2]
PPfiles2 <- osf_ls_files(PPtoken2,n_max=300)

pfilelist <- dir(paste(sourcedir,'Psychopy',sep=''),pattern='*.csv')
indices <- NULL
counter <- 0
for (n in 1:length(pfilelist)){
if (pmatch(participant,pfilelist[n],nomatch=0)){counter <- counter + 1
indices[counter] <- n}
}
pfilelist <- pfilelist[indices]
for (n in 1:length(pfilelist)){
  if (!pmatch(pfilelist[n],Pyfiles$name,nomatch=0)){
    osf_upload(Pytoken,paste(sourcedir,'Psychopy/',pfilelist[n],sep=''))
  }
}

for (n in 1:3){
  if (file.exists(paste(sourcedir,'EEG/',participant,'_S',n,'_EEG.gz',sep=''))){
  if (!pmatch(paste(participant,'_S',n,'_EEG.gz',sep=''),EEGfiles$name,nomatch=0)){
    osf_upload(EEGtoken,paste(sourcedir,'EEG/',participant,'_S',n,'_EEG.gz',sep=''))
  }}
}

for (n in 1:3){
  if (file.exists(paste(sourcedir,'Eyetracking/',participant,'/',participant,'_S',n,'_info.csv',sep=''))){
  if (!pmatch(paste(participant,'_S',n,'_info.csv',sep=''),PPfiles$name,nomatch=0)){
    osf_upload(PPtoken,paste(sourcedir,'Eyetracking/',participant,'/',participant,'_S',n,'_info.csv',sep=''))
  }}
  if (file.exists(paste(sourcedir,'Eyetracking/',participant,'/',participant,'_S',n,'_pupil_positions.csv',sep=''))){
  if (!pmatch(paste(participant,'_S',n,'_pupil_positions.csv',sep=''),PPfiles$name,nomatch=0)){
    osf_upload(PPtoken,paste(sourcedir,'Eyetracking/',participant,'/',participant,'_S',n,'_pupil_positions.csv',sep=''))
  }}  
}

if (file.exists(paste(sourcedir,'EEG/',participant,'.zip',sep=''))){
  if (!pmatch(paste(participant,'.zip',sep=''),EEGfiles2$name,nomatch=0)){
    osf_upload(EEGtoken2,paste(sourcedir,'EEG/',participant,'.zip',sep=''))
  }}

for (n in 1:3){
  if (file.exists(paste(sourcedir,'Eyetracking/',participant,'/',participant,'_S',n,'.zip',sep=''))){
  if (!pmatch(paste(participant,'_S',n,'.zip',sep=''),PPfiles2$name,nomatch=0)){
    osf_upload(PPtoken2,paste(sourcedir,'Eyetracking/',participant,'/',participant,'_S',n,'.zip',sep=''))
  }}
}
