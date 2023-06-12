# packagelist <- c('remotes','tictoc','utils')
# missingpackages <- packagelist[!packagelist %in% installed.packages()[,1]]
# if (length(missingpackages)>0){install.packages(missingpackages)}
# if (packagelist[!'osfr' %in% installed.packages()[,1]]){remotes::install_github("centerforopenscience/osfr")}
# packagelist <- c(packagelist,'osfr')
# toinstall <- packagelist[which(!packagelist %in% (.packages()))]
# invisible(lapply(toinstall,library,character.only=TRUE))

library(osfr)

osf_auth(token = 'PRIVATETOKENGOESHERE')
osfproject <- osf_retrieve_node("x8u4v")
componentlist <- osf_ls_nodes(osfproject)

p <- getwd()

for (s in 1:30){
  
  subj <- paste0('P',s+100)
  
  subjdir <- paste0('/Users/danbaker/Desktop/Pupildata/',subj)
  if (!file.exists(subjdir)){dir.create(subjdir)}
  
  eegfiles <- dir('/Volumes/DANLAB/Eyetracking/Federico/Exp1/EEG/Zipped',pattern=paste0(subj,'_*'),full.names=TRUE)
  file.copy(eegfiles,subjdir)
  
  setwd(subjdir)
  
  d <- dir(subjdir)
  for (n in 1:length(d)){file.rename(d[n],paste0(unlist(strsplit(d[n],'[.]'))[1],'.csv.gz'))}
  
  PPfiles <- dir('/Volumes/DANLAB/Eyetracking/Federico/Exp1/Psychopy/CSV/',pattern=paste0(subj,'_*'),full.names=TRUE)
  file.copy(PPfiles,subjdir)
  pupilfiles <- dir('/Volumes/DANLAB/Eyetracking/Federico/Exp1/Eyetracking/CSV/',pattern=paste0(subj,'_*'),full.names=TRUE)
  file.copy(pupilfiles,subjdir)
  
  d <- dir(subjdir)
  
  tar(tarfile=paste0(subj,'.tar'),files=d,compression='gzip')
  setwd(p)
  
  datafilesonOSF <- osf_ls_files(osfproject,n_max=300)
  
  if (!pmatch(paste0(subj,'.tar'),datafilesonOSF$name,nomatch=0)){
    osf_upload(osfproject,paste0(subjdir,'/',subj,'.tar'),progress=TRUE)
  }
  
  unlink(subjdir,recursive=TRUE)
  
}




osfproject <- osf_retrieve_node("mvpyz")
componentlist <- osf_ls_nodes(osfproject)

p <- getwd()

for (s in 1:12){
  
  subj <- paste0('P',s+150)
  
  subjdir <- paste0('local/TFdata/',subj)

  setwd(subjdir)
  
  d <- dir()
  
  tar(tarfile=paste0(subj,'.tar'),files=d,compression='gzip')
  setwd(p)
  
  datafilesonOSF <- osf_ls_files(osfproject,n_max=300)
  
  if (!pmatch(paste0(subj,'.tar'),datafilesonOSF$name,nomatch=0)){
    osf_upload(osfproject,paste0(subjdir,'/',subj,'.tar'),progress=TRUE)
  }
  
  # unlink(subjdir,recursive=TRUE)
  
}


osfproject <- osf_retrieve_node("sjymf")
componentlist <- osf_ls_nodes(osfproject)

p <- getwd()

for (s in 1:12){
  
  subj <- paste0('P',s+200)
  
  subjdir <- paste0('/Users/db900/Documents/git/PupillometryEEG/local/rawdata/',subj)
  if (!file.exists(subjdir)){dir.create(subjdir)}
  
  setwd(subjdir)
  
  d <- dir(subjdir)
  
  tar(tarfile=paste0(subj,'.tar'),files=d,compression='gzip')
  setwd(p)
  
  datafilesonOSF <- osf_ls_files(osfproject,n_max=300)
  
  if (!pmatch(paste0(subj,'.tar'),datafilesonOSF$name,nomatch=0)){
    osf_upload(osfproject,paste0(subjdir,'/',subj,'.tar'),progress=TRUE)
  }
  
  # unlink(subjdir,recursive=TRUE)
  
}
