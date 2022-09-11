
# packagelist <- c('remotes','tictoc','utils')
# missingpackages <- packagelist[!packagelist %in% installed.packages()[,1]]
# if (length(missingpackages)>0){install.packages(missingpackages)}
# if (packagelist[!'osfr' %in% installed.packages()[,1]]){remotes::install_github("centerforopenscience/osfr")}
# packagelist <- c(packagelist,'osfr')
# toinstall <- packagelist[which(!packagelist %in% (.packages()))]
# invisible(lapply(toinstall,library,character.only=TRUE))



d <- dir('~/Desktop/P101/')
p <- getwd()
setwd('~/Desktop/P101/')
tar(tarfile='P101.tar',files=d,compression='gzip')
setwd(p)

osf_auth(token = '3dEYuhZNmwWbG3xuhRIiVRo0T2oniOkEP3Ip8i1LPG4PNjeTRln54eGNAG7oyTT9xozWwJ')
osfproject <- osf_retrieve_node("x8u4v")
componentlist <- osf_ls_nodes(osfproject)
datafilesonOSF <- osf_ls_files(osfproject,n_max=300)

if (!pmatch('P101.tar',datafilesonOSF$name,nomatch=0)){
  osf_upload(osfproject,'~/Desktop/P101/P101.tar',progress=TRUE)
}

