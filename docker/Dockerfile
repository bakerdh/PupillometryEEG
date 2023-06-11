FROM rocker/verse:4.2.2

RUN apt-get update -qq
RUN apt-get -y install git

RUN install2.r –error \
--deps TRUE \
rmarkdown \
bookdown \
knitr \
remotes \
tictoc \
R.matlab \
grImport \
png \
tiff \
pals \
ez \
gtools \
signal \
boot \
quickpsy \
rstan \
coda \
osfr

