# Different rules for binocular combination of luminance flicker in cortical and subcortical pathways (code)

These materials are a computationally reproducible version of the paper:

Segala, F.G., Bruno, A., Aung, M.T., Wade, A.R. & Baker, D.H. (2023). Different rules for binocular combination of luminance flicker in cortical and subcortical pathways, *eLife*, https://doi.org/10.7554/eLife.87048

The file manuscript.Rmd is an R markdown file that will perform all analyses and figure creation, and produce a pdf version of the manuscript.

The full repository can be downloaded (cloned), and will automatically download the required packages and data files, depending on the level of analysis requested. If any data files are missing the code will attempt to download them from the OSF repository for this project:
http://doi.org/10.17605/OSF.IO/TBEMA

The 'docker' directory contains a Dockerfile and instructions for making a local computationally reproducible version of the analysis. In addition, the Docker environment is set up to run automatically on a remote server via Github Actions, each time a change is made (i.e. on a 'commit' to the repo). The output document is then posted back to the main repository (manuscript.pdf). If you want to make changes to the analysis and have these build automatically, you can fork the repository into your own account.

![autobuild](https://github.com/bakerdh/PupillometryEEG/workflows/autobuild/badge.svg)