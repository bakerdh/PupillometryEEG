Instructions for Docker

Docker is a system that permits preservation of the full computational environment used in an analysis. The enclosed Dockerfile can be used to reconstruct this on your own machine, by following the instructions below. Reading online documentation to familiarise yourself with Docker is also recommended if you have not used it before.

Note that at the time of writing this does not work properly on Apple Silicon devices (M1 and M2 processors built on the ARM architecture) because the rocker project do not yet make frozen builds available for this architecture. However if you edit the first line of the Dockerfile to read:

FROM rocker/rstudio:latest-daily

it will run with the most recent build of RStudio instead of version 4.2.2. This seems likely to work for the forseeable future. However this build may not have all the required packages for building the pdf file.


Instructions:

1. First you need to download and install Docker desktop from: https://www.docker.com/

2. Next download the Dockerfile from the project GitHub repository and place it in a subdirectory of your root user directory called 'docker'

3. Run the docker file by entering the following in a terminal:
docker build -t babymind -f docker/Dockerfile .

4. Once the image has finished building, launch a container in Docker desktop by clicking the play icon

5. Call the container something sensible (i.e. babymind), set the host port to 8787, and add an Environment variable called PASSWORD with the value: pass

6. The container should then appear in the Containers section of Docker desktop, and can be accessed by clicking on the port number (8787:8787), which will launch it in your web browser. Once the container launches, you can login to the RStudio instance with username: rstudio, password: pass

7. Click on the Terminal tab in the RStudio window and enter:
git clone https://github.com/bakerdh/babymind.git
This will download the repository to a folder called Bilingual-switching-ageing

8. You can then open the project and the markdown file, and run the code
