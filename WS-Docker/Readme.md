# Introduction
In this workshop we introduce Docker (<https://www.docker.com>).
Docker is a tool to easily create, deploy and share linux-based images for virtual machine.
All documentation and manual can be found on <https://docs.docker.com>.
This tutorial was written on a macOS environment.


# Installation
1. Download and install docker for your platform (<http://docker.com>)
2. Launch the installed docker application
3. Open settings (or Preferences). In panel *"Advanced"*, set the desired number of CPUs and memory allocated to the VMs


# Quick start
1. In a terminal, run  
   `docker run -it ubuntu bash`  
   *This download the latest minimal ubuntu image from "docker hub" (a repository of VM images). It then creates a new container from this image, and run an interactive shell in it. You can browse available images on <https://hub.docker.com> or <https://store.docker.com>.*
2. In another terminal, run  
   `docker run -d -e USER=username -e PASSWORD=passwd -p 8080:8787 bioconductor/release_core`  
   *This download the latest official bioconductor image from docker hub, and run the Rstudio server installed in it. The server listen on port 8787 and is remapped to port 8080 on your host machine. You can access the service with a web-browser at the address <http://localhost:8080> and login with the given username and password.*
3. You need an HTML web server ? Just run  
   `docker run -d -p 80:80 -v $(pwd):/usr/share/nginx/html:ro -d nginx`  
   *This maps the current working directory $(pwd) of the host machine onto /usr/share/nginx/html/ in the guest. It maps as well port 80 of the guest onto port 80 of the host and run nginx web server. Just put an html page in the working directory of the host machine, and open your web-browser to the address <http://localhost>*


# Quick docker reference manual
Detailed documentation and manual can be found on <https://docs.docker.com>. Here are useful commands.

1. List running containers: `docker ps`
2. Stop a container: `docker stop <container-name>`
3. List all containers (running and stopped containers): `docker ps -a`
4. Restart a container: `docker restart <container-name>`  
   *This will restore all modified files in the container*
5. Delete a container: `docker rm <container-name>`  
   *Warning, this will erase all created files or any change done in the container*
6. List downloaded images: `docker images`
7. Save modification done in a container as a new image `docker commit <container-name>`
8. Run a command in a running container: `docker exec <container-name> <command>`  
   e.g. `docker exec -it <container-name> bash`
9. Give a name to your image: `docker tag <name>`




# Build a custom image
In this example, we show how to create a custom image with SAMTOOLS and BWA installed inside.

1. Create an empty directory named `ngs/`
2. Create inside a text file named `Dockerfile`, with the content given below.
3. run `docker build -t ngs ./`  
   *Staring from an ubuntu image, it runs the successive commands as if you were typing them into a shell. After each command an intermediate image is created to avoid restarting compilation from scratch when a change in the Dockerfile is introduced. The full documentation for Dockerfile is here <https://docs.docker.com/engine/reference/builder/>*. 

```
FROM ubuntu

RUN apt-get update && apt-get install -y \
    git \
    python \
    make \
    curl \
    bzip2 \
    g++ \
    libz-dev \
    wget \
    libncurses-dev \
    unzip \
    ftp \
    vim


#-#-#-#-#-#-#-#-#-#-#-#-#
# Install SAMTOOLS
#-#-#-#-#-#-#-#-#-#-#-#-#
RUN curl -kL https://github.com/samtools/htslib/releases/download/1.3.2/htslib-1.3.2.tar.bz2 | tar -C /tmp -jxf - && \
    cd /tmp/htslib-1.3.2 && make && make install && \
    rm -rf /tmp/htslib-1.3.2
RUN curl -kL https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2 | tar -C /tmp -jxf - && \
    cd /tmp/samtools-1.3.1 && make && make install && \
    rm -rf /tmp/samtools-1.3.1
RUN curl -kL https://github.com/samtools/bcftools/releases/download/1.3.1/bcftools-1.3.1.tar.bz2 | tar -C /tmp -jxf - && \
    cd /tmp/bcftools-1.3.1 && make && make install && \
    rm -rf /tmp/bcftools-1.3.1



#-#-#-#-#-#-#-#-#-#-#-#-#
# Install BWA
#-#-#-#-#-#-#-#-#-#-#-#-#
RUN curl -kL http://netix.dl.sourceforge.net/project/bio-bwa/bwa-0.7.15.tar.bz2 | tar -C /tmp -jxf - && \
    cd /tmp/bwa-0.7.15 && make && find /tmp/bwa-0.7.15/ -type f -executable -exec mv '{}' /usr/local/bin/ ';' && \
    rm -rf /tmp/bwa-0.7.15/

ENTRYPOINT ["/bin/bash"]
```

## Test your image
Once the image is build you can test it:
```
docker run -it ngs
>bwa
>samtools
```

## Update Dockerfile
Appending command to your Dockerfile doesn't require complete rebuilding. 
For example, add the following command in your Dockerfile to install STAR in your image.
```
RUN curl -kL https://github.com/alexdobin/STAR/archive/2.5.2b.tar.gz | tar -C /tmp -zxf - && \ 
    mv /tmp/STAR-2.5.2b/bin/Linux_x86_64_static/* /usr/local/bin/ && \
    rm -rf /tmp/STAR-2.5.2b/
```
The building process, that have cached intermediate images, will only redo the necessary steps:  
`docker build -t ngs ./`


## Publish your image on hub.docker.com
1. Create an account on docker-hub <https://hub.docker.com>
2. Run `docker login`
3. Run `docker push`
4. Anyone can now run your image by typing:
   `docker run -it login/ngs`

## Automatic building
Alternatively, you have also the possibility to create a GitHub repository containing your Dockerfile, and setup your docker repository to automatically build your image when you push changes to GitHub. In this case the building process is done on docker servers instead of your local machine, but your are limited in compilation time.


# docker-compose

In docker philosophy, an image provides a service and you compose several services to implement your App. 
A typical configuration consists in a service providing the database (e.g. mongoDB), and web-servers (e.g. nginx) implementing the web-pages that query the database. The database-server and the web-server are on the same local (virtual) network, and only the web-service is exposed to outside.

Here, we will see how we can compose the Rstudio-server we saw in [quick start](#quick-start) with an HTTPS proxy to implement a security layer. Indeed, the *bioconductor/release_core* image we saw [earlier](#quick-start) contains the Open Source Edition of RStudio-server. This edition is free, but the traffic over HTTP is not encrypted i.e. the passwords and the commands you send in your browser are visible to a sniffer and can be intercepted or modified. For security reason it is not recommended to use it over internet, its usage should be restricted to a local network. Encrypted traffic over HTTPS, is however possible with the Professional Edition (see <https://www.rstudio.com/products/rstudio/#Server>).

Using docker, we can also easily compose the Rstudio-server service with an HTTPS-proxy to implement the missing security layer. The Rstudio-server listen for HTTP traffic on port 8787. The HTTPS-proxy will be in charge to forward and translate the encrypted HTTPS-traffic it listen on port 443 to port 8787 of Rstudio-server. The two services will be put on the same virtual network, but only the encrypted traffic of HTTPS-proxy (port 443) will be exposed to outside.


To implement this App, create into an empty directory the file `docker-compose.yml` with the following content:
```
version: '2'
services:
  rstudio:
    image: bioconductor/release_core
    environment:
     - USER=rstudio
     - PASSWORD=rstudio
    volumes:
     - ./:/export

  https:
   image: yajo/https-proxy
   environment:
     - PORT=8787
   ports:
     - "443:443"
   links:
     - rstudio:www
```

Run your App: `docker-compose up`

You can now access the App in your web-browser at this address <https://localhost> or anywhere from the web using the IP of your host. Note that we have mapped the working directory of the host to `/export` on the guest. Then when working within the App, files located into /export will be shared between host and guest and will persist even if the container is deleted.

The detailed documentation for `docker-compose` is here <https://docs.docker.com/compose/>.

> **Note:** even when you have a single service to compose, a `docker-compose.yml` is very useful to store parameters for the `docker` command.




