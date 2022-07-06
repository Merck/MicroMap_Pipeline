# Set the base image
FROM r-base:4.0.2

# File Author / Maintainer
MAINTAINER Cory White <cory.white@merck.com>
RUN mkdir /output/
WORKDIR /Prox_Labeling

COPY ./bin/ /Prox_Labeling/bin/
COPY ./data/ /Prox_Labeling/data/
COPY ./README.md /Prox_Labeling
COPY ./docker_files/ /Prox_Labeling/docker_files/

WORKDIR docker_files/

#RUN wget https://www.imagemagick.org/download/ImageMagick.tar.gz
RUN tar -zxvf ImageMagick.tar.gz
WORKDIR /Prox_Labeling/docker_files/ImageMagick-7.1.0-39/
RUN ./configure
RUN make
RUN make install
RUN ldconfig /usr/local/lib
WORKDIR ../

RUN tar -zxvf ghostscript-9.56.1.tar.xz
WORKDIR /Prox_Labeling/docker_files/ghostscript-9.56.1/
RUN ./configure
RUN make
RUN make install

WORKDIR ../

RUN Rscript /Prox_Labeling/docker_files/install.r

WORKDIR ../

RUN rm -r docker_files/

WORKDIR /Prox_Labeling/bin/

#ENTRYPOINT /Prox_Labeling/bin/Run_All.sh

LABEL	dock.img.name="MicroMap_Proximity_Labeling"   \
        dock.img.description="Image for Merck Proximity Labeling Pipeline"   \
        \
        dock.maintainer.name="Cory White"                                     \
        dock.maintainer.isid="<ISID>"                                          \
        dock.maintainer.email="cory.white@merck.com"                         \
        dock.maintainer.division="ESC"                            \
        dock.maintainer.team="ESC - Systems Biology"                                    \
        \
        dock.docker.run='docker run --rm -it /bin/bash'                          \
        dock.docker.run.runs-forever='false'                                     \
        \
        dock.schema-version="0.1" \