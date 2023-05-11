FROM pytorch/pytorch:1.11.0-cuda11.3-cudnn8-runtime

LABEL authors="FeTS_Admin <admin@fets.ai>"

RUN apt-get update -y

RUN apt-get install wget zip unzip software-properties-common gcc g++ make -y

RUN apt-get update -y && add-apt-repository ppa:deadsnakes/ppa && apt update -y && apt install python3.7 python3.7-venv python3.7-dev python3-setuptools -y

# We will do git pull on the FeTS_Front-End master, because that is the repo using which the base image is made
# We will not do compiles on the PR because the idea is that the Xenial build will check the build status of
# the PR in any case.

ARG VERSION=2.0.0

# download installer
RUN wget https://fets.projects.nitrc.org/FeTS_${VERSION}_Installer.bin && chmod +x FeTS_${VERSION}_Installer.bin

# install FeTS and remove installer
RUN yes yes | ./FeTS_${VERSION}_Installer.bin --target ./FeTS_${VERSION} -- --cudaVersion 11 && rm -rf ./FeTS_${VERSION}_Installer.bin

# set up the docker for GUI
ENV QT_X11_NO_MITSHM=1
ENV QT_GRAPHICSSYSTEM="native"

# define entry point
ENTRYPOINT ["/FeTS_\${VERSION}/bin/install/appdir/usr/bin/BraTSPipeline"]