FROM pytorch/pytorch:1.11.0-cuda11.3-cudnn8-runtime

LABEL authors="FeTS_Admin <admin@fets.ai>"

RUN apt-get update -y

RUN apt-get install wget zip unzip yes -y

RUN add-apt-repository ppa:deadsnakes/ppa && apt update -y && apt install python3.6 python3.6-venv python3.6-dev python3-setuptools -y

# We will do git pull on the CBICA/CaPTk master, because that is the repo using which the base image is made
# We will not do compiles on the PR because the idea is that the Xenial build will check the build status of
# the PR in any case.

ARG VERSION=0.0.9

# download installer
RUN wget https://fets.projects.nitrc.org/FeTS_${VERSION}_Installer.bin && chmod +x FeTS_${VERSION}_Installer.bin

# install FeTS and remove installer
RUN yes yes | ./FeTS_${VERSION}_Installer.bin --target ./FeTS_${VERSION} && rm -rf ./FeTS_${VERSION}_Installer.bin

# set up environment and install correct version of pytorch
RUN cd FeTS_${VERSION}/squashfs-root/usr/bin/OpenFederatedLearning && \
    rm -rf ./venv && python3.6 -m venv ./venv && \
    ./venv/bin/pip install torch==1.7.1+cu110 torchvision==0.8.2+cu110 torchaudio==0.7.2 -f https://download.pytorch.org/whl/torch_stable.html && \
    make install_openfl && \
    ./venv/bin/pip install opencv-python==4.2.0.34 && \
    make install_fets && \
    ./venv/bin/pip install -e ./submodules/fets_ai/Algorithms/GANDLF && \
    cd ../LabelFusion && \
    rm -rf venv && python3.6 -m venv ./venv && \
    ./venv/bin/pip install -e .

# set up the docker for GUI
ENV QT_X11_NO_MITSHM=1
ENV QT_GRAPHICSSYSTEM="native"

# define entry point
ENTRYPOINT ["/FeTS_\${VERSION}/bin/install/appdir/usr/bin/FeTS_CLI"]