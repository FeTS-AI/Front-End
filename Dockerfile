FROM ghcr.io/fets-ai/fetstool_docker_dependencies:0.0.2.gpu

LABEL authors="FeTS_Admin <admin@fets.ai>"

RUN apt-get update && apt-get update --fix-missing


#general dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    apt-utils \
    sudo \
    libssl-dev \
    make \
    gcc-5 \
    g++-5 
RUN apt-get update && apt-get install -y \
    wget \
    git \
    liblapack-dev \
    unzip \
    tcl \
    tcl-dev 
RUN apt-get update && apt-get install -y \
    tk \
    tk-dev \
    libgl1-mesa-dev \
    libxt-dev \
    libmpc-dev \
    libmpfr-dev 
RUN apt-get update && apt-get install -y \
    libgmp-dev \
    dos2unix \
    doxygen \
    libubsan0 \
    libcilkrts5 

# installing CMake
RUN rm -rf /usr/bin/cmake; \
    wget https://cmake.org/files/v3.12/cmake-3.12.4-Linux-x86_64.sh; \
    mkdir /opt/cmake; \
    sh cmake-3.12.4-Linux-x86_64.sh --prefix=/opt/cmake --skip-license; \
    ln -s /opt/cmake/bin/cmake /usr/bin/cmake; \
    rm -rf cmake-3.12.4-Linux-x86_64.sh

# setting up the build environment
ARG GIT_LFS_SKIP_SMUDGE=1
ARG PKG_FAST_MODE=1
ARG PKG_COPY_QT_LIBS=1
ENV GIT_LFS_SKIP_SMUDGE=$GIT_LFS_SKIP_SMUDGE
ENV PKG_FAST_MODE=$PKG_FAST_MODE
ENV PKG_COPY_QT_LIBS=$PKG_COPY_QT_LIBS

# cloning CaPTk
RUN if [ ! -d "`pwd`/CaPTk" ] ; then git clone "https://github.com/CBICA/CaPTk.git" CaPTk; fi 
RUN cd CaPTk &&  git pull; \
    git submodule update --init && mkdir bin

RUN cd CaPTk/bin && echo "=== Starting CaPTk Superbuild ===" && \
    if [ ! -d "`pwd`/qt" ] ; then wget https://github.com/CBICA/CaPTk/raw/master/binaries/qt_5.12.1/linux.zip -O qt.zip; fi ; \
    cmake -DCMAKE_INSTALL_PREFIX=./install_libs -DQT_DOWNLOAD_FORCE=OFF -Wno-dev .. && make -j$(nproc) && rm -rf qt.zip && cd .. && mkdir Front-End

RUN pwd && ls -l

WORKDIR /Front-End

COPY . .

RUN pwd && ls -l

## C++ build
RUN mkdir bin && cd bin && cmake -DCMAKE_INSTALL_PREFIX="./install/appdir/usr" -DITK_DIR="/CaPTk/bin/ITK-build" -DDCMTK_DIR="/CaPTk/bin/DCMTK-build" -DBUILD_TESTING=OFF .. && make -j$(nproc) && make install/strip 

## Python package installation
RUN cd bin/install/appdir/usr/bin/ && python3.8 -m venv ./venv && ./venv/bin/pip install --upgrade pip wheel && ./venv/bin/pip install torch==1.13.1+cpu torchvision==0.14.1+cpu torchaudio==0.13.1 --extra-index-url https://download.pytorch.org/whl/cpu && ./venv/bin/pip install -e . && ./venv/bin/pip install setuptools-rust Cython scikit-build scikit-learn openvino-dev==2023.0.1 && ./venv/bin/pip install -e .

# set up the docker for GUI
ENV LD_LIBRARY_PATH=/CaPTk/bin/qt/5.12.1/lib:$LD_LIBRARY_PATH
ENV PATH=/Front-End/bin/install/appdir/usr/bin/:$PATH
ENV QT_X11_NO_MITSHM=1
ENV QT_GRAPHICSSYSTEM="native"

RUN echo "Env paths\n" && echo $PATH && echo $LD_LIBRARY_PATH

# define entry point
ENTRYPOINT ["/Front-End/bin/install/appdir/usr/bin/FeTS_CLI_Segment"]
