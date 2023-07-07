FROM ubuntu:18.04
MAINTAINER MLPerf MLBox Working Group

RUN apt-get update && \
	apt-get install -y --no-install-recommends \
	software-properties-common \
	cmake \
	python3-dev \
	build-essential \
	qtbase5-dev \
	libqt5*-dev \
	qtwayland5 \
	qtwebengine5-dev \
	qml-module-qtwebview \
	git \
	libxt-dev  \
	doxygen \
	curl && \
	rm -rf /var/lib/apt/lists/*

# TODO: Point to the main repo/branch once this is merged
RUN git clone --branch fets_2.0 https://github.com/FeTS-AI/Front-End

WORKDIR /Front-End

# Build binaries
RUN mkdir bin

WORKDIR /Front-End/bin

RUN cmake ..

RUN make -j4

ENV PATH=/opt/qt/5.11.2/gcc_64/bin:/opt/qt/5.11.2/gcc_64/libexec:$PATH

ENV CMAKE_PREFIX_PATH=/work/CaPTk/bin/ITK-build:/work/CaPTk/bin/DCMTK-build:/opt/qt/5.11.2/gcc_64/lib/cmake/Qt5:$CMAKE_PREFIX_PATH

RUN cmake -DCMAKE_INSTALL_PREFIX="./install/appdir/usr" -DITK_DIR="/work/CaPTk/bin/ITK-build" -DDCMTK_DIR="/work/CaPTk/bin/DCMTK-build" -DBUILD_TESTING=OFF ..

RUN make -j4

RUN make install/strip

RUN find /Front-End/bin/install/appdir/usr/bin -type f \( -perm -u=x -o -type l \) -exec cp -P {} /usr/bin \;

WORKDIR /

RUN add-apt-repository ppa:deadsnakes/ppa -y && apt-get update

RUN apt-get install python3 -y

RUN curl -fSsL -O https://bootstrap.pypa.io/pip/3.6/get-pip.py && \
	python3 get-pip.py && \
	rm get-pip.py

COPY ./requirements.txt /project/requirements.txt 

RUN pip3 install --upgrade pip

RUN pip3 install -r /project/requirements.txt

ENV LANG C.UTF-8

COPY . /project

#############################################################
# TODO: This should be removed once everything gets merged
RUN rm -rf /Front-End

RUN git clone --branch reformat-prepare-script https://github.com/aristizabal95/FeTS-AI-Front-End
##############################################################

COPY /Front-End/src/applications/PrepareDataset.py /project/stages/PrepareDataset.py

COPY /Front-End/src/applications/CreateCSVForDICOMs.py /project/stages/CreateCSVForDICOMs.py

ENV QT_X11_NO_MITSHM=1
ENV QT_GRAPHICSSYSTEM="native"

RUN rm -rf /Front-End

ENTRYPOINT ["python3", "/project/mlcube.py"]