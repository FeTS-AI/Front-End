FROM ghcr.io/fets-ai/fetstool_docker_dependencies AS fets_base

LABEL authors="FeTS_Admin <admin@fets.ai>"

RUN apt-get update && apt-get update --fix-missing && apt-get install -y libnss3 libnspr4 libxcursor1 libxcursor-dev libasound2 libdbus-1-dev libglfw3-dev libgles2-mesa-dev ffmpeg libsm6 libxext6 python3.8 python3.8-venv python3.8-dev python3-setuptools

ENV PATH=/CaPTk/bin/qt/5.12.1/bin:/CaPTk/bin/qt/5.12.1/libexec:$PATH
ENV CMAKE_PREFIX_PATH=/CaPTk/bin/ITK-build:/CaPTk/bin/DCMTK-build:/CaPTk/bin/qt/5.12.1/lib/cmake/Qt5:$CMAKE_PREFIX_PATH

RUN pwd && ls -l

WORKDIR /Front-End

COPY src src

COPY CMakeLists.txt README.txt LICENSE .

COPY cmake_modules cmake_modules

COPY data data

COPY docs_sources docs_sources

RUN pwd && ls -l

## C++ build
RUN mkdir bin && cd bin && cmake -DCMAKE_INSTALL_PREFIX="./install/appdir/usr" -DITK_DIR="/CaPTk/bin/ITK-build" -DDCMTK_DIR="/CaPTk/bin/DCMTK-build" -DBUILD_TESTING=OFF .. && make -j$(nproc) && make install/strip 

## Python package installation
RUN cd bin/install/appdir/usr/bin/ && python3.8 -m venv ./venv && ./venv/bin/pip install --upgrade pip wheel && ./venv/bin/pip install torch==1.13.1+cpu torchvision==0.14.1+cpu torchaudio==0.13.1 --extra-index-url https://download.pytorch.org/whl/cpu && ./venv/bin/pip install -e . && ./venv/bin/pip install setuptools-rust Cython scikit-build scikit-learn openvino-dev==2022.1.0 && ./venv/bin/pip install -e .

# set up the docker for GUI
ENV LD_LIBRARY_PATH=/CaPTk/bin/qt/5.12.1/lib:$LD_LIBRARY_PATH
ENV PATH=/Front-End/bin/install/appdir/usr/bin/:$PATH
ENV QT_X11_NO_MITSHM=1
ENV QT_GRAPHICSSYSTEM="native"

RUN echo "Env paths\n" && echo $PATH && echo $LD_LIBRARY_PATH

# define entry point
ENTRYPOINT ["/Front-End/bin/install/appdir/usr/bin/venv/bin/python", "/Front-End/bin/install/appdir/usr/bin/PrepareDataset.py"]

FROM fets_base AS data_prep

RUN find /Front-End/bin/install/appdir/usr/bin -type f \( -perm -u=x -o -type l \) -exec cp -P {} /usr/bin \;

WORKDIR /

RUN apt-get install software-properties-common curl -y && \
    add-apt-repository ppa:deadsnakes/ppa -y && apt-get update && \
    apt-get install python3.8 python3.8-distutils -y && \
    apt-get remove --purge python3.6 -y && \
    apt autoremove -y && \
    apt-get install python3.8-distutils -y && \
    rm -fr /usr/bin/python /usr/bin/python3 /usr/bin/pip /usr/bin/pip3 && \
    ln -s /usr/bin/python3.8 /usr/bin/python && ln -s /usr/bin/python3.8 /usr/bin/python3 && \
    ln -s /usr/bin/pip3.8 /usr/bin/pip && ln -s /usr/bin/pip3.8 /usr/bin/pip3

RUN curl -fSsL -O https://bootstrap.pypa.io/get-pip.py -o get-pip.py && \
    python3.8 get-pip.py && \
    rm get-pip.py

COPY ./mlcubes/data_preparation/project/requirements.txt /project/requirements.txt 

RUN pip install --upgrade pip

RUN cd /Front-End/bin/install/appdir/usr/bin/ && pip install --upgrade pip wheel && pip install torch==1.13.1+cpu torchvision==0.14.1+cpu torchaudio==0.13.1 --extra-index-url https://download.pytorch.org/whl/cpu && pip install -e . && pip install setuptools-rust Cython scikit-build scikit-learn openvino-dev==2022.1.0 && pip install -e .

RUN pip install -r /project/requirements.txt

ENV LANG C.UTF-8

COPY ./mlcubes/data_preparation/project /project

RUN cp /Front-End/src/applications/*.py /project/stages/

ENTRYPOINT ["python", "/project/mlcube.py"]
