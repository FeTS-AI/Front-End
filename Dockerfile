FROM ghcr.io/fets-ai/fetstool_docker_dependencies:0.0.2.gpu

LABEL authors="FeTS_Admin <admin@fets.ai>"

RUN apt-get update && apt-get update --fix-missing && apt-get install -y libnss3 libnspr4 libxcursor1 libxcursor-dev libasound2 libdbus-1-dev libglfw3-dev libgles2-mesa-dev ffmpeg libsm6 libxext6 python3.8 python3.8-venv python3.8-dev python3-setuptools

ENV PATH=/workspace/CaPTk/bin/qt/5.12.1/bin:/workspace/CaPTk/bin/qt/5.12.1/libexec:$PATH
ENV CMAKE_PREFIX_PATH=/workspace/CaPTk/bin/ITK-build:/workspace/CaPTk/bin/DCMTK-build:/workspace/CaPTk/bin/qt/5.12.1/lib/cmake/Qt5:$CMAKE_PREFIX_PATH

RUN pwd && ls -l && ls -l workspace/

WORKDIR /Front-End

COPY . .

RUN pwd && ls -l

## C++ build
RUN mkdir bin && cd bin && cmake -DCMAKE_INSTALL_PREFIX="./install/appdir/usr" -DITK_DIR="/workspace/CaPTk/bin/ITK-build" -DDCMTK_DIR="/workspace/CaPTk/bin/DCMTK-build" -DBUILD_TESTING=OFF .. && make -j$(nproc) && make install/strip 

## Python package installation
RUN cd bin/install/appdir/usr/bin/ && python3.8 -m venv ./venv && ./venv/bin/pip install --upgrade pip wheel && ./venv/bin/pip install torch==1.13.1+cpu torchvision==0.14.1+cpu torchaudio==0.13.1 --extra-index-url https://download.pytorch.org/whl/cpu && ./venv/bin/pip install -e . && ./venv/bin/pip install setuptools-rust Cython scikit-build scikit-learn openvino-dev==2023.0.1 && ./venv/bin/pip install -e .

### put together a data example that is already aligned and ready to invoke the brain extraction and tumor segmentation

# set up the docker for GUI
ENV LD_LIBRARY_PATH=/CaPTk/bin/qt/5.12.1/lib:$LD_LIBRARY_PATH
ENV PATH=/Front-End/bin/install/appdir/usr/bin/:$PATH
ENV QT_X11_NO_MITSHM=1
ENV QT_GRAPHICSSYSTEM="native"

RUN echo "Env paths\n" && echo $PATH && echo $LD_LIBRARY_PATH

# define entry point
ENTRYPOINT ["/Front-End/bin/install/appdir/usr/bin/FeTS_CLI_Segment"]
