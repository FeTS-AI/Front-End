FROM ghcr.io/fets-ai/fetstool_docker_dependencies

LABEL authors="FeTS_Admin <admin@fets.ai>"

RUN apt-get update && apt-get update --fix-missing && apt-get install -y libnss3 libnspr4 libxcursor1 libxcursor-dev libasound2 libdbus-1-dev libglfw3-dev libgles2-mesa-dev ffmpeg libsm6 libxext6 python3.8 python3.8-venv python3.8-dev python3-setuptools

ENV PATH=/CaPTk/bin/qt/5.12.1/bin:/CaPTk/bin/qt/5.12.1/libexec:$PATH
ENV CMAKE_PREFIX_PATH=/CaPTk/bin/ITK-build:/CaPTk/bin/DCMTK-build:/CaPTk/bin/qt/5.12.1/lib/cmake/Qt5:$CMAKE_PREFIX_PATH

RUN pwd && ls -l

WORKDIR /Front-End

COPY . .

RUN pwd && ls -l

RUN mkdir bin && cd bin && cmake -DCMAKE_INSTALL_PREFIX="./install/appdir/usr" -DITK_DIR="/CaPTk/bin/ITK-build" -DDCMTK_DIR="/CaPTk/bin/DCMTK-build" -DBUILD_TESTING=OFF .. && make -j$(nproc) && make install/strip && cd ./install/appdir/usr/bin/ && python3.8 -m venv ./venv && ./venv/bin/pip install --upgrade pip setuptools wheel setuptools-rust Cython scikit-build scikit-learn openvino-dev==2022.1.0 && ./venv/bin/pip install torch==1.13.1+cpu torchvision==0.14.1+cpu torchaudio==0.13.1 --extra-index-url https://download.pytorch.org/whl/cpu && ./venv/bin/pip install -e .

# set up the docker for GUI
ENV LD_LIBRARY_PATH=/CaPTk/bin/qt/5.12.1/lib:$LD_LIBRARY_PATH
ENV PATH=/Front-End/bin/install/appdir/usr/bin/:$PATH
ENV QT_X11_NO_MITSHM=1
ENV QT_GRAPHICSSYSTEM="native"

RUN echo "Env paths\n" && echo $PATH && echo $LD_LIBRARY_PATH

# define entry point
ENTRYPOINT ["/Front-End/bin/install/appdir/usr/bin/BraTSPipeline"]
