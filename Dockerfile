FROM ghcr.io/fets-ai/fetstool_docker_dependencies

LABEL authors="FeTS_Admin <admin@fets.ai>"

RUN apt-get update && apt-get update --fix-missing && apt-get install -y libnss3 libnspr4 libxcursor1 libxcursor-dev libasound2 libdbus-1-dev libglfw3-dev libgles2-mesa-dev

ENV PATH=/CaPTk/bin/qt/5.12.1/bin:/CaPTk/bin/qt/5.12.1/libexec:$PATH
ENV CMAKE_PREFIX_PATH=/CaPTk/bin/ITK-build:/CaPTk/bin/DCMTK-build:/CaPTk/bin/qt/5.12.1/lib/cmake/Qt5:$CMAKE_PREFIX_PATH

RUN pwd && ls -l

WORKDIR /Front-End

COPY . .

RUN pwd && ls -l

RUN mkdir bin && cd bin && cmake -DCMAKE_INSTALL_PREFIX="./install/appdir/usr" -DITK_DIR="/CaPTk/bin/ITK-build" -DDCMTK_DIR="/CaPTk/bin/DCMTK-build" -DBUILD_TESTING=OFF .. && make -j$(nproc) && make install/strip
# install FeTS and remove installer
RUN yes yes | ./FeTS_${VERSION}_Installer.bin --target ./FeTS_${VERSION} -- --cudaVersion 11 && rm -rf ./FeTS_${VERSION}_Installer.bin

ENV PATH=./FeTS_${VERSION}/squashfs-root/usr/bin/:$PATH
ENV LD_LIBRARY_PATH=./FeTS_${VERSION}/squashfs-root/usr/lib/:$LD_LIBRARY_PATH

# set up environment and install correct version of pytorch
RUN cd ./FeTS_${VERSION}/squashfs-root/usr/bin/OpenFederatedLearning && \
    rm -rf ./venv && python3.7 -m venv ./venv && \
    ./venv/bin/pip install torch==1.7.1+cu110 torchvision==0.8.2+cu110 torchaudio==0.7.2 -f https://download.pytorch.org/whl/torch_stable.html 

RUN cd ./FeTS_${VERSION}/squashfs-root/usr/bin/OpenFederatedLearning && \
    ./venv/bin/pip install wheel && \
    ./venv/bin/pip install scikit-build && \
    ./venv/bin/pip install SimpleITK==1.2.4 && \
    ./venv/bin/pip install protobuf==3.17.3 && \
    ./venv/bin/pip install opencv-python==4.2.0.34 && \
    ./venv/bin/pip install python-gdcm

RUN cd ./FeTS_${VERSION}/squashfs-root/usr/bin/OpenFederatedLearning && \
    ./venv/bin/pip install setuptools --upgrade && \
    make install_openfl && \
    make install_fets && \
    ./venv/bin/pip install -e ./submodules/fets_ai/Algorithms/GANDLF && \
    cd ../LabelFusion && \
    rm -rf venv && python3.7 -m venv ./venv && \
    ./venv/bin/pip install -e .

# set up the docker for GUI
ENV QT_X11_NO_MITSHM=1
ENV QT_GRAPHICSSYSTEM="native"

# define entry point
ENTRYPOINT ["/Front-End/bin/install/appdir/usr/bin/BraTSPipeline"]
