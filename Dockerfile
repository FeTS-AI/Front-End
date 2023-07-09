FROM ghcr.io/fets-ai/fetstool_docker_dependencies

LABEL authors="FeTS_Admin <admin@fets.ai>"

RUN apt-get update && apt-get update --fix-missing

# RUN apt-get install git

# RUN echo "running ls -l" && ls -l \
#     ls -l /CaPTk \
#     ls -l /CaPTk/bin \
#     ls -l /bin

ENV PATH=/CaPTk/bin/qt/5.12.1/bin:/CaPTk/bin/qt/5.12.1/libexec:$PATH
ENV CMAKE_PREFIX_PATH=/CaPTk/bin/ITK-build:/CaPTk/bin/DCMTK-build:/CaPTk/bin/qt/5.12.1/lib/cmake/Qt5:$CMAKE_PREFIX_PATH

WORKDIR /Front-End

COPY . .

RUN echo "running ls -l" && ls -l && pwd

RUN mkdir bin && cd bin \
    cmake -DCMAKE_INSTALL_PREFIX="./install/appdir/usr" -DITK_DIR="/CaPTk/bin/ITK-build" -DDCMTK_DIR="/CaPTk/bin/DCMTK-build" -DBUILD_TESTING=OFF .. \
    make -j$(nproc) \
    make install/strip

# RUN bash ./buildscript.sh

# set up the docker for GUI
ENV QT_X11_NO_MITSHM=1
ENV QT_GRAPHICSSYSTEM="native"

# define entry point
ENTRYPOINT ["/Front-End/bin/install/appdir/usr/bin/BraTSPipeline"]
