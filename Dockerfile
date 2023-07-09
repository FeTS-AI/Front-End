FROM ghcr.io/fets-ai/fetstool_docker_dependencies

LABEL authors="FeTS_Admin <admin@fets.ai>"

RUN apt-get update && apt-get update --fix-missing

RUN apt-get install git

RUN echo "running ls -l" && ls -l && pwd

WORKDIR /Front-End

COPY . .

RUN echo "running ls -l" && ls -l && pwd

RUN bash ./buildscript.sh

# set up the docker for GUI
ENV QT_X11_NO_MITSHM=1
ENV QT_GRAPHICSSYSTEM="native"

# define entry point
ENTRYPOINT ["/Front-End/bin/install/appdir/usr/bin/BraTSPipeline"]
