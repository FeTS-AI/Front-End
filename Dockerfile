FROM cbica/captk_centos7:devtoolset-4_superbuild

LABEL authors="FeTS_Admin <admin@fets.ai>"

RUN yum update -y

RUN yum install git

RUN echo "running ls -l" && ls -l && pwd

WORKDIR /Front-End

COPY . .

RUN echo "running ls -l" && ls -l && pwd

# set up the docker for GUI
ENV QT_X11_NO_MITSHM=1
ENV QT_GRAPHICSSYSTEM="native"

# define entry point
ENTRYPOINT ["/Front-End/bin/install/appdir/usr/bin/BraTSPipeline"]