FROM cbica/captk_centos7:devtoolset-4_superbuild

LABEL authors="FeTS_Admin <admin@fets.ai>"

RUN yum update -y

RUN yum install git

RUN --mount=type=bind,target=.

RUN echo "running ls -l" && ls -l

# We will do git pull on the FeTS-AI/Front-End master, because that is the repo using which the base image is made
# We will not do compiles on the PR because the idea is that the Xenial build will check the build status of
# the PR in any case.
RUN cd Front-End; \ 
    git pull origin master

RUN cd Front-End; \
    mkdir bin; \
    cmake -DITK_DIR=../../CaPTK/bin/ITK-build -DDCMTK_DIR=../../CaPTK/bin/DCMTK-build -DCMAKE_INSTALL_PREFIX="./install/appdir/usr" -DBUILD_TESTING=OFF ..; \
    make && make install/strip; 

# set up the docker for GUI
ENV QT_X11_NO_MITSHM=1
ENV QT_GRAPHICSSYSTEM="native"

# define entry point
ENTRYPOINT ["/work/Front-End/bin/install/appdir/usr/bin/BraTSPipeline"]