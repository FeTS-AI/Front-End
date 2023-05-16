#!/bin/bash

ls -lt; ls -lt ~/work; docker run --rm --entrypoint="mkdir bin; cmake -DITK_DIR=~/work/CaPTK/bin/ITK-build -DDCMTK_DIR=~/work/CaPTK/bin/DCMTK-build -DCMAKE_INSTALL_PREFIX="./install/appdir/usr" -DBUILD_TESTING=OFF ..; make && make install/strip;" -v "$(pwd)"/Front-End:/Front-End ghcr.io/fets-ai/front-end:latest