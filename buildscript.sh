
mkdir -p ./bin
cd ./bin

echo "Exporting environment variables"
export PATH=/opt/qt/5.11.2/gcc_64/bin:/opt/qt/5.11.2/gcc_64/libexec:$PATH
export CMAKE_PREFIX_PATH=/work/CaPTk/bin/ITK-build:/work/CaPTk/bin/DCMTK-build:/opt/qt/5.11.2/gcc_64/lib/cmake/Qt5:$CMAKE_PREFIX_PATH

echo "Running CMake"
cmake -DCMAKE_INSTALL_PREFIX="./install/appdir/usr" -DITK_DIR="/work/CaPTk/bin/ITK-build" -DDCMTK_DIR="/work/CaPTk/bin/DCMTK-build" -DBUILD_TESTING=OFF ..

echo "Running make + install"
make -j4 && make install/strip
