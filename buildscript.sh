
mkdir -p ./bin
cd ./bin

echo "Testing locations"
echo "work and contents:"
ls -lt /work
ls -lt /work/CaPTk
ls -lt /work/CaPTk/bin
ls -lt /work/CaPTk/bin/ITK-build
ls -lt /work/CaPTk/bin/DCMTK-build
ls -lt /work/CaPTk/bin/qt/5.12.1/bin

echo "Exporting environment variables"
export PATH=/work/CaPTk/bin/qt/5.12.1/bin:/work/CaPTk/bin/qt/5.12.1/libexec:$PATH
export CMAKE_PREFIX_PATH=/work/CaPTk/bin/ITK-build:/work/CaPTk/bin/DCMTK-build:/work/CaPTk/bin/qt/5.12.1/lib/cmake/Qt5:$CMAKE_PREFIX_PATH

echo "Running CMake"
cmake -DCMAKE_INSTALL_PREFIX="./install/appdir/usr" -DITK_DIR="/work/CaPTk/bin/ITK-build" -DDCMTK_DIR="/work/CaPTk/bin/DCMTK-build" -DBUILD_TESTING=OFF ..

echo "Running make + install"
make -j4 && make install/strip
