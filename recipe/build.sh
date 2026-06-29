#!/bin/bash

set -e

# Run autoreconf to generate configure script
autoreconf -fi

# fix tools build failed
sed -i 's/\(LDADD = .*libAFEPack\.so\)$/\1 -lopenblas/' tools/Makefile.am

# Configure the project
./configure --prefix="${PREFIX}"

# Build
make -j${CPU_COUNT:-1}

# Install
make install
