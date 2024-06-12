#!/bin/bash

# go to io directory
cd /io

# install pip
yum -y update
# yum install python3-pip python3-devel -y

PYTHON_VERSIONS=("cp36-cp36m" "cp37-cp37m" "cp38-cp38" "cp39-cp39" "cp310-cp310" "cp311-cp311" "cp312-cp312")

for PYTHON_VERSION in "${PYTHON_VERSIONS[@]}"; do
    export PATH="/opt/python/${PYTHON_VERSION}/bin:$PATH"
    
    python3 --version
    pip3 --version

    # upgrade pip and install build and wheel package
    python3 -m pip install --upgrade pip
    python3 -m pip install build

    ls -alhr

    # Run the build command
    python3 -m build --wheel

done

# manylinux version
auditwheel show ./dist/*
auditwheel repair ./dist/*
