#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
SYNCDIR="${DIR}"

APPDIR="Library/Mathematica/Applications"
PKG="MDS"
TO=${HOME}/${APPDIR}/${PKG}

dirs="Kernel"
files="MDS.m
LICENSE
README.md"

echo "installing to ${TO}"

if [[ ! -d ${TO} ]] ; then
    mkdir ${TO}
fi

for f in $dirs
do
  echo "Processing $f..."

  ln -sf ${SYNCDIR}/$f ${TO}/$f 
done

for f in $files
do
  echo "Processing $f..."

  ln -sf ${SYNCDIR}/$f ${TO}/$f 
done
