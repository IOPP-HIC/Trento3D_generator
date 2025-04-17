#!/usr/bin/env bash


rm -fr trento3d_subnucleon
git clone --depth=1 git@github.com:xyw2016/trento3D_subnucleon.git -b master trento3d_subnucleon
commitHash="dded85bbe1a06eb93abcfbe406cf598226d22894"
cd trento3d_subnucleon 
git checkout $commitHash
rm -fr .git
cd ..
