#!/usr/bin/env bash


rm -fr trento3d_subnucleon
git clone --depth=1 https://github.com/xyw2016/trento3D_subnucleon.git -b master trento3d_subnucleon
commitHash="d6913d93053fef8729e988b10a1a706d6e17fa8c"
cd trento3d_subnucleon 
git checkout $commitHash
rm -fr .git
cd ..
