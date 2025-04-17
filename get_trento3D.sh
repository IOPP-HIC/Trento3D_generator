#!/usr/bin/env bash


rm -fr trento3d_subnucleon
git clone --depth=1 https://github.com/xyw2016/trento3d_KW.git -b ncoll trento3d_subnucleon
commitHash="d0e4bed9138616432db11eb52caa900e269b453c"
cd trento3d_subnucleon 
git checkout $commitHash
rm -fr .git
cd ..
