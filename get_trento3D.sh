#!/usr/bin/env bash


rm -fr trento3d_subnucleon
git clone --depth=1 https://github.com/Duke-QCD/trento3d-1.0.git  -b master trento3d_subnucleon
commitHash="7affda6752318f6fdf148ba70413936453bb85cf"
cd trento3d_subnucleon 
git checkout $commitHash
rm -fr .git
cd ..
