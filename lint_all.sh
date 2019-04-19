#!/usr/bin/env sh
echo_green(){
    echo -e "\e[1;32m$1\e[0m"
}

set -e

cur_dir=$(pwd)
echo --- Running lint ---
cd ${cur_dir}/task4/
make lint


echo_green 'Lint pass!'
