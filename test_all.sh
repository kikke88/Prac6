#!/usr/bin/env sh
echo_green(){
    echo -e "\e[1;32m$1\e[0m"
}

set -e

cur_dir=$(pwd)
echo --- Running test ---
cd ${cur_dir}/task4/
make test

echo_green 'All tests pass!'
