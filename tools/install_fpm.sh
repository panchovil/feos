#!/bin/sh

install() {
    wget https://github.com/fortran-lang/fpm/releases/download/v0.6.0/fpm-0.6.0-linux-x86_64 -O fpm
    chmod +x fpm
}

run_test() {
    ./fpm test
}

case $1 in
    "install")  install;;
    "test") run_test;;
esac
