#!/bin/sh

# Install script for fpm
wget https://github.com/fortran-lang/fpm/releases/download/v0.6.0/fpm-0.6.0-linux-x86_64 -O fpm
chmod +x fpm
mv -p fpm ~/.local/bin
