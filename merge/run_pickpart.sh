#!/bin/bash

g++ pickpart.C $(root-config --libs --cflags) -g -o pickpart.exe || exit 1
./pickpart.exe $1 $2 $3
rm pickpart.exe