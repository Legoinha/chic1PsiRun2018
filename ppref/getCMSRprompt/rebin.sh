#!/bin/bash

g++ rebin.cc $(root-config --libs --cflags) -I"../../includes/" -g -o rebin.exe || exit 1

./rebin.exe

rm rebin.exe 
