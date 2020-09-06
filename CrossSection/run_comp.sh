#!/bin/bash

g++ comp_results.C $(root-config --libs --cflags) -I$HOME -I"../includes/" -g -o comp_results.exe || exit 1

./comp_results.exe 

rm comp_results.exe
