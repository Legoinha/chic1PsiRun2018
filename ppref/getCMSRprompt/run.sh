#!/bin/bash

g++ getCMSRprompt.cc $(root-config --libs --cflags) -I"../../includes/" -g -o getCMSRprompt.exe || exit 1

./getCMSRprompt.exe

rm getCMSRprompt.exe
