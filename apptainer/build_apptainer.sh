#!/bin/bash
apptainer build --sandbox --fix-perms --fakeroot test.sandbox dbgc.def
apptainer build dbgc.sif test.sandbox/
rm -rf test.sandbox/
