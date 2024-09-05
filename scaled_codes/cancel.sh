#!/usr/bin/env bash

for jobid in {316761..318365}
do
    llcancel $jobid
done