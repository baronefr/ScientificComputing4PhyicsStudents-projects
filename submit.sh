#!/bin/bash -x

for _dir in "$1-"*; do
    # this code compresses only the first match of $1-*
    tar -cvf "$1-handout.tar" "$_dir"
    break
done


