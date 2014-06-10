#! /bin/bash

for i in *
do
# echo "Item $((l++)) : $i"

if [ -d $i ]; then
  cd $i
    echo "$PWD"
    astyle --style=allman -s3 -C -S -K  *.cpp 
    astyle --style=allman -s3 -C -S -K  *.h
  cd ..
fi

done
