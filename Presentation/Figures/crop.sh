#!/bin/bash
folder=.
for i in $(ls $folder | grep .pdf); do
    filename="${i%.*}"
    pdfcrop $folder/$filename.pdf $folder/$filename.pdf
done
