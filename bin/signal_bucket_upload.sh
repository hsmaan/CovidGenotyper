#!/bin/bash

gsutil rm -r gs://signal-umap/*

cd ../data

mkdir latest

cp signal_aligned* latest
cp signal_dist* latest
cp signal_meta* latest 
cp signal_umap* latest

gsutil cp -r latest gs://signal-umap/

rm -r latest 

