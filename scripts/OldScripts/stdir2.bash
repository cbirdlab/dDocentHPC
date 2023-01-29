#!/bin/bash

#this script makes the birdlab directory structure for processing ddRAD and ezRAD data

#to run, copy this script into your project directory and type the following at the command prompt
#bash stdir.sh

mkdir assembly
mkdir assembly/FASTQC
mkdir assembly/logfiles
mkdir assembly/trimreports
mkdir fastq
mkdir FASTQC
mkdir highReadNum
mkdir lowReadNum
mkdir mapping
mkdir mapping/FASTQC
mkdir mapping/logfiles
mkdir mapping/trimreports
mkdir removed_seqs
mkdir filtering

