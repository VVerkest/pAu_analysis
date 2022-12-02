#!/bin/bash
for ppdAu in pp dAu; 
do
    for which in 2 3 4;
    do
       echo ./calc_NchHT_perEA_$ppdAu.cc $which
       ./calc_NchHT_perEA_$ppdAu.cc $which
    done;
done;
