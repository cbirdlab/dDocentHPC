#!/bin/bash

#count reads of any file ending in fq.gz in the present directory
ls *fq.gz | parallel --no-notice "echo -n {}' ' && zgrep -c '^\+$' {}" | tr " " "\t" | sort -g -k 2 > NumFqGzReads.dat

#visualize read counts
cut -f2 NumFqGzReads.dat > NumFqGzReads

gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale
unset label
set logscale x
set title "Histogram of Reads Per Individual Per Direction"
set ylabel "Number of Individuals x Direction"
set xlabel "Depth of Coverage"
xmax="`tail -1 NumFqGzReads`"
xmin="`head -1 NumFqGzReads`"
binwidth=(xmax-xmin)/30
bin(x,width)=width*floor(x/width) + binwidth/2.0
plot 'NumFqGzReads' using (bin($1,binwidth)):(1.0) smooth freq with boxes
pause -1
EOF

gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale 
unset label
set title "Scatter plot of coverage per individual x direction"
set ylabel "NumReads"
set xlabel "Individual x Direction"
set logscale y
plot 'NumFqGzReads' pt "*" 
pause -1
EOF

rm NumFqGzReads