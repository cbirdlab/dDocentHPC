#!/bin/bash

rm namelist*
cp ../configs/config.6.rad ./config.6.rad_test
sed -i -e 's/^2\(\t\+Cutoff1\)/test\1/' \
	-e 's/^2\(\t\+Cutoff2\)/basic\1/' \
	-e 's/^PE\t/PE\t/' \
	-e 's/^75\(\t\+trimmomatic MINLEN\)/33\1/' \
	-e 's/^40\(\t\+Number of Processors\)/10\1/' \
	-e 's/^50\(\t\+bwa mem \-T\)/0\1/' config.6.rad_test
wait
bash ../dDocentHPC_dev.bash mkBAM config.6.rad_test
