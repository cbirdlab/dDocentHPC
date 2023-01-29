#!/bin/bash

awk '
	{
		if (NR == 1) {
			# save first file name in file to variable
			filename = $0
			filename_len = length($0)
			n = filename_len
		} else {
			next_filename = $0
			next_filename_len = length($0)
			if (n > next_filename_len) {
				n = next_filename_len
			} 
			for (i = 1; i <= n; i++) {
				if (substr(filename, filename_len - i + 1, i) != substr(next_filename, next_filename_len - i + 1)) {
					common_suffix = substr(filename, filename_len - i + 2)
					break
				}
				if (i = n) {
					common_suffix = substr(filename, filename_len - i + 1)
				}
			}
		#print i, common_suffix
		}
	}
	END {
		print common_suffix
	}
' 
