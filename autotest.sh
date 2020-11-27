#!/bin/bash
make clean ip=1
make ip=1 > make_log.txt 2>make_errors.txt
if [[ -s make_errors.txt ]]; then
	echo "Build failed"
	cat make_log.txt
	cat make_errors.txt
else
	./bmark benchmark Benchmark
fi
