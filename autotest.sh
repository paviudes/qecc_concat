# Script to make and run an executable.
# Sometimes the script may return an error message
# In those cases, run the following command on the shell: sed -i -e 's/\r$//' automake.sh
#
#!/bin/bash
make clean >/dev/null 2>&1
make db=0 > make.log 2>&1
make_errors=$(cat make.log | grep -i "error")
make_warns=$(cat make.log | grep -i "warning")
if [[ ! -z ${make_errors} ]]; then
	echo "\033[4mBuild failed after producing the following output.\033[0m"
	echo "xxxxxxxx"
	cat make.log
	echo "xxxxxxxx"
else
	if [[ ! -z ${make_warns} ]]; then
		echo "\033[4mBuild passed with warnings, after producing the following output.\033[0m"
		echo "xxxxxxxx"
		cat make.log
		echo "xxxxxxxx"
		echo "Running bmark anyways."
	else
		echo "Make succesful."
	fi
	echo "xxxxxxxx"
	./bmark benchmark Benchmark
fi