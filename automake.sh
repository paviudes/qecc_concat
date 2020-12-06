# Script to make and copy the shared object file to chflow.
# Sometimes the script may return an error message
#
# In those cases, run the following command on the shell: sed -i -e 's/\r$//' automake.sh

#!/bin/bash
make clean
make > make_log.txt 2>make_errors.txt
if [[ -s make_errors.txt ]]; then
	echo "Build failed"
	cat make_log.txt
	cat make_errors.txt
else	
	cat make_log.txt
	echo "Copying bmark.so to chflow/src/simulate"
	cp bmark.so ./../chflow/src/simulate/
	ls -lt ./../chflow/src/simulate/bmark.so
fi
