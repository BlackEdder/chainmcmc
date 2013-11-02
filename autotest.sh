#!/bin/bash
# Immediately run once
ctest

trap ctest INT

function watch_tests() {
while : 
do
	#file=`inotifywait -q -e CREATE bin/ --format %f`
	file=`inotifywait -q -e ATTRIB bin/ --format %f`
	if [[ "$file" =~ test_* ]]
	then
		sleep 2
		ctest
	fi
done
}

watch_tests
