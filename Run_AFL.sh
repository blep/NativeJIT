#!/bin/sh
#
# This script:
# - restore AFL state from afl-resume-state.tar.gz
# - run afl-fuzz
# - delete afl-crashes/id=* files
# - delete afl-hangs/id=* files
# - copy inputs file causing crashes to afl-crashes/afl-hangs (after renaming them for Windows compatibility)
# - backup AFL state to afl-resume-state.tar.gz

echo "Run American Fuzzy Lop fuzzer"
rm -rf /tmp/afl-work/
mkdir /tmp/afl-work/
tar xzf afl-resume-state.tar.gz -C /tmp/afl-work
echo "Press Ctrl+C to interrupt AFL. Re-run this script to resume"
sleep 2
# Allow stopping afl-fuzz with Ctrl+C but resume the execution of this script
trap ' ' INT
afl-fuzz -i afl-inputs/ -o /tmp/afl-work build-make/Examples/Parser/Parser
echo "Crash test case inputs can be found in afl-crashes"
mkdir afl-crashes
rm afl-crashes/id=*
mkdir afl-hangs
rm afl-hangs/id=*
./afl_safe_copy.py /tmp/afl-work/crashes afl-crashes
tar czf afl-resume-state.tar.gz -C /tmp/afl-work .
