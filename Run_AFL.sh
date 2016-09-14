#!/bin/sh
#
# This script:
# - restore AFL state from afl-resume-state.tar.gz
# - run afl-fuzz
# - delete afl/crashes/id=* files
# - delete afl/hangs/id=* files
# - copy inputs file causing crashes to afl/crashes/afl/hangs (after renaming them for Windows compatibility)
# - backup AFL state to afl-resume-state.tar.gz

EXE_PATH=build-afl-make/Examples/Parser/Parser

echo "Run American Fuzzy Lop fuzzer"
if [ -d /tmp/afl-work/ ]; then
    rm -rf /tmp/afl-work/*
else
    mkdir /tmp/afl-work/
fi
if [ -f afl-resume-state.tar.gz ]; then
    tar xzf afl-resume-state.tar.gz -C /tmp/afl-work
    # when resuming, pass '-' for test inputs
    AFL_INPUT="-i-"
else
    AFL_INPUT="-i afl/inputs/"
fi
echo "Press Ctrl+C to interrupt AFL. Re-run this script to resume"
sleep 2
# Allow stopping afl-fuzz with Ctrl+C but resume the execution of this script
trap ' ' INT
echo "Running afl-fuzz $AFL_INPUT -o /tmp/afl-work -- $EXE_PATH"
# Use token dictionnary to generate more sensible input
# See /home/afl/dictionaries/README.dictionaries
afl-fuzz $AFL_INPUT -o /tmp/afl-work -x afl/dictionnary/parser_tokens.dict -- $EXE_PATH

echo "Minifying test case size..."
ls /tmp/afl-work/crashes |grep -v 'README' | while read testpath ; do
  tmin_cmd="afl-tmin -i /tmp/afl-work/crashes/$testpath -o /tmp/afl-work/crashes/$testpath -- $EXE_PATH"
  echo $tmin_cmd
  $tmin_cmd
done

echo "Crash/Hang test case inputs can be found in afl/crashes and afl/hangs"
if [ ! -d afl/crashes ]; then
    mkdir afl/crashes
fi
echo > afl/crashes/id=dummy
rm afl/crashes/id=*
if [ ! -d afl/hangs ]; then
    mkdir afl/hangs
fi
echo > afl/hangs/id=dummy
rm afl/hangs/id=*
./afl_safe_copy.py /tmp/afl-work/crashes afl/crashes
./afl_safe_copy.py /tmp/afl-work/hangs afl/hangs
tar czf afl-resume-state.tar.gz -C /tmp/afl-work .
