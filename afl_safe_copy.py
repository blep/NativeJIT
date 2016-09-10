#!/usr/bin/env python3
import argparse
import glob
import os.path
import shutil


def main():
    parser = argparse.ArgumentParser(description='''\
Copy crash files renaming them for Windows filename compatibility.

American Fuzzy Lop fuzzer crash test filenames contains ':' which is not legal on Windows. This script
copy the file from source_dir to target_dir and replace ':' by '='.
''')
    parser.add_argument('source_dir', metavar='SOURCE_DIR', help='source directory')
    parser.add_argument('target_dir', metavar='TARGET_DIR', help='target directory')
    args = parser.parse_args()
    for p in glob.glob(os.path.join(args.source_dir, '*')):
        safe_filename = os.path.basename(p).replace(':', '=')
        shutil.copyfile( p, os.path.join( args.target_dir, safe_filename ) )

if __name__ == '__main__':
    main()
