#!/usr/local/bin/python3

import os
import re
import subprocess
import sys

def run(args):
	directory = os.getcwd() if len(args) < 2 else args[1]

	matches = [i for i in os.listdir(directory) if re.match('^.*_R1(_001)?\.fastq\.gz', i)]
	possible_pairs = [(os.path.join(directory, i), os.path.join(directory, re.sub(r'(.+)_R1(.+)', r'\1_R2\2', i))) for i in matches]

	pairs = []
	for pair in possible_pairs:
		if not os.path.isfile(pair[1]):
			print("No paired reads found for:", pair[0], "skipping...")
		else:
			pairs.append(pair)

	command = ['parallel', '--xapply', '-j2', './pipeline.py', ':::']
	command += [pair[0] for pair in pairs]
	command.append(':::')
	command += [pair[1] for pair in pairs]

	subprocess.call(command)

if __name__ == '__main__':
	run(sys.argv)