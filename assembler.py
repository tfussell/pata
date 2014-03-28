#!/usr/local/bin/python3

def run(args):
	import os
	directory = os.getcwd() if len(args) < 2 else args[1]

	import re
	matches = [i for i in os.listdir(directory) if re.match('^.*_R1(_001)?\.fastq\.gz', i)]
	possible_pairs = [(os.path.join(directory, i), os.path.join(directory, re.sub(r'(.+)_R1(.+)', r'\1_R2\2', i))) for i in matches]

	pairs = []
	for pair in possible_pairs:
		if not os.path.isfile(pair[1]):
			print("No paired reads found for:", pair[0], "skipping...")
		else:
			pairs.append(pair)

	import subprocess
	subprocess.call(['parallel', '--gnu', '-j2', './pipeline.py', ':::'] + ['{},{}'.format(*pair) for pair in pairs])

if __name__ == '__main__':
	import sys
	run(sys.argv)