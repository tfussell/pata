#!/usr/local/bin/python3

import os
import re
import shutil
import subprocess
import sys

trim_galore_bin = '/Users/thomas/Downloads/trim_galore_v0.3.3/trim_galore'
adapter1, adapter2 = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC', 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'
abyss_bin = '/Users/thomas/Downloads/abyss-1.3.7/bin/abyss-pe'
#trans_abyss_dir = '~/Desktop/PyPipeline/Software/trans-ABySS-v.1.4.4'
trans_abyss_dir = '/Users/thomas/Downloads/trans-ABySS-v.1.4.4'
trans_abyss_env = os.path.join(trans_abyss_dir, 'setup.sh')
trans_abyss_fem_bin = os.path.join(trans_abyss_dir,'wrappers/abyss-ta-filter')
trans_abyss_rmdups_bin = os.path.join(trans_abyss_dir, 'wrappers/abyss-rmdups-iterative')

def trim_galore(working_directory, file1, file2, out_file):
	command = trim_galore_bin + ' -a {} -a2 {} --paired {} {}'.format(adapter1, adapter2, file1, file2)
	print(command)
	out = subprocess.check_output(command, stderr=subprocess.STDOUT, shell=True, cwd=working_directory)
	out_file.write(out.decode('utf8'))

	os.remove(os.path.join(working_directory, file1 + '_trimming_report.txt'))
	os.remove(os.path.join(working_directory, file2 + '_trimming_report.txt'))

	if os.path.splitext(file1)[1] == '.gz':
		new_file1 = os.path.splitext(os.path.splitext(file1)[0])[0] + '_val_1.fq.gz'
		new_file2 = os.path.splitext(os.path.splitext(file2)[0])[0] + '_val_2.fq.gz'
	else:
		new_file1 = os.path.splitext(file1)[0] + '_val_1.fq'
		new_file2 = os.path.splitext(file2)[0] + '_val_2.fq'
		
	return new_file1, new_file2

def abyss(working_directory, file1, file2, kmer, root_name, kmer_directory, out_file):
	command = 'source {}; {} -C k{} s=150 n=5 k={} in=\'{} {}\' name={}'.format(trans_abyss_env, abyss_bin, kmer, kmer, file1, file2, root_name)
	subprocess.check_output(command, shell=True, cwd=working_directory)
	out_file.write(out.decode('utf8'))

def trans_abyss_fem(working_directory, file1, file2, kmer, root_name, kmer_directory, out_file):
	command = 'source {}; {} -i {} -k {} -n {} -o {} -l 151'.format(trans_abyss_env, trans_abyss_fem_bin, kmer_directory, kmer, root_name, kmer_directory)
	subprocess.check_output(command, shell=True, cwd=working_directory)
	out_file.write(out.decode('utf8'))

def trans_abyss_rmdups(working_directory, root_name, out_file):
	command = trans_abyss_rmdups_bin + '-i trans_abyss_fem -n {} -o ./trans_abyss_fem'.format(root_name)
	subprocess.check_output(command, shell=True, cwd=working_directory)
	out_file.write(out.decode('utf8'))

def run(args):
	file1, file2 = args[1].split(',')
	working_directory = os.path.dirname(file1)
	file1, file2 = os.path.basename(file1), os.path.basename(file2)
	
	root_name = os.path.basename(re.sub(r'(.+)(_S[0-9]+)(_L[0-9]+)_R1(_[0-9]+)?\.fastq\.gz', r'\1', file1))
	out_file = open(os.path.join(working_directory, root_name + '.log'), 'w')

	file1, file2 = trim_galore(working_directory, file1, file2, out_file)

	ta_directory = os.path.join(working_directory, root_name + '_trans_abyss')
	if not os.path.isdir(ta_directory):
		os.mkdir(ta_directory)

	for kmer in range(51, 202, 10):
		kmer_directory = os.path.join(ta_directory, 'k{}'.format(kmer))
		if not os.path.isdir(kmer_directory):
			os.mkdir(kmer_directory)
		abyss(working_directory, file1, file2, kmer, root_name, kmer_directory, out_file)
		trans_abyss_fem(working_directory, file1, file2, kmer, root_name, kmer_directory, out_file)
	trans_abyss_rmdups(working_directory, root_name, out_file)

	os.move('{}/{}-contigs.fa {}/{}-contigs.fa'.format(ta_directory, root_name, working_directory, root_name))
	#shutil.rmtree(ta_directory)

if __name__ == '__main__':
	if len(sys.argv) < 2:
		print("Usage: pipeline.py file1,file2")
	else:
		run(sys.argv)