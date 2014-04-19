#!/usr/local/bin/python3

import os
import re
import shutil
import subprocess
import sys

script_path = os.path.dirname(os.path.realpath(__file__))
trim_galore_bin = os.path.join(script_path, 'software', 'trim_galore_v0.3.1', 'trim_galore')
adapter1, adapter2 = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC', 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'
abyss_bin = 'abyss-pe'
trans_abyss_dir = os.path.join(script_path, 'software', 'trans-ABySS-v1.4.4')
trans_abyss_env = os.path.join(trans_abyss_dir, 'setup.sh')
trans_abyss_fem_bin = os.path.join(trans_abyss_dir,'wrappers/abyss-ta-filter')
trans_abyss_rmdups_bin = os.path.join(trans_abyss_dir, 'wrappers/abyss-rmdups-iterative')

def trim_galore(working_directory, file1, file2, out_file):
    command = '/usr/bin/perl ' + trim_galore_bin + ' -a {} -a2 {} --paired {} {}'.format(adapter1, adapter2, file1, file2)
    print('\'' + command + '\'')
    out = subprocess.check_output(command, stderr=subprocess.STDOUT, shell=True, cwd=working_directory)
    out_file.write(out.decode('utf8', 'ignore'))

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
    command = 'source {}; {} -C {} s=150 n=5 k={} name={} in=\'../../{} ../../{}\''.format(trans_abyss_env, abyss_bin, kmer_directory, kmer, root_name, file1, file2)
    print(command)
    try:
        out = subprocess.check_output(command, stderr=subprocess.STDOUT, shell=True, cwd=working_directory)
        out_file.write(out.decode('utf8', 'ignore'))
    except subprocess.CalledProcessError as e:
        pass

def trans_abyss_fem(working_directory, file1, file2, kmer, root_name, kmer_directory, out_file):
    command = 'source {}; {} -i {} -k {} -n {} -o {} -l 151'.format(trans_abyss_env, trans_abyss_fem_bin, kmer_directory, kmer, root_name, kmer_directory)
    print(command)
    out = subprocess.check_output(command, stderr=subprocess.STDOUT, shell=True, cwd=working_directory)
    out_file.write(out.decode('utf8', 'ignore'))

def trans_abyss_rmdups(working_directory, ta_directory, root_name, out_file):
    command = trans_abyss_rmdups_bin + ' -i {} -n {} -o {}'.format(ta_directory, root_name, ta_directory)
    print(command)
    out = subprocess.check_output(command, stderr=subprocess.STDOUT, shell=True, cwd=working_directory)
    out_file.write(out.decode('utf8', 'ignore'))

def get_read_length(filename):
    if os.path.splitext(filename)[1] == '.gz':
        import gzip
        file = gzip.open(filename)
    else:
        file = open(filename)
    line = file.readline()
    line = file.readline()
    return len(line.strip())

def run(file1, file2):
    working_directory = os.path.dirname(file1)
    file1, file2 = os.path.basename(file1), os.path.basename(file2)

    read_length = get_read_length(os.path.join(working_directory, file1))    
    root_name = os.path.basename(re.sub(r'(.+)(_S[0-9]+)(_L[0-9]+)_R1(_[0-9]+)?\.fastq(\.gz)?', r'\1', file1))
    out_file = open(os.path.join(working_directory, root_name + '.log'), 'w')

    print('Assembling {}'.format(root_name))

    try:
        file1, file2 = trim_galore(working_directory, file1, file2, out_file)
    except subprocess.CalledProcessError as e:
        print('Error in trim_galore: ' + str(e.output))
        return

    ta_directory = root_name + '_trans_abyss'
    if not os.path.isdir(os.path.join(working_directory, ta_directory)):
        os.mkdir(os.path.join(working_directory, ta_directory))

    upper_bounds = {150 : 131, 250 : 201, 300 : 251}
    try:
        upper_bound = upper_bounds[read_length]
    except:
        out_file.write('Error: unknown read length, {}, stopping assembly\n'.format(read_length))
        print('Error during assembly of {}: unknown read length, {}, stopping assembly'.format(root_name, read_length))
        return

    out_file.write('Detected read length: {}\n'.format(read_length))
    out_file.write('Assembling with k-mer in range 51 to {} (multiples of 10)\n'.format(upper_bound))

    for kmer in range(51, upper_bound + 1, 10):
        kmer_directory = os.path.join(ta_directory, 'k{}'.format(kmer))
        if not os.path.isdir(os.path.join(working_directory, kmer_directory)):
            os.mkdir(os.path.join(working_directory, kmer_directory))
        abyss(working_directory, file1, file2, kmer, root_name, kmer_directory, out_file)
        if not os.path.isfile(os.path.join(working_directory, kmer_directory, root_name + '-contigs.fa')):
            out_file.write('No contigs file created for k-mer {}, not attempting any higher k-mer\n'.format(kmer))
            shutil.rmtree(os.path.join(working_directory, kmer_directory))
            break
        trans_abyss_fem(working_directory, file1, file2, kmer, root_name, kmer_directory, out_file)
    trans_abyss_rmdups(working_directory, ta_directory, root_name, out_file)

    os.rename('{}/{}-contigs.fa'.format(os.path.join(working_directory, ta_directory), root_name), '{}/{}-contigs.fa'.format(working_directory, root_name))
    shutil.rmtree(os.path.join(working_directory, ta_directory))
    os.remove(os.path.join(working_directory, file1))
    os.remove(os.path.join(working_directory, file2))
    print('Finished {}'.format(root_name))

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: assemble.py file1 file2")
    else:
        run(sys.argv[1], sys.argv[2])
