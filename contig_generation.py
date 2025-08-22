import sys
import string
import os
import subprocess
import re
import textwrap
import time


file1 = sys.argv[1]
file2 = sys.argv[2]
file3 = sys.argv[3]
file4 = sys.argv[4]
file5 = sys.argv[5]
res = ""
count = 1
fai_list = {}

start_time = time.time()

def format_header(header, pos):
    match = re.match(r'>(\w+):(\d+)-(\d+)', header)
    if match:
        chrom = (match.groups())[0]
        return f'>{chrom}_contig_{pos}'
    return header


with open(file1, 'r') as fasta, open(file2, 'r') as vcf, open(file3, 'w') as new, open(file4, 'w') as chainfile, open(file5, 'r') as fai:
        
        while string != "":
            string = fasta.readline()

            if string == "":
                break

            new.write(string)


        string = "start"

        for line in fai:
            columns = line.strip().split()
            if len(columns) >= 2:
                key = columns[0]
                value = columns[1]
                fai_list[key] = value


        while string != "":
            string = vcf.readline()

            if string == "":
                break

            tmp = string
            full_sequence = []

            if tmp[0] != '#':
                
                tmp = tmp.split()
                chr = tmp[0]
                pos = int(tmp[1])
                pos_start = pos - 455
                pos_start = pos_start if pos_start >= 0 else 0
                pos_end = pos + 524

                cmd = f"samtools faidx {file1} {chr}:{pos_start}-{pos - 1}"
                ps = subprocess.run(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT, text=True)
                for line in ps.stdout.splitlines():
                    if line.startswith('>'):
                        line = format_header(line, pos)
                        contig_header = line
                        new.write(f'{line}\n')
                    else:
                        full_sequence.append(line.strip())

                full_sequence.append(tmp[4].strip())
                
                cmd = f"samtools faidx {file1} {chr}:{pos + 1}-{pos_end}"
                ps = subprocess.run(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT, text=True)
                for line in ps.stdout.splitlines():
                    if line.startswith('>'):
                        continue
                    else:
                        full_sequence.append(line.strip())
        
                combined_sequence = ''.join(full_sequence)

                wrapped_sequence = textwrap.fill(combined_sequence, width=70)
                new.write(wrapped_sequence)
                new.write('\n')

                chainfile.write(f'chain {len(combined_sequence)} {contig_header[1:]} {len(combined_sequence)} + 0 {len(combined_sequence) - 1} {contig_header[1:5]} {fai_list.get(contig_header[1:5])} + {pos_start - 1} {pos_end} {count}\n')
                count += 1
                chainfile.write(f'456 {len(tmp[4]) - len(tmp[3])}\t0\n{len(combined_sequence) - 1 - 456 - (len(tmp[4]) - len(tmp[3]))}\n\n')
        
        for chrom, chrom_len in fai_list.items():
            chainfile.write(f'chain {chrom_len} {chrom} {chrom_len} + 0 {chrom_len} {chrom} {chrom_len} + 0 {chrom_len} {count}\n')
            count += 1
            chainfile.write(f'{chrom_len}\n\n')


end_time = time.time() 
elapsed_time = end_time - start_time 
print(f"Программа выполнялась {elapsed_time:.2f} секунд")