#!/usr/bin/env python3

import sys,os,subprocess,shutil,statistics
from Bio import SeqIO
from collections import defaultdict

clusters = sys.argv[1]
reffile = sys.argv[2]


#split_ref
tmp_ref_dir = '.tmprefdir'
if os.path.isdir(tmp_ref_dir):
    shutil.rmtree(tmp_ref_dir)
os.makedirs(tmp_ref_dir)
ref = SeqIO.parse(reffile, "fasta")
fout_reflist = []
for rec in ref:
    SeqIO.write([rec], os.path.join(tmp_ref_dir,"{}.fa".format(rec.id)), 'fasta')
    fout_reflist.append(rec.id)
refseqpaths = [os.path.join(tmp_ref_dir,f) for f in os.listdir(tmp_ref_dir)]
reffilelist = os.path.join(tmp_ref_dir,"reflist.txt")
with open(reffilelist, 'w') as rfl:
    rfl.write('\n'.join(refseqpaths))

#run fastANI
total_output = {}
for cluster in os.listdir(clusters):
    if not cluster.startswith('cluster'):
        continue
    clusterd_path = os.path.join(clusters, cluster)
    cluster_seq_filepaths = [os.path.join(clusterd_path, '1_contigs', f) for f in os.listdir(os.path.join(clusterd_path, '1_contigs'))]
    quefilelist = os.path.join(tmp_ref_dir,"quelist.txt")
    with open(quefilelist, 'w') as qfl:
        qfl.write('\n'.join(cluster_seq_filepaths))
    outputfile = os.path.join(tmp_ref_dir,"output.txt")
    subprocess.run(['fastANI', '--ql', quefilelist, '--rl', reffilelist, '-o', outputfile, '-t', '16'])
    output_d = defaultdict(list)
    with open(outputfile) as oh:
        for line in oh:
            data = line.strip().split()
            output_d[os.path.basename(data[1])[:-3]].append(float(data[2]))
    total_output[cluster] = output_d

with open("final_output.tsv", 'w') as fout:
    fout.write("Cluster\t{}\n".format('\t'.join(fout_reflist)))
    for cluster in sorted(total_output.keys()):
        outline = []
        outline.append(cluster)
        for seq in fout_reflist:
            if seq in list(total_output[cluster]):
                outline.append('{:.2f}'.format(statistics.mean(total_output[cluster][seq])))
            else:
                outline.append('na')
        fout.write('\t'.join(outline)+'\n')
