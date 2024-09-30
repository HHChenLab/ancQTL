#!/usr/bin/env python3
import numpy as np
import argparse
import pandas
import gzip
import os
from multiprocessing import Pool
import functools
parser = argparse.ArgumentParser(description="Transfer RFMix2 output to ancQTL local ancestry input",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("input", type=str, nargs='+', help="RFMix2 msp.tsv output")
parser.add_argument("-g", "--gtf", help="gtf/gff file for gene information")
parser.add_argument("-c", "--chr", help="chromosome")
parser.add_argument("-o", "--output", default = './', help="output directory")
parser.add_argument('-d', '--dist', default=1000000, type=int, help='distance(bp) for TSS')
parser.add_argument('-t', '--threads', default=1, type=int, help='number of threads to use')
args = vars(parser.parse_args())
print(args)
#quit()
def read_gtf(gtffile):
    if gtffile.split('[.]')[-1] == '.gz':
        return gzip.open(gtffile, 'rb')
    else:
        return open(gtffile, 'r')
geneTSS = {}
with read_gtf(args['gtf']) as gtf:
    for line in gtf:
        if line[0] == '#':
            continue
        line = line.strip().split('\t')
        if args['chr'] is None or line[1] == args['chr']:
            if line[2] == 'gene':
                ensID = dict(map(lambda x: (x.strip().split(' ')[0], x.strip().split(' ')[1].replace('"','')), line[8].strip(';').split(';')))['gene_id']
                if line[0] not in geneTSS:
                    geneTSS[line[0]] = {}
                geneTSS[line[0]][ensID] = min(int(line[3]), int(line[4]))
print('{0} genes from {1}'.format(sum(map( lambda x: len(geneTSS[x]), geneTSS.keys())), args['gtf']))
def ind_anc(iid, pdt):
    if (pdt[iid]==list(pdt[iid])[0]).all():
        return str(list(pdt[iid])[0])
    else:
        return 'NA'
def ind_dose(iid, pdt, genelen, anclist, ancdict):
    ancdos = list( map( lambda x: '{:.3f}'.format(sum(np.array(pdt[iid].map(ancdict[x]))*genelen)/genelen.sum()), anclist))
    return ';'.join(ancdos)
inputfile = args['input']
if not os.path.isdir(args['output']):
    os.system('mkdir {}'.format(args['output']))
cpu = int(args['threads'])
dist = int(args['dist'])
for rfmixin in inputfile:
    if not os.path.isfile(rfmixin):
        print('Cannot find inputfile: {}'.format(rfmixin))
        continue
    print('reading {}'.format(rfmixin))
    outfile = rfmixin.split('/')[-1]
    with open(rfmixin, 'r') as rfmixhead, open(args['output']+'/'+outfile.replace('msp.tsv','ancQTL.tsv'), 'w') as ansout, open(args['output']+'/'+outfile.replace('msp.tsv','ancQTL.dosage.tsv'), 'w') as dosout:
        header = rfmixhead.readline()
        ansout.write('{}\n'.format(header.strip()))
        anc_dict = {}
        anc = list(header.strip().split(': ')[1].split())
        dosout.write('#Dosage: {}\n'.format(';'.join(map(lambda x: x.split('=')[0], anc))))
        for anc_i in list(range(0, len(anc))):
            anc_dict[anc[anc_i].split('=')[0]] = {}
            for i in list(range(0, len(anc))):
                if i == anc_i:
                    anc_dict[anc[anc_i].split('=')[0]][i] =1
                else:
                    anc_dict[anc[anc_i].split('=')[0]][i] =0
        rfmix = pandas.read_table(rfmixin, header=0 , skiprows=1)
        ids = list(rfmix.columns)[6:]
        allanc = list( map( lambda x: x.split('=')[0], anc))
        print(ids[0:10])
        ansout.write('Gene\tTSS\t{}\n'.format('\t'.join(ids)))
        dosout.write('Gene\tTSS\t{}\n'.format('\t'.join(ids)))
        for CHROM in set(rfmix['#chm']):
            if CHROM not in geneTSS:
                continue
            for ensID in geneTSS[CHROM]:
                tss = geneTSS[CHROM][ensID]
                tss_str = tss - dist
                tss_end = tss + dist
                if tss_str < int(list(rfmix.spos)[1]) or tss_end > int(list(rfmix.epos)[-1]):
                    ansout.write('{0}\t{1}\t{2}\n'.format(ensID, tss, '\t'.join(list(map(lambda x: 'NA', ids)))))
                    dosout.write('{0}\t{1}\t{2}\n'.format(ensID, tss, '\t'.join(list(map(lambda x: 'NA', ids)))))
                else:
                    gene = rfmix[(rfmix['epos']>tss_str) & (rfmix['spos']<tss_end)]
                    gene_len = np.array(gene['epos']) - np.array(gene['spos'])
                    gene_len[0] = np.array(gene['epos'])[0] - tss_str
                    gene_len[-1] = tss_end - np.array(gene['spos'])[-1]
                    get_anc = functools.partial( ind_anc, pdt = gene)
                    get_dose = functools.partial( ind_dose, pdt = gene, genelen = gene_len, anclist = allanc, ancdict = anc_dict)
                    #print('\t'.join(list(map( lambda x: str(ind_anc(x, gene)), ids[0:10]))))
                    #print('\t'.join(list(map( lambda x: str(ind_dose(x, gene, gene_len, anc, anc_dict)), ids[0:10]))))
                    p = Pool(cpu)
                    ansout.write('{0}\t{1}\t{2}\n'.format(ensID, tss, '\t'.join(list(p.map( get_anc, ids)))))
                    dosout.write('{0}\t{1}\t{2}\n'.format(ensID, tss, '\t'.join(list(p.map( get_dose, ids)))))
                    p.close()
                    p.join()
