
import argparse
import pandas as pd
from Bio import SeqIO


#This file takes in a vcf file and a template genome
#It gives the template all the reference alleles to create a reference genome

def ref_to_fasta(vcf_file, template_file, out_file):
    #Read in the vcf file
    vcf = pd.read_csv(vcf_file, sep = '\t', comment = '#', header = None)
    vcf.rename(columns = {0: 'CHROM', 1: 'POS', 3: 'REF', 4: 'ALT'}, inplace = True)

    #Read in the template genome
    template = SeqIO.read(template_file, 'fasta')
    template_genome = list(template.seq)
    


    for index, row in vcf.iterrows():
        #Get the position and the reference allele
        pos = row['POS']
        ref = row['REF']
        #Change the reference allele in the template genome
        template_genome[pos-1] = ref

    #Write the new genome to a fasta file
    with open(out_file, 'w') as f:
        f.write('>' + template.id + '\n')
        f.write(''.join(template_genome) + '\n')



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert a vcf file to a fasta file')
    parser.add_argument('vcf', type=str, help='The vcf file with the reference alleles')
    parser.add_argument('template', type=str, help='The template genome to put mutations on')
    parser.add_argument('--output', type=str, help='The output fasta file', default='reference_genome.fasta')
    args = parser.parse_args()

    ref_to_fasta(args.vcf, args.template, args.output)

