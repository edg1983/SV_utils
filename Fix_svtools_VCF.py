'''
Fix formatting issues and float GQ values produced by svtool pipeline
- Float GQ values are rounded to int
- Dot values in sample fields are converted to zero
- Poor variants lacking proper FORMAT fields are removed

Author: Edoardo Giacopuzzi
'''

from collections import OrderedDict
import argparse
import gzip
import sys, os

def checkFile(file):
    if not os.path.isfile(file):
        sys.exit("CRITICAL! " + file + " do not exists!")
    
def parseFORMAT(format_string, sample_string):
    vcfFORMAT = format_string.split(":")
    vcfSAMPLE = sample_string.split(":")
    values = OrderedDict(zip(vcfFORMAT, vcfSAMPLE))
    return values

def updateDotFORMAT(sample_values, new_string):
    updated = 0
    keys = [x for x in sample_values.keys() if x != "GT"]
    for k in keys:
        if sample_values[k] == '.':
            updated = 1
            sample_values[k] = new_string
    return updated, sample_values

def roundFORMATfield(sample_values, field):
    updated = 0
    if field in sample_values.keys() and sample_values[field] != '.':
        updated = 1
        sample_values[field] = round(float(sample_values[field]))
    return updated, sample_values

def tokenize(line,sep):
    line = line.rstrip('\n')
    line = line.split(sep)
    return line

parser = argparse.ArgumentParser(
    description='Round GQ values to make them integers, so they can be processed dy bcftools')
parser.add_argument('-v','--vcf', action='store', required=True,
                   help='VCF file with structural vars (.vcf / .vcf.gz)')
parser.add_argument('-o','--out', action='store', required=True,
                   help='output file name')                  
args = parser.parse_args()

checkFile(args.vcf)
if args.vcf.endswith("vcf.gz"):
    vcf = gzip.open(args.vcf,"rt")
elif args.vcf.endswith("vcf"):
    vcf = open(args.vcf,"r")

outfile = open(args.out, "w+")

line = vcf.readline()
while line.startswith('#'):
    outfile.write(line)
    line = vcf.readline()

nvars = 0
ngenos = 0
gq_rounded = 0
dot_changed = 0
output_vars = 0
while line:
    nvars += 1
    line = tokenize(line, "\t")
    format_string = line[8]
    #When FORMAT contains only GT or GT:CN, it likely represent a poor call from Lumpy
    #These calls can not be processed further so they are discarded
    if format_string != "GT" and format_string != "GT:CN":
        output_vars += 1
        for i in range(9, len(line)):
            ngenos += 1
            sample_values = parseFORMAT(format_string,line[i])

            #svtools generate float GQ values, but most tools requires int
            #so we round GQ values when they are float 
            is_updated, new_values = roundFORMATfield(sample_values,'GQ')
            gq_rounded += is_updated

            #Calls with no genotype or very low support result in . values in GQ and other fields
            #This breaks processing with downstream tools, so . values are converted to zero
            is_updated, new_values = updateDotFORMAT(sample_values,"0")
            dot_changed += is_updated

            line[i] = ":".join([str(x) for x in new_values.values()])            
        outfile.write("\t".join(line) + "\n")
    line = vcf.readline()

outfile.close()
vcf.close()   
print(nvars, " variants in input VCF")
print(output_vars, " variants with proper format")
print(ngenos, " total genotypes in output vars")
print(gq_rounded, " GQ values rounded")
print(dot_changed, " dot values converted to zero")