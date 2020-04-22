#####################################################################################################################
## Code for annotating variants and parsing allele associated calculations on per variant and per sample basis.     #
## Written by P.Mankoo                                                                                              #
#####################################################################################################################

import sys
import requests
import json
import os 
import subprocess
from ToolConfig import snpeff_jar, snpeff_config_file
import logging
import time
start_time = time.time()

#################################################
print("Start Variant Annotations using SnpEff\n")
#################################################

def run_snpeff(target_vcf):
    
    outfile_name = target_vcf[:-4]
    
    snpeff_input = target_vcf
    snpeff_output = outfile_name + ".snpeff.vcf"
  
    annotate_cmd = """java -Xmx2g -jar %s hg19 -hgvs -c %s -v %s > %s""" % (snpeff_jar, snpeff_config_file, snpeff_input, snpeff_output)
    
    logging.debug( annotate_cmd ) 
    result = subprocess.call(annotate_cmd, shell=True, stdout = open("log_out.txt", "w"), stderr = open("log_err.txt", "w"))
    
    if result > 0:
        logging.critical( annotate_cmd )
        logging.critical( "SnpEff failed. Is everything OK with your SnpEff installation??" )
        sys.exit()
    
    return snpeff_output

# Uncomment the following line to run snpEff yourself!
#run_snpeff('Challenge_data.vcf')

######################################################################################
print("SnpEff Annotations written in file Challenge_data.snpeff.vcf.\n")
######################################################################################

#########################################################################################################################################################
print("Let's extract allele frequencies from Exac REST API")
print("Running each variant, with added complexity of parsing through multiallelic variants through EXAC is time consuming.")
print("Therefore, this step is converted to BULK API version and resulting allele frequencies stored in a list to be printed later.\n")
#########################################################################################################################################################

inputvariants_for_exac=[]
VCF = open('Challenge_data.vcf', "r")
for line in VCF:
    if line.startswith("#CHROM"):
        # First parse column headers into keys to store each column
        fields = line.strip("# \n").split("\t")
    elif not line.startswith("#"):
        # Second extract data associated with header columns
        variables = line.strip().split("\t")
        toks = dict(zip(fields, variables))

        # Remember to loop over multi allelic variants and add to the list of variants to be passed to EXAC!!
        multialleles = toks["ALT"].split(",")
        for i in range(0,len(multialleles)):
            inputvariants_for_exac.append(r'\"' + "-".join((toks["CHROM"], toks["POS"], toks["REF"], multialleles[i])) + r'\"')


lengthvar=len(inputvariants_for_exac)

# The length of this list is longer than the variants in the original vcf file referring to many multiallelic variants.
# Later we will also print out the number of alleles per variant to make sure we can QC the output appropriately.
# It is good to test the input(s) to be sent to API before running the API

# print(inputvariants_for_exac[0],inputvariants_for_exac[6000])

exac_allele_freq = {}
url = 'http://exac.hms.harvard.edu/rest/bulk/variant'
for i in range(0, lengthvar, 2000):
    data = r'"[' + ",".join(inputvariants_for_exac[i:i+2000]) + r']"'
    # It is a good idea to check the progress while running APIs.
    if i % 2000 == 0: print("At Step",i)
    r = requests.post(url, eval(data))
    exac_annotations = json.loads(r.text)
    for line in exac_annotations:
        position = "-".join(line.split('-')[0:4])
        if 'allele_freq' in exac_annotations[line]['variant']:
            allele_freq = exac_annotations[line]['variant']['allele_freq']
        else:
            allele_freq = 'NA'
        exac_allele_freq[position] = allele_freq

print(" Done with Exac API runs.")
print(" Now let's parse the allele frequencies requested.\n")
        
# Check the total number of variants taking into account multiallelic variants as well.
# print(len(exac_allele_freq))

######################################################################################
#
print(" Now we are ready to compute the following required calculations")
print(" 1. Type of variation (Substitution, Insertion, Silent, Intergenic, etc.) If there are multiple possibilities, annotate with the most deleterious possibility.")
print(" 2. Depth of sequence coverage at the site of variation.")
print(" 3. Number of reads supporting the variant.")
print(" 4. Percentage of reads supporting the variant versus those supporting reference reads.")
print(" 5. Allele frequency of variant from Broad Institute ExAC Project API\n")
#
######################################################################################

print("   As explained in http://snpeff.sourceforge.net/SnpEff_manual.html,")
print("   When multiple effects are available, they are sorted first by Effect_Impact in SNpEff\n")

orig_stdout = sys.stdout
f = open('Challenge_data_snpeff_modified.vcf', 'w')
sys.stdout = f

VCF = open("Challenge_data.snpeff.vcf", "r")
countf = 0 # to grab the format as a dict (Could vary depending on file)
for line in VCF:
    if line.startswith("#") and not line.startswith("#CHROM"):
        print (line, end='')
    elif line.startswith("#CHROM"):
        fields = line.strip("# \n").split("\t")
        print (line.strip() + "\t" + "NumAltAlleles" + "\t" + "MostDeleteriousPerVariant" + "\t" + "EXAC_AF_PERALLELE" + "\t" + "normalResults" + "\t" + "vaf5Results")
    elif not line.startswith("#"):
        variables = line.strip().split("\t")
        toks = dict(zip(fields, variables))

        # 1. The most deleterious impact on a per allele basis is extracted into a separate column.
        #    * SnpEff sorts the variants with most deleterious impact to appear first on a per allele bases.
        #    * Result is extracted into a separate column to highlight taking into account multiallelic variants.
        #    * Need to loop over mulltiallelic variants
        snpeffimp = toks["INFO"].split("ANN=")[1].split(";")[0]
        multialleles = toks["ALT"].split(",")
        if len(multialleles) > 1:
            snpeffimp2 = ",".join(snpeffimp.split(",")[0:len(multialleles)])
#            print("snpeffimp2=",snpeffimp2)            
        else :
            snpeffimp2 = snpeffimp.split(",")[0]        

################################################################################################
# Our challenge vcf contains a normal sample and perhaps an associated diseased sample annotated as vaf5.
#        2. FORMAT/DP is the depth of sequence coverage at the site of variation.
#        3. FORMAT/AO reports the number of reads supporting the variant.
#           FORMAT/RO reports the number of reads supporting the reference.
#        4. FORMAT/AO divided by FORMAT/DP to provide percentage of reads supporting the variant versus those supporting reference reads.
#           This ratio will provide the fraction of reads supporting the variants and not the reference reads. It should be 1 where RO=0.
#        5. Allele frequency of variant from Broad Institute ExAC Project API
#           Now we can easily extract allele frequency stored in a list object from BULK query ran earlier in the code.
#        The questions do not explicitly mention the calculations as "Overall" or on a "Per Sample" basis.
#           We are reportign calculations on a per sample basis.
#        In future, the script can be automated to run over all the samples. Need additional time to add extra functionalities.
#           Currently, decided to code on a per sample basis given the time restriction.
#
################################################################################################
        # Some definitions and initiations for calculating on a per sample bases.
        countallelesnormal = [];
        countallelesvaf5 = [];
        pctsuppvnormal = [];
        pctsuppvvaf5 = [];
        afexac = [];
        if countf == 0:
            formtoks = toks["FORMAT"].split(":")
            countf += 1
        #create dicts for computing counts of alleles on a per sample bases
        samplenormal = dict(zip(formtoks,toks["normal"].split(":")))
        samplevaf5 = dict(zip(formtoks,toks["vaf5"].split(":")))
        #reference allele counts
        normaltotaldepth = float(samplenormal["DP"])
        vaf5totaldepth = float(samplevaf5["DP"])
        #alternate allele counts
        samplenao = samplenormal["AO"].split(",")
        samplevao = samplevaf5["AO"].split(",")

        # Let's not forget to loop over the multiallelic variants.
        multialleles = toks["ALT"].split(",")
        for i in range(0,len(multialleles)):
            # Need the read information for each allele per sample
            nalt=float(samplenao[i]);
            countallelesnormal.append(str(nalt))
            pctsuppvnormal.append("%.4f" % (nalt/normaltotaldepth))

            # repeat for vaf5. 
            valt=float(samplevao[i])
            countallelesvaf5.append(str(valt))
            pctsuppvvaf5.append("%.4f" % (valt/vaf5totaldepth))
            location = toks["CHROM"]+"-"+toks["POS"]+"-"+toks["REF"]+"-"+multialleles[i]
            av=exac_allele_freq[location]
            if av == "NA":
                av = "NA"
            else:
                av = format(float(av), '.4f')
            afexac.append(av)

        # reporting the total number of alleles per variant can be very helpful for QCing multiallelic rows
        # snpeffimp2 variable reports the most deleterious impact on a per allele bases
        print(line.strip() + "\t" + str(len(multialleles)) + "\t" + snpeffimp2 + "\t" + \
              ",".join(afexac) + "\t" + \
              ":".join(["TotalDepth="+str(normaltotaldepth),"AllelicCounts="+",".join(countallelesnormal),"PctAlleleVsTotal="+",".join(pctsuppvnormal)]) + "\t" + \
            ":".join(["TotalDepth="+str(vaf5totaldepth), "AllelicCounts="+",".join(countallelesvaf5), "PctAlleleVsTotal="+",".join(pctsuppvvaf5)]))

sys.stdout = orig_stdout
f.close()

################################################################################################################
#
print("Congratulations, the output file Challenge_data_snpeff_modified.vcf is now available for your review")
#
################################################################################################################
print("The Program took --- %s seconds ---from start to finish" % (time.time() - start_time))
#
################################################################################################################
