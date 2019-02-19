
import urllib
import os
import csv
from Bio import Entrez
import math
from openpyxl import load_workbook
import urllib2
from urllib2 import URLError, HTTPError
from Bio.Blast.Applications import NcbiblastxCommandline
import xml.etree.ElementTree as ET
import sys

from sklearn.metrics import roc_curve, auc, precision_recall_curve
from sklearn.metrics import matthews_corrcoef, confusion_matrix

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

"""

ta Contig_Acc_version (extrahera ut hela contig)
sen ta 100 baser runtom

sen ta start, end
get sekvens

blast mot kromosom (enklare?) / gen
titta om vi get _en_ hit,

sen har vi poisitionen i next field


"""



def check_mvp():

    matches = 0
    sequence = ''
    prev_chr = -1
    tot_var = 0

    with open('/scratch/MVP_scores_hg19.txt_old', 'r') as file:

        file.readline()

        for line in file:

            line = line.strip()
            fields = line.split('\t')

            chr_no = fields[0]
            position = fields[1]
            ref = fields[4]
            alt = fields[5]

            nuc_ref = fields[2]

            variant_line = "%s:%s:%s:%s" % (chr_no, position, ref, alt)
            #print(variant_line)

            if prev_chr != chr_no:

                print('On chr %s' % chr_no)
                print('%d / %d' % (matches, tot_var))
                with open('/scratch/hg19/chr' + chr_no + '.fa') as hg19f:
                    hg19f.readline()
                    sequence = hg19f.read()
                    sequence = sequence.strip()
                    sequence = sequence.replace('\n', '')

            
            position = int(position)
            #print(position)
            #print(nuc_ref)
            #print(sequence[position-1])

            if nuc_ref == sequence[position-1]:
                matches += 1

            tot_var += 1

            prev_chr = chr_no

    print('%d / %d' % (matches, tot_var))




def blast_compare_path():

    Entrez.email = "alexander.kvist@outlook.com"

    pathogenic_variants = []
    prev_id = ''
    sequence = ''
    p_handle = None
    total_match = 0
    total_variants = 0
    #go through pathogenic test variants
    with open('../../Pipeline/Data/Pathogenic_PON-P2_test_data.csv', 'r') as csvfile:
        #skip first row
        csvfile.readline()
        creader = csv.reader(csvfile, delimiter='\t')
        for row in creader:


            genomic_field = row[6]
            hgsv_field = row[8]

            if not genomic_field or not hgsv_field:
                print('Skip')
                continue

            chr_no = genomic_field.split('.')[0]
            chr_no = int(chr_no[3:])
            chr_no = str(chr_no)
            if chr_no == "23": chr_no = 'X'
            position = genomic_field.split('.')[-1]
            position = position[0:-2]

            ref = hgsv_field[2]
            alt = hgsv_field[-1]
            ref_nucl = genomic_field[-2]

            variant_line = "%s:%s:%s:%s:%s" % (chr_no, position, ref, alt, ref_nucl)
            
            strand = row[-3]
            print(strand)

            #ok we need the sequence pointed to by the genomic identifier
            genomic_id = genomic_field.split(':')[0]

            #we need to download these into some folder
            genomic_folder = "/scratch/p2_genomic/"

            #ok open the genomic seq
            with open(genomic_folder + genomic_id + '.fasta') as file:
                file.readline()
                sequence = file.read()
                sequence = sequence.strip()
                sequence = sequence.replace('\n', '')


            position = int(position)
            #then we want the sequence +-100 from the mutation position
            sequence_window = sequence[position-100:position+100]

            #save the sequence to a temporary file..
            tmp_file = 'TEMP_FILE.fasta'
            with open(tmp_file, 'w') as f:
                f.write('>\n%s\n' % sequence_window)

            print('blasting..')

            fasta_query = tmp_file
            out_name = 'TMP_out.txt'
            blastx_cline = NcbiblastxCommandline(cmd="blastn", query=fasta_query, db="/scratch/hg19/blastdb/chr" + chr_no + ".fa", evalue=0.001, 
                                                 outfmt=7, out=out_name)
            print(blastx_cline)
            stdout, stderr = blastx_cline()

            #then we open the result.. we want fields 8 and 9
            hits = []
            with open(out_name, 'r') as res_f:
                for line in res_f:
                    if not line.startswith('#'):
                        hits.append(line.strip())

            if len(hits) > 1:
                print('More than one hit, oops?')

            if not hits:
                print('No hits?')
                continue

            sbjct_start = hits[0].split('\t')[8]
            print('Subject start: ' + sbjct_start)

            blast_pos = int(sbjct_start) + 99

            #now we need to grab the MVP prediction for this chromosome and that blast_pos
            seek_line = chr_no + '\t' + str(blast_pos)
            mvp_hits = []
            with open('/scratch/MVP_scores_hg19.txt_old', 'r') as file:

                file.readline()

                for line in file:
                    if not line.startswith(seek_line) and len(mvp_hits) > 0:
                        break
                    elif not line.startswith(seek_line):
                        continue
                    else:
                        print('Found MVP prediction')
                        mvp_hits.append(line)
                        

            #set up a fasta like format? with variant as > and MVP hits as lines
            variant_line += ":" + genomic_id
            with open('mvp_p2_pathogenic.fasta', 'a') as of:

                print('Variant is \n' + variant_line)
                of.write(">%s:%s\n" % (variant_line, 'pathogenic'))
                for mvp_line in mvp_hits:
                    mvp_line.strip()
                    mvp_chr_no = mvp_line.split('\t')[0]
                    mvp_position = mvp_line.split('\t')[1]
                    mvp_ref = mvp_line.split('\t')[4]
                    mvp_alt = mvp_line.split('\t')[5]
                    mvp_ref_nucl = mvp_line.split('\t')[2]
                    mvp_score = mvp_line.split('\t')[-1]
                    mvp_variant_line = "%s:%s:%s:%s:%s:%s" % (mvp_chr_no, mvp_position, mvp_ref, mvp_alt, mvp_ref_nucl, mvp_score)
                    print('MVP prediction is \n' + mvp_variant_line)
                    of.write(mvp_line)


            print('\n')


def get_mvp_predictions():

    #we should set up a file with the variant line and the MVP prediction scores appended at the end
    mvp_prediction_file = 'mvp_predictions.txt'

    Entrez.email = "alexander.kvist@outlook.com"

    prev_id = ''
    sequence = ''
    p_handle = None
    total_match = 0
    total_variants = 0

    
    #go through pathogenic
    #go through neutral test variants saved in file
    with open('p2_variants_mvp_running_path.txt', 'r') as file:

        for line in file:

            #skip empty lines
            if len(line) <= 1:
                continue

            skip_variant = False

            print(line)
            fields = line.split(':')

            genomic_id = fields[1]
            chr_no = fields[0]
            #23 is X
            if chr_no == "23": chr_no = "X"
            position = fields[2]
            ref = fields[3]
            alt = fields[4]
            ref_nucl = fields[5]
            variant_line = "%s:%s:%s:%s:%s" % (chr_no, position, ref, alt, ref_nucl)

            #let's see if we have the variant (without chr no and ref nucl etc in the file already..)
            #r_seek_line = "%s:%s:%s" % (chr_no, position, ref, alt, ref_nucl)
            with open(mvp_prediction_file, 'r') as r_f:
                for line in r_f:
                    line = line.strip()
                    if variant_line in line:
                        print('Already have; next')
                        skip_variant = True

            if skip_variant:
                continue

            #we need to download these into some folder
            genomic_folder = "/scratch/p2_genomic/"

            #ok open the genomic seq
            with open(genomic_folder + genomic_id + '.fasta') as s_f:
                s_f.readline()
                sequence = s_f.read()
                sequence = sequence.strip()
                sequence = sequence.replace('\n', '')


            position = int(position)
            #then we want the sequence +-100 from the mutation position
            sequence_window = sequence[position-100:position+100]

            #save the sequence to a temporary file..
            tmp_file = 'TEMP_FILE.fasta'
            with open(tmp_file, 'w') as f:
                f.write('>\n%s\n' % sequence_window)

            print('blasting..')

            fasta_query = tmp_file
            out_name = 'TMP_out.txt'
            blastx_cline = NcbiblastxCommandline(cmd="blastn", query=fasta_query, db="/scratch/hg19/blastdb/chr" + chr_no + ".fa", evalue=0.001, 
                                                 outfmt=7, out=out_name)
            print(blastx_cline)
            stdout, stderr = blastx_cline()

            #then we open the result.. we want fields 8 and 9
            hits = []
            with open(out_name, 'r') as res_f:
                for line in res_f:
                    if not line.startswith('#'):
                        hits.append(line.strip())

            if len(hits) > 1:
                print('More than one hit, oops?')

            if not hits:
                print('No hits?')
                continue

            sbjct_start = hits[0].split('\t')[8]
            print('Subject start: ' + sbjct_start)

            blast_pos = int(sbjct_start) + 99

            #now we need to grab the MVP prediction for this chromosome and that blast_pos
            seek_line = chr_no + '\t' + str(blast_pos)
            mvp_hits = []
            with open('/scratch/MVP_scores_hg19.txt_old', 'r') as file:

                file.readline()

                for line in file:
                    if not line.startswith(seek_line) and len(mvp_hits) > 0:
                        break
                    elif not line.startswith(seek_line):
                        continue
                    else:
                        print('Found MVP prediction')
                        mvp_hits.append(line)
                        

            match_mvp_prediction = 'missing'
            print('Variant is \n' + variant_line)
            for mvp_line in mvp_hits:
                mvp_chr_no = mvp_line.split('\t')[0]
                mvp_position = mvp_line.split('\t')[1]
                mvp_ref = mvp_line.split('\t')[4]
                mvp_alt = mvp_line.split('\t')[5]
                mvp_ref_nucl = mvp_line.split('\t')[2]
                mvp_variant_line = "%s:%s:%s:%s:%s" % (mvp_chr_no, mvp_position, mvp_ref, mvp_alt, mvp_ref_nucl)
                print('MVP prediction is \n' + mvp_variant_line)

                #we only want the mutation with the same ref and alt AA and same ref nucleotide
                if mvp_ref == ref and mvp_alt == alt and mvp_ref_nucl == ref_nucl:
                    match_mvp_prediction = mvp_line.split('\t')[-2] + ':' + mvp_line.split('\t')[-1]
                    print('Matching MVP prediction is %s' % match_mvp_prediction.strip())

            #append prediction along with true label at end
            with open(mvp_prediction_file, 'a') as pred_f:
                pred_f.write(variant_line + ':' + match_mvp_prediction.strip() + ':pathogenic\n')

            print('\n')


    #go through neutral test variants saved in file
    with open('p2_variants_mvp_running_neutral.txt', 'r') as file:

        for line in file:

            #skip empty lines
            if len(line) <= 1:
                continue

            skip_variant = False

            print(line)
            fields = line.split(':')

            genomic_id = fields[1]
            chr_no = fields[0]
            position = fields[2]
            ref = fields[3]
            alt = fields[4]
            ref_nucl = fields[5]
            variant_line = "%s:%s:%s:%s:%s" % (chr_no, position, ref, alt, ref_nucl)

            with open(mvp_prediction_file, 'r') as r_f:
                for line in r_f:
                    line = line.strip()
                    if variant_line in line:
                        print('Already have; next')
                        skip_variant = True

            if skip_variant:
                continue

            #we need to download these into some folder
            genomic_folder = "/scratch/p2_genomic/"

            #ok open the genomic seq
            with open(genomic_folder + genomic_id + '.fasta') as s_f:
                s_f.readline()
                sequence = s_f.read()
                sequence = sequence.strip()
                sequence = sequence.replace('\n', '')


            position = int(position)
            #then we want the sequence +-100 from the mutation position
            sequence_window = sequence[position-100:position+100]

            #save the sequence to a temporary file..
            tmp_file = 'TEMP_FILE.fasta'
            with open(tmp_file, 'w') as f:
                f.write('>\n%s\n' % sequence_window)

            print('blasting..')

            fasta_query = tmp_file
            out_name = 'TMP_out.txt'
            blastx_cline = NcbiblastxCommandline(cmd="blastn", query=fasta_query, db="/scratch/hg19/blastdb/chr" + chr_no + ".fa", evalue=0.001, 
                                                 outfmt=7, out=out_name)
            print(blastx_cline)
            stdout, stderr = blastx_cline()

            #then we open the result.. we want fields 8 and 9
            hits = []
            with open(out_name, 'r') as res_f:
                for line in res_f:
                    if not line.startswith('#'):
                        hits.append(line.strip())

            if len(hits) > 1:
                print('More than one hit, oops?')

            if not hits:
                print('No hits?')
                continue

            sbjct_start = hits[0].split('\t')[8]
            print('Subject start: ' + sbjct_start)

            blast_pos = int(sbjct_start) + 99

            #now we need to grab the MVP prediction for this chromosome and that blast_pos
            seek_line = chr_no + '\t' + str(blast_pos)
            mvp_hits = []
            with open('/scratch/MVP_scores_hg19.txt_old', 'r') as file:

                file.readline()

                for line in file:
                    if not line.startswith(seek_line) and len(mvp_hits) > 0:
                        break
                    elif not line.startswith(seek_line):
                        continue
                    else:
                        print('Found MVP prediction')
                        mvp_hits.append(line)
                        

            match_mvp_prediction = 'missing'
            print('Variant is \n' + variant_line)
            for mvp_line in mvp_hits:
                mvp_chr_no = mvp_line.split('\t')[0]
                mvp_position = mvp_line.split('\t')[1]
                mvp_ref = mvp_line.split('\t')[4]
                mvp_alt = mvp_line.split('\t')[5]
                mvp_ref_nucl = mvp_line.split('\t')[2]
                mvp_variant_line = "%s:%s:%s:%s:%s" % (mvp_chr_no, mvp_position, mvp_ref, mvp_alt, mvp_ref_nucl)
                print('MVP prediction is \n' + mvp_variant_line)

                #we only want the mutation with the same ref and alt AA and same ref nucleotide
                if mvp_ref == ref and mvp_alt == alt and mvp_ref_nucl == ref_nucl:
                    match_mvp_prediction = mvp_line.split('\t')[-2] + ':' + mvp_line.split('\t')[-1]
                    print('Matching MVP prediction is %s' % match_mvp_prediction.strip())

            #append prediction along with true label at end
            with open(mvp_prediction_file, 'a') as pred_f:
                pred_f.write(variant_line + ':' + match_mvp_prediction.strip() + ':neutral\n')

            print('\n')





    print('DONE')



def map_variants():


    all_variants = []

    #go through neutral variants
    neutral_variants = []
    prev_id = ''
    sequence = ''
    p_handle = None
    total_match_neutral = 0
    total_variants_neutral = 0
    chr_no = -1
    inversion = False

    #can we read variants and continue at end?
    running_file = 'p2_variants_mvp_running_neutral.txt'
    skip_variant = False

    #go through pathogenic test variants
    with open('../Pipeline/Data/Neutral_PON-P2_test_data.csv', 'r') as csvfile:
        #skip first row
        csvfile.readline()
        creader = csv.reader(csvfile, delimiter='\t')
        for row in creader:

            skip_variant = False

            genomic_id = row[1]
            hgsv_field = row[21]
            position_field = row[2]     #this is 0-based

            if not genomic_id or not hgsv_field or  len(position_field) <= 1:
                print('Skip')
                continue

            
            position = int(position_field) + 1  #make it 1-based

            ref = hgsv_field[2]
            alt = hgsv_field[-1]
            ref_nucl = row[11]

            


            #let's see if we have the variant (without chr no and ref nucl etc in the file already..)
            r_seek_line = "%s:%s:%s:%s" % (genomic_id, position, ref, alt)
            with open(running_file, 'r') as r_f:
                for line in r_f:
                    if r_seek_line in line:
                        print('Already have; next')
                        skip_variant = True

            if skip_variant:
                continue

            if genomic_id != prev_id:

                #we want to find the chromosome number
                url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=" + genomic_id + "&rettype=native&retmode=xml"
                try:
                    response = urllib2.urlopen(url)
                except HTTPError as e:
                    print('NCBI: The server could not fulfill the request.')
                    print('NCBI: Error code: ', e.code)
                    GOresult = None
                    sys.exit(0)
                except URLError as e:
                    print('NCBI: We failed to reach a server.')
                    print('NCBI: Reason: ', e.reason)
                    GOresult = None
                    sys.exit(0)
                else:
                    result = response.read()
                    response.close()

                with open('TEMPxml.xml', 'w') as x_f:
                    x_f.write(result)


                tree = ET.parse('TEMPxml.xml')
                root = tree.getroot()
                for child in root.iter('Seqdesc_title'):
                    print(child.text)
                    #let's hope it's always as easy as to grab fourth field on split whitespace
                    chr_no = child.text.split(' ')[3]

            print('Chr no: ' + chr_no)
            variant_line = "%s:%s:%s:%s:%s:%s:%s" % (chr_no, genomic_id, position, ref, alt, ref_nucl, 'neutral')
            print(variant_line)


            #we need to download these into some folder
            genomic_folder = "/scratch/p2_genomic/"
            if not os.path.exists(genomic_folder + genomic_id + '.fasta'):

                #we can try instead something like 
                #https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_000005.9&rettype=fasta&retmode=text
                print('Downloading ' + genomic_id)
                url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=" + genomic_id + "&rettype=fasta&retmode=text"

                try:
                    response = urllib2.urlopen(url)
                except HTTPError as e:
                    print('NCBI: The server could not fulfill the request.')
                    print('NCBI: Error code: ', e.code)
                    GOresult = None
                    sys.exit(0)
                except URLError as e:
                    print('NCBI: We failed to reach a server.')
                    print('NCBI: Reason: ', e.reason)
                    GOresult = None
                    sys.exit(0)
                else:
                    result = response.read()
                    response.close()

                with open(genomic_folder + genomic_id + '.fasta', 'w') as out:
                    out.write(result)

            #ok open the genomic seq
            with open(genomic_folder + genomic_id + '.fasta', 'r') as file:
                file.readline()
                sequence = file.read()
                sequence = sequence.strip()
                sequence = sequence.replace('\n', '')


            position = int(position)

            #So it looks like some reference nucleotides refer to the -strand, but otherwise everything is a match. So if not match, swap for complement, and keep that.
            #check for inversion only once, and then keep it like that for this id
            if genomic_id != prev_id:
                if ref_nucl != sequence[position-1]: 
                    print('Seems like -strand')
                    inversion = True
                else:
                    inversion = False

            if inversion:
                if ref_nucl == 'A':
                    ref_nucl = 'T'
                elif ref_nucl == 'T':
                    ref_nucl = 'A'
                elif ref_nucl == 'G':
                    ref_nucl = 'C'
                elif ref_nucl == 'C':
                    ref_nucl = 'G'


            print('Ref: %s, on position: %s, match: %s' % (ref_nucl, sequence[position-1], str(ref_nucl == sequence[position-1])))
            if ref_nucl == sequence[position-1]:
                total_match_neutral += 1
            total_variants_neutral += 1

            variant_line = "%s:%s:%s:%s:%s:%s:%s:%s" % (chr_no, genomic_id, position, ref, alt, ref_nucl, 'neutral', str(inversion))
            print(variant_line)
            all_variants.append(variant_line)
            #write as we go in case of errors
            with open(running_file, 'a') as r_f:
                r_f.write(variant_line + '\n')

            prev_id = genomic_id


    print('Found ' + str(total_match_neutral) + ' out of ' + str(total_variants_neutral))

    #first we go through test variants and grab the contig_acc_version

    Entrez.email = "alexander.kvist@outlook.com"

    pathogenic_variants = []
    prev_id = ''
    sequence = ''
    p_handle = None
    total_match_path = 0
    total_variants_path = 0
    chr_no = -1

    #can we read variants and continue at end?
    running_file = 'p2_variants_mvp_running_path.txt'
    skip_variant = False

    #go through pathogenic test variants
    with open('../Pipeline/Data/Pathogenic_PON-P2_test_data.csv', 'r') as csvfile:
        #skip first row
        csvfile.readline()
        creader = csv.reader(csvfile, delimiter='\t')
        for row in creader:


            genomic_field = row[6]
            hgsv_field = row[8]

            if not genomic_field or not hgsv_field:
                print('Skip')
                continue

            chr_no = genomic_field.split('.')[0]
            chr_no = int(chr_no[3:])
            chr_no = str(chr_no)
            position = genomic_field.split('.')[-1]
            position = position[0:-2]

            ref = hgsv_field[2]
            alt = hgsv_field[-1]
            ref_nucl = genomic_field[-2]

            #ok we need the sequence pointed to by the genomic identifier
            genomic_id = genomic_field.split(':')[0]
            #print(genomic_id)

            #let's see if we have the variant (without chr no and ref nucl etc in the file already..)
            r_seek_line = "%s:%s:%s:%s" % (genomic_id, position, ref, alt)
            with open(running_file, 'r') as r_f:
                for line in r_f:
                    if r_seek_line in line:
                        print('Already have; next')
                        skip_variant = True

            if skip_variant:
                continue

            
            

            variant_line = "%s:%s:%s:%s:%s:%s:%s" % (chr_no, genomic_id, position, ref, alt, ref_nucl, 'pathogenic')
            print(variant_line)
            strand = row[-3]
            #print(strand)


            #we need to download these into some folder
            genomic_folder = "/scratch/p2_genomic/"
            if not os.path.exists(genomic_folder + genomic_id + '.fasta'):

                #we can try instead something like 
                #https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_000005.9&rettype=fasta&retmode=text
                print('Downloading ' + genomic_id)
                url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=" + genomic_id + "&rettype=fasta&retmode=text"

                try:
                    response = urllib2.urlopen(url)
                except HTTPError as e:
                    print('NCBI: The server could not fulfill the request.')
                    print('NCBI: Error code: ', e.code)
                    GOresult = None
                    sys.exit(0)
                except URLError as e:
                    print('NCBI: We failed to reach a server.')
                    print('NCBI: Reason: ', e.reason)
                    GOresult = None
                    sys.exit(0)
                else:
                    result = response.read()
                    response.close()

                with open(genomic_folder + genomic_id + '.fasta', 'w') as out:
                    out.write(result)

            #ok open the genomic seq
            with open(genomic_folder + genomic_id + '.fasta') as file:
                file.readline()
                sequence = file.read()
                sequence = sequence.strip()
                sequence = sequence.replace('\n', '')


            position = int(position)
            #then we want the sequence +-100 from the mutation position
            sequence_window = sequence[position-100:position+100]


            print('Ref: %s, on position: %s, match: %s' % (ref_nucl, sequence[position-1], str(ref_nucl == sequence[position-1])))
            if ref_nucl == sequence[position-1]:
                total_match_path += 1
            total_variants_path += 1

            with open(running_file, 'a') as r_f:
                r_f.write(variant_line + '\n')

            all_variants.append(variant_line)


    print('Found ' + str(total_match_path) + ' out of ' + str(total_variants_path))


    #then write to file
    with open('p2_test_mvp_format', 'w') as f:
        for variant in all_variants:
            f.write(variant + '\n')

    print('Summary:')
    print('Found ' + str(total_match_neutral) + ' out of ' + str(total_variants_neutral) + ' neutral')
    print('Found ' + str(total_match_path) + ' out of ' + str(total_variants_path) + ' pathogenic')





def get_mvp_metrics():


    #so we open the MVP prediction file, grab the MVP sigmoid out (the -2 fields) and compare to true label (-1 field)
    #set up a vector of predictions and true labels, get ROC

    predictions = []
    labels = []
    tot_variants = 0
    tot_missing = 0
    with open('mvp_predictions.txt', 'r') as file:

        for line in file:

            tot_variants += 1
            line = line.strip()
            mvp_pred = line.split(':')[-2]
            if mvp_pred == 'missing':
                tot_missing += 1
                continue
            true_label = line.split(':')[-1]
            predictions.append(float(mvp_pred))
            if true_label == 'neutral':
                labels.append(0)
            else:
                labels.append(1)

    #print(predictions)
    #print(labels)
    predictions = np.array(predictions)
    labels = np.array(labels)

    print('Missing %d out of %d predictions' % (tot_missing, tot_variants))
    #then we get the ROC and AUC

    #Get test metrics
    c0_indices = labels == 0
    c1_indices = labels == 1
    num_neg = sum(c0_indices.astype(int))
    num_pos = sum(c1_indices.astype(int))

    #threshold is 0.5 to assign class
    predictions_copy = (predictions > 0.35).astype(int)

    #uhhh
    match = 0
    for i in range(0, len(labels)):
        if labels[i] == predictions_copy[i]: match += 1
    print(match)
    print(len(labels))
    print( float(match) / float(len(labels)))

    confmat = confusion_matrix(labels, predictions_copy)
    TN = confmat[0][0]
    FN = confmat[1][0]
    TP = confmat[1][1]
    FP = confmat[0][1]
    PPV = float(TP)/(TP+FP)
    NPV = float(TN)/(TN+FN)
    sens = float(TP)/(TP+FN)
    spec = float(TN)/(FP+TN)
    acc = float((TP+TN))/(TP+TN+FP+FN)
    mcc = matthews_corrcoef(labels, predictions_copy)
    nmcc = float((1+mcc))/2
    OPM = float(( (PPV+NPV)*(sens+spec)*(acc+nmcc) )) / 8
    bacc = ( float(TP) / num_pos + float(TN) / num_neg ) / 2
    accmetrics_test = [PPV,NPV,sens,spec,acc,mcc,OPM, bacc]
    colNames = ['PPV', 'NPV', 'Sens', 'Spec', 'Acc', 'MCC', 'OPM', 'BACC']
    df = pd.DataFrame([accmetrics_test], columns=colNames)
    print("##### On test set ######")
    print(df)

    fpr = dict()
    tpr = dict()
    roc_auc = dict()

    # Compute micro-average ROC curve and ROC area
    fpr["test"], tpr["test"], _ = roc_curve(labels, predictions, pos_label=1)
    roc_auc["test"] = auc(fpr["test"], tpr["test"])

    plt.figure()
    lw = 2
    plt.plot(fpr["test"], tpr["test"], color='darkorange',
             lw=lw, label='ROC test (area = %0.2f)' % roc_auc["test"])
    plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
    plt.grid()
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC')
    plt.legend(loc="lower right")
    plt.savefig('MVP_ROC')






def map_pathogenic_variants():
    #first we go through test variants and grab the contig_acc_version

    Entrez.email = "alexander.kvist@outlook.com"

    pathogenic_variants = []
    prev_id = ''
    sequence = ''
    p_handle = None
    total_match_path = 0
    total_variants_path = 0
    chr_no = -1

    #can we read variants and continue at end?
    running_file = 'p2_PATHVARS_mvp.txt'
    skip_variant = False

    #go through pathogenic test variants
    with open('../Pipeline/Data/Pathogenic_PON-P2_test_data.csv', 'r') as csvfile:
        #skip first row
        csvfile.readline()
        creader = csv.reader(csvfile, delimiter='\t')
        for row in creader:


            genomic_field = row[6]
            hgsv_field = row[8]

            if not genomic_field or not hgsv_field:
                print('Skip')
                continue

            chr_no = genomic_field.split('.')[0]
            chr_no = int(chr_no[3:])
            chr_no = str(chr_no)
            position = genomic_field.split('.')[-1]
            position = position[0:-2]

            ref = hgsv_field[2]
            alt = hgsv_field[-1]
            ref_nucl = genomic_field[-2]

            #ok we need the sequence pointed to by the genomic identifier
            genomic_id = genomic_field.split(':')[0]
            #print(genomic_id)

            #let's see if we have the variant (without chr no and ref nucl etc in the file already..)
            r_seek_line = "%s:%s:%s:%s" % (genomic_id, position, ref, alt)
            with open(running_file, 'r') as r_f:
                for line in r_f:
                    if r_seek_line in line:
                        print('Already have; next')
                        skip_variant = True

            if skip_variant:
                continue

            
            

            variant_line = "%s:%s:%s:%s:%s:%s:%s" % (chr_no, genomic_id, position, ref, alt, ref_nucl, 'pathogenic')
            print(variant_line)
            strand = row[-3]
            #print(strand)


            #we need to download these into some folder
            genomic_folder = "/scratch/p2_genomic/"
            if not os.path.exists(genomic_folder + genomic_id + '.fasta'):

                #we can try instead something like 
                #https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_000005.9&rettype=fasta&retmode=text
                print('Downloading ' + genomic_id)
                url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=" + genomic_id + "&rettype=fasta&retmode=text"

                try:
                    response = urllib2.urlopen(url)
                except HTTPError as e:
                    print('NCBI: The server could not fulfill the request.')
                    print('NCBI: Error code: ', e.code)
                    GOresult = None
                    sys.exit(0)
                except URLError as e:
                    print('NCBI: We failed to reach a server.')
                    print('NCBI: Reason: ', e.reason)
                    GOresult = None
                    sys.exit(0)
                else:
                    result = response.read()
                    response.close()

                with open(genomic_folder + genomic_id + '.fasta', 'w') as out:
                    out.write(result)

            #ok open the genomic seq
            with open(genomic_folder + genomic_id + '.fasta') as file:
                file.readline()
                sequence = file.read()
                sequence = sequence.strip()
                sequence = sequence.replace('\n', '')


            position = int(position)
            #then we want the sequence +-100 from the mutation position
            sequence_window = sequence[position-100:position+100]


            print('Ref: %s, on position: %s, match: %s' % (ref_nucl, sequence[position-1], str(ref_nucl == sequence[position-1])))
            if ref_nucl == sequence[position-1]:
                total_match_path += 1
            total_variants_path += 1

            with open(running_file, 'a') as r_f:
                r_f.write(variant_line + '\n')

            all_variants.append(variant_line)


    print('Found ' + str(total_match_path) + ' out of ' + str(total_variants_path))


    #then write to file
    with open('p2_test_mvp_format', 'w') as f:
        for variant in all_variants:
            f.write(variant + '\n')

    print('Summary:')
    print('Found ' + str(total_match_neutral) + ' out of ' + str(total_variants_neutral) + ' neutral')
    print('Found ' + str(total_match_path) + ' out of ' + str(total_variants_path) + ' pathogenic')




def map_id_varline():

    #actually -- open the mvp predictions, append the identifier and mutation in the format we're using at start
    map_handle = open("tested_variants.txt", "w")

    with open('mvp_predictions.txt', 'r') as pred_file:

        for line in pred_file:

            line = line.strip()
            #get position, ref aa, alt aa, make sure we only get one hit
            mvp_position = line.split(':')[1]
            mvp_ref_aa = line.split(':')[2]
            mvp_alt_aa = line.split(':')[3]
            mvp_ref_nucl = line.split(':')[4]
            label = line.split(':')[-1]

            if label == 'neutral':

                #go through neutral test variants check matches
                matches = 0
                matching_mutation = "missing"
                with open('../../Pipeline/Data/Neutral_PON-P2_test_data.csv', 'r') as csvfile:
                    #skip first row
                    csvfile.readline()
                    creader = csv.reader(csvfile, delimiter='\t')
                    for row in creader:

                        genomic_id = row[1]
                        hgsv_field = row[21]
                        position_field = row[2]     #this is 0-based
                        if not genomic_id or not hgsv_field or  len(position_field) <= 1:
                            #print('Skip')
                            continue
                        
                        position = int(position_field) + 1  #make it 1-based
                        position = str(position)
                        ref = hgsv_field[2]
                        alt = hgsv_field[-1]
                        ref_nucl = row[11]
                        refseq_field = row[14]

                        if mvp_position == position and mvp_ref_aa == ref and mvp_alt_aa == alt:
                            matches += 1
                            matching_mutation = refseq_field + ':' + hgsv_field[2:]

                        

                #print("%s : matches %d" % (line, matches))
                if matches > 1:
                    print("%s : MORE matches %d" % (line, matches))
                    sys.exit(0)
                else:
                    map_handle.write(matching_mutation + ":" + line + "\n")
                        
            elif label == 'pathogenic':


                #go through pathogenic test variants check matches
                matches = 0
                matching_mutation = "missing"
                with open('../../Pipeline/Data/Pathogenic_PON-P2_test_data.csv', 'r') as csvfile:
                    #skip first row
                    csvfile.readline()
                    creader = csv.reader(csvfile, delimiter='\t')
                    for row in creader:

                        genomic_field = row[6]
                        hgsv_field = row[8]

                        position_field = row[2]     #this is 0-based
                        if not genomic_field or not hgsv_field:
                            #print('Skip')
                            continue
                        
                        position = genomic_field.split('.')[-1]
                        position = position[0:-2]

                        ref = hgsv_field[2]
                        alt = hgsv_field[-1]
                        ref_nucl = genomic_field[-2]
                        uniprot_field = row[7]

                        if mvp_position == position and mvp_ref_aa == ref and mvp_alt_aa == alt:
                            matches += 1
                            matching_mutation = uniprot_field + ':' + hgsv_field[2:]

                        

                #print("%s : matches %d" % (line, matches))
                if matches > 1:
                    print("%s : MORE matches %d" % (line, matches))
                    sys.exit(0)
                else:
                    map_handle.write(matching_mutation + ":" + line + "\n")


if __name__ == '__main__':


    #blast_compare_path()
    #get_mvp_predictions()
    #get_mvp_metrics()

    map_id_varline()


    

    