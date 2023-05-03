#!/usr/bin/env python
# coding: utf-8

# In[ ]:


### Boas Pucker ###
### Shakunthala Natarajan ###
### shakunthalabiotechnology@gmail.com ###
__usage__= """
            python3 TwinFinder.py
            mandatory:
            --gff <path to gff3file of organism 1> <path to gff3 file of organism 2>
            --fasta <path to whole genome FASTA file of organism 1> <path to whole genome FASTA file of organism 2>
            -- out <path to output folder>
            --tmp <path to temporary folder>
            
            optional:
            --cpu <number of threads to run BLAST (integer)> <default - 10>
            --name <name of organism 1> <name of organism 2> <default - org1 org2>
            --similarity <minimum similarity value that must be there between duplicated genes> <type - integer or float> < default - 80>
            --continuity <mode of working of the script> <normal> or <fragmented> <default - normal>
            --flankq <number of flanking genes to be considered in the query organism> <type - integer> <default - 6>
            --flankr <number of flanking genes to be considered in the reference organism> <type - integer> <default - 6>
            --cutoff <minimum cutoff syntney ratio to be fulfilled for an identified tandem group of genes and their corresponding  
                      anchor gene in the reference organism to be considered as syntenic> <type - float> <between 0 and 1> <default - 0.5>
            """
######### imports #########
import re, os, sys, subprocess
import pandas as pd
import numpy 
from numpy import mean
from operator import itemgetter
from datetime import datetime
from collections import Counter
from matplotlib import pyplot
import copy
import seaborn as sns
##### end of imports #####

print('Welcome to TwinFinder!!! Your analysis will begin shortly...')

"""
Part 1: Code to obtain peptide and cds FASTA files from GFF3 and Whole genome FASTA files
"""

def load_sequences( multiple_fasta_file ):
    """! @brief load candidate gene IDs from file """
    sequences = {}
    with open( multiple_fasta_file ) as f:
        header = f.readline()[1:].strip()
        if " " in header:
            header = header.split(' ')[0]
            if "\t" in header:
                header = header.split('\t')[0]
        elif "\t" in header:
                header = header.split('\t')[0]
        seq = []
        line = f.readline()
        while line:
            if line[0] == '>':
                    sequences.update( { header: "".join( seq ).upper() } )
                    header = line.strip()[1:]
                    if " " in header:
                        header = header.split(' ')[0]
                        if "\t" in header:
                            header = header.split('\t')[0]
                    elif "\t" in header:
                        header = header.split('\t')[0]
                    seq = []
            else:
                seq.append( line.strip() )
            line = f.readline()
        sequences.update( { header: "".join( seq ).upper() } )
    return sequences


def load_transcript_information_from_gff3( gff3_input_file ):
    """! @brief load all transcript information from gff3 file """
    # --- load all data from file --- #
    information = []
    with open( gff3_input_file, "r" ) as f:
        line = f.readline()
        while line:
            if line[0] != '#':
                parts = line.strip().split('\t')
                if len( parts ) > 5:
                    if parts[2].upper() == 'CDS':
                        if len( parts[-1] ) > len( "Parent=" ):
                            if ";" in parts[-1]:
                                parent = False
                                subparts = parts[-1].split(';')
                                for subp in subparts:
                                    if "Parent=" in subp:
                                        parent = subp.replace( "Parent=", "" )
                                    elif "transcript_id " in subp:
                                        parent = subp.split('"')[1]
                                if parent:
                                    information.append( { 	'chr': parts[0],
                                                                        'start': int( parts[3] ),
                                                                        'end': int( parts[4] ),
                                                                        'orientation': parts[6],
                                                                        'parent': parent
                                                                    } )
                                else:
                                    print ("no parent detected - " + line)
                            else:
                                parent=False
                                if "Parent=" in parts[-1]:
                                    parent=str(parts[-1]).replace("Parent=","")
                                if parent:
                                    information.append( { 	'chr': parts[0],
                                                                        'start': int( parts[3] ),
                                                                        'end': int( parts[4] ),
                                                                        'orientation': parts[6],
                                                                        'parent': parent
                                                                    } )
                                print ("only one field - " + line)
            line = f.readline()
    # --- sort data by parent --- #
    sorted_data = {}
    for each in information:
        try:
            sorted_data[ each['parent'] ].append( each )
        except KeyError:
            sorted_data.update( { each['parent']: [ each ] } )
    final_data = []
    for key in sorted_data.keys():
        if sorted_data[ key ][0] ['orientation'] == '+':
            final_data.append( sorted( sorted_data[ key ], key=itemgetter('start') ) )
        else:
            final_data.append( sorted( sorted_data[ key ], key=itemgetter('start') )[::-1] )
    return final_data

def construct_CDS_file( transcript_info, CDS_file, assembly ):
    """! @brief construct file with all sequences for translation """
    with open( CDS_file, "w" ) as out:
        for transcript in transcript_info:
            seq = []
            revcomp_status = False
            if transcript[0]['orientation'] == '-':
                revcomp_status = True
            for part in transcript:
                if revcomp_status:
                    seq.append( revcomp( assembly[ part['chr'] ][ part['start']-1:part['end'] ] ) )
                else:
                    seq.append( assembly[ part['chr'] ][ part['start']-1:part['end'] ] )
            out.write( '>' + transcript[0]['parent'].replace("Parent=", "") + '\n' + "".join( seq ) + '\n' )


def revcomp( seq ):
    """! @brief constructs revcomp """
    new_seq = []
    dictionary = { 'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N' }
    for nt in seq:
        try:
            new_seq.append( dictionary[ nt ] )
        except KeyError:
            print ("ERROR: " + nt)
            new_seq.append( "N")
    return ''.join( new_seq[::-1] )


def translate( seq, genetic_code ):
    """! @brief translates the given nucleotide sequence into peptide and splits at each star (stop codon) """
    seq = seq.upper()
    peptide = []
    for i in range( int( len( seq ) / 3.0 ) ):
        codon = seq[i*3:i*3+3]
        try:
            peptide.append( genetic_code[ codon ] )
        except:
            peptide.append( "*" )
    peptide =  "".join( peptide )
    if sum( [ peptide[0] != "M", "*" in peptide[:-1] ] ) > 0:
        peptide2 = []
        for i in range( int( ( len( seq )-1 ) / 3.0 ) ):
            codon = seq[1+i*3:1+i*3+3]
            try:
                peptide2.append( genetic_code[ codon ] )
            except:
                peptide2.append( "*" )
        peptide2 =  "".join( peptide2 )
        if sum( [ peptide2[0] != "M", "*" in peptide2[:-1] ] ) > 0:
            peptide3 = []
            for i in range( int( ( len( seq )-2 ) / 3.0 ) ):
                codon = seq[2+i*3:2+i*3+3]
                try:
                    peptide3.append( genetic_code[ codon ] )
                except:
                    peptide3.append( "*" )
            peptide3 =  "".join( peptide3 )
            pep_options = []
            if '*' in peptide:
                pep_options.append( { 'seq': peptide.split('*')[0], 'stopps': peptide.count('*'), 'len': len( peptide.split('*')[0] ) } )
            else:
                pep_options.append( { 'seq': peptide, 'stopps': peptide.count('*'), 'len': len( peptide ) } )
            if "*" in peptide2:
                pep_options.append( { 'seq': peptide2.split('*')[0], 'stopps': peptide2.count('*'), 'len': len( peptide2.split('*')[0] ) } )
            else:
                pep_options.append( { 'seq': peptide2, 'stopps': peptide2.count('*'), 'len': len( peptide2 ) } )
            if "*" in peptide3:
                pep_options.append( { 'seq': peptide3.split('*')[0], 'stopps': peptide3.count('*'), 'len': len( peptide3.split('*')[0] ) } )
            else:
                pep_options.append( { 'seq': peptide3, 'stopps': peptide3.count('*'), 'len': len( peptide3 ) } )
            return sorted( pep_options, key=itemgetter( 'len' ) )[-1]['seq']
        else:
            pep_options = []
            if '*' in peptide:
                pep_options.append( { 'seq': peptide.split('*')[0], 'stopps': peptide.count('*'), 'len': len( peptide.split('*')[0] ) } )
            else:
                pep_options.append( { 'seq': peptide, 'stopps': peptide.count('*'), 'len': len( peptide ) } )
            if "*" in peptide2:
                pep_options.append( { 'seq': peptide2.split('*')[0], 'stopps': peptide2.count('*'), 'len': len( peptide2.split('*')[0] ) } )
            else:
                pep_options.append( { 'seq': peptide2, 'stopps': peptide2.count('*'), 'len': len( peptide2 ) } )
            return sorted( pep_options, key=itemgetter( 'len' ) )[-1]['seq']
    else:
        return peptide
    
def transeq( input_file, output_file, strict_start, strict_end ):
    """! @brief run translation of coding sequences """	
    genetic_code = {'CTT': 'L', 'ATG': 'M', 'AAG': 'K', 'AAA': 'K', 'ATC': 'I', 'AAC': 'N', 'ATA': 'I', 'AGG': 'R', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'ACA': 'T', 'AGA': 'R', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'ACG': 'T', 'CCG': 'P', 'AGT': 'S', 'CAG': 'Q', 'CAA': 'Q', 'CCC': 'P', 'TAG': '*', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CCA': 'P', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': '*', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'GAG': 'E', 'TCG': 'S', 'TTA': 'L', 'GAC': 'D', 'TCC': 'S', 'GAA': 'E', 'TCA': 'S', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'TTC': 'F', 'GTT': 'V', 'GCT': 'A', 'ACC': 'T', 'TTG': 'L', 'CGT': 'R', 'TAA': '*', 'CGC': 'R'}
    sequences = load_sequences( input_file )
    with open( output_file, "w" ) as out:
        keys = sequences.keys()
        for key in keys:
            seq = sequences[ key ]
            if len( seq ) > 9:
                peptide = translate( seq, genetic_code )
                out.write( '>' + key + '\n' + peptide + '\n' )
            else:
                print (key + " - too short!")
t1=datetime.now()
print("Peptide and CDS FASTA files are ready"+':'+'\t'+str(t1))
"""
Part 2: To perform Reciprocal blast using the generated peptide FASTA files of organisms 1 and 2 and arrive at 
a consolidated list of genes showing reciprocal blast hit
"""
#defining the output destinations for the forward and the reverse blast steps

def gene_multiplications (output_dir, tmp_dir, name1, name2, pep_output_file_1, pep_output_file_2, gff3_input_file_1, gff3_input_file_2, cpu, threshold, flank_org1, flank_org2, cutoff, continuity):
    blast_output=os.path.join(output_dir, name1+name2+'blast')
    os.makedirs(blast_output, exist_ok=True)
    result=os.path.join(output_dir, 'RESULTS')
    os.makedirs(result, exist_ok=True)
    no_trans_pep1=os.path.join(tmp_dir, name1+'_no_transcripts.pep.fasta')
    no_trans_pep2=os.path.join(tmp_dir, name2+'_no_transcripts.pep.fasta')
    fwdblast_out=os.path.join(blast_output, 'fwd-blast.tsv')
    fwdblast_besthits=os.path.join(blast_output, 'fwd-blast-results.tsv')
    database1=name1+"_db"
    database2=name2+"_db"
    orthologs=os.path.join(result, 'Orthologs_threshold_'+str(threshold)+'.tsv')
    self_blast=os.path.join(blast_output, name1+'_self_blast.tsv')
    multiplications=os.path.join(tmp_dir, name1+'_tandem_multiplications_temporary.tsv')
    tandem_multiplications=os.path.join(result, name1+'_highly_similar_tandem_multiplications_threshold_'+str(threshold)+'.tsv')
    non_syntelog_tandems=os.path.join(result, name1+'_non_syntelog_tandem_multiplications_threshold_'+str(threshold)+'.tsv')
    syntelog_tandem_file=os.path.join(result, name1+'_syntelog_tandem_multiplications_threshold_'+str(threshold)+'.tsv')
    histogram=os.path.join(result,name1+'_mutual_gene_similarities_histogram_'+str(threshold))
    import numpy
    
    
    """
    cleaning the peptide files to remove alternative transcripts
    """
    
    if not os.path.isfile( no_trans_pep1 ):
        transcripts_per_gene1 = {}
        with open( gff3_input_file_1, "r" ) as f:
            line = f.readline()
            while line:
                if line[0] != "#":
                    parts = line.strip().split('\t')
                    if len( parts ) > 2:
                        if (parts[2] == "transcript") or (parts[2] == "mRNA"):
                            partsnew=parts[-1].strip().split(';')
                            for each in partsnew:
                                if 'Parent=' in each:
                                    partsnew1=str(each).replace("Parent=", "")
                            for every in partsnew:
                                if 'ID=' in every:
                                    partsnew0=str(every).replace("ID=", "")
                            try:
                                transcripts_per_gene1[ partsnew1 ].append( partsnew0)
                            except KeyError:
                                transcripts_per_gene1.update( { partsnew1: [ partsnew0 ] } )
                line = f.readline()
        gene_names1 = transcripts_per_gene1.keys()
        with open( pep_output_file_1, "r" ) as f:
            with open( no_trans_pep1, "w" ) as out:
                line=f.readlines()
                pep_dict1={}
                for each in line:
                    if '>' in each:
                        ind=line.index(each)
                        pep_dict1.update({str(each).replace('>','').replace('\n',''):line[ind+1]})
                transcripts1=pep_dict1.keys()
                for gene in gene_names1:
                    trans_length1=[]
                    for trans in transcripts_per_gene1[gene]:
                        for each in transcripts1:
                            if trans==each:
                                if len(transcripts_per_gene1[gene])<2:
                                    out.write('>'+str(gene)+'\n'+str(pep_dict1[each]))
                                else:
                                    trans_length1.append(pep_dict1[each])
                    if len(trans_length1)==0:
                        pass
                    else:
                        seq = max(trans_length1, key=len)
                        out.write('>'+str(gene)+"\n"+str(seq))
    if not os.path.isfile( no_trans_pep2 ):
        transcripts_per_gene2 = {}
        with open( gff3_input_file_2, "r" ) as f:
            line = f.readline()
            while line:
                if line[0] != "#":
                    parts = line.strip().split('\t')
                    if len( parts ) > 2:
                        if (parts[2] == "transcript") or (parts[2] == "mRNA"):
                            partsnew=parts[-1].strip().split(';')
                            for each in partsnew:
                                if 'Parent=' in each:
                                    partsnew1=str(each).replace("Parent=", "")
                            for every in partsnew:
                                if 'ID=' in every:
                                    partsnew0=str(every).replace("ID=", "")
                            try:
                                transcripts_per_gene2[ partsnew1 ].append( partsnew0 )
                            except KeyError:
                                transcripts_per_gene2.update( { partsnew1: [ partsnew0] } )
                line = f.readline()
        gene_names2 = transcripts_per_gene2.keys()
        with open( pep_output_file_2, "r" ) as f:
            with open( no_trans_pep2, "w" ) as out:
                line=f.readlines()
                pep_dict2={}
                for each in line:
                    if '>' in each:
                        ind=line.index(each)
                        pep_dict2.update({str(each).replace('>','').replace('\n',''):line[ind+1]})
                transcripts2=pep_dict2.keys()
                for gene in gene_names2:
                    trans_length2=[]
                    for trans in transcripts_per_gene2[gene]:
                        for each in transcripts2:
                            if trans==each:
                                if len(transcripts_per_gene2[gene])<2:
                                    out.write('>'+str(gene)+'\n'+str(pep_dict2[each]))
                                else:
                                    trans_length2.append(pep_dict2[each])
                    if len(trans_length2)==0:
                        pass
                    else:
                        seq = max(trans_length2, key=len)
                        out.write('>'+str(gene)+"\n"+str(seq))
    t2=datetime.now()
    print("Peptide FASTA files without alternate transcripts are ready for BLASTp:"+'\t'+str(t2))
    """
    generating command line arguments for forward and self blast using blastp 
    and execution of blastp by interacting with the os
    """
    if not os.path.isfile( fwdblast_out ):
        #making blast database of organism 2
        pdata2=subprocess.Popen(args= "makeblastdb -in " + no_trans_pep2 + " -dbtype prot -parse_seqids" + " -out " + database2, shell=True)
        pdata2.communicate()
        pforward = subprocess.Popen( args= "blastp -query " + no_trans_pep1 + " -db " + database2 + " -evalue 1e-5" + " -out " + fwdblast_out + " -outfmt \'6 qseqid sseqid pident qcovs qlen slen length evalue bitscore\'" + " -num_threads " + str(cpu), shell=True )
        pforward.communicate()
        t3=datetime.now()
        print ("Forward blast done!:"+"\t"+str(t3))
    
    if not os.path.isfile( self_blast ):
        #making blast database of organism 1
        pdata1=subprocess.Popen(args= "makeblastdb -in " + no_trans_pep1 + " -dbtype prot -parse_seqids" + " -out " + database1, shell=True)
        pdata1.communicate()
        #self blast of organism 1
        pself = subprocess.Popen( args= "blastp -query " + no_trans_pep1 + " -db " + database1 + " -evalue 1e-5" + " -out " + self_blast + " -outfmt \'6 qseqid sseqid pident qcovs qlen slen length evalue bitscore\'" + " -num_threads " + str(cpu), shell=True )
        pself.communicate()
        tself=datetime.now()
        print ("Self blast of organism 1 done!:"+"\t"+str(tself))
        
    if not os.path.isfile( fwdblast_besthits):  
        #to get the best hits from forward blast 
    
        number_of_hits = 1
    
        fwd_best_hits = {}
        with open( fwdblast_out, "r" ) as f:
            line = f.readline()
            while line:
                parts = line.strip().split('\t')
                try:
                    fwd_best_hits[ parts[0] ].append( { 'score': float( parts[-1] ), 'line': line } )
                except:
                    fwd_best_hits.update( { parts[0]: [ { 'score': float( parts[-1] ), 'line': line } ] } )
                line = f.readline()
        fwd_final_sorted_hits = {}
        for key in sorted( fwd_best_hits.keys() ):
            data = sorted( fwd_best_hits[ key ], key=itemgetter('score') )[::-1]
            fwd_final_sorted_hits.update( { key: data } )
        with open( fwdblast_besthits, "w" ) as out:
            for key in sorted( fwd_final_sorted_hits.keys() ):
                if len( fwd_final_sorted_hits[ key ] ) > number_of_hits:
                    for each in fwd_final_sorted_hits[ key ][ : number_of_hits ]:
                        out.write( each[ 'line' ] )
                else:
                    for each in fwd_final_sorted_hits[ key ]:
                        out.write( each[ 'line' ] )
    
        t5=datetime.now()
        print ("Best hits from forward blast are ready:"+"\t"+str(t5))            
    """
    To correlate the genomic positions of the genes from the GFF3 file of organism 1 
    """
    
    #making a list of coding genes in organism1
    coding_genes1=[]
    with open(no_trans_pep1,'r') as f:
        line=f.readline()
        while line:
            if '>' in line:
                coding_genes1.append(str(line).replace('>','').replace('\n',''))
            line=f.readline()
    #making a list of coding genes in organism2
    coding_genes2=[]
    with open(no_trans_pep2,'r') as f:
        line=f.readline()
        while line:
            if '>' in line:
                coding_genes2.append(str(line).replace('>','').replace('\n',''))
            line=f.readline()
    #making dictionary of lists for organism1
    gene_pos_per_chr_1_unordered = {}
    genes_org1_unordered=[]
    with open( gff3_input_file_1, "r" ) as f:
        line = f.readline()
        while line:
            if line[0] != "#":
                parts = line.strip().split('\t')
                if len( parts ) > 2:
                    if parts[2] == "gene":
                        partsnew=parts[-1].strip().split(';')
                        for each in partsnew:
                            if 'ID=' in each:
                                partsnew1=str(each).replace("ID=","")
                
                        genes_org1_unordered.append(partsnew1)
                        try:
                            gene_pos_per_chr_1_unordered[ parts[0] ].append( partsnew1 )
                        except KeyError:
                            gene_pos_per_chr_1_unordered.update( { parts[0]: [ partsnew1 ] } )
            line = f.readline()
    contig_names_1 = gene_pos_per_chr_1_unordered.keys()
    
    #sorting the dictionary of lists to arrange the genes in each sequence
    gene_pos_per_chr_1_ordered={}
    for each in contig_names_1:
        gene_pos_per_chr_1_ordered.update({each:sorted(gene_pos_per_chr_1_unordered[each],key=lambda x: int("".join([i for i in x if i.isdigit()])))})
        
    #creating a copy of the ordered dictionary of lists
    gene_pos_per_chr_1=copy.deepcopy(gene_pos_per_chr_1_ordered)    
        
    #sorting the gene list in org1 and creating its copy
    genes_org1_ordered=sorted(genes_org1_unordered,key=lambda x: int("".join([i for i in x if i.isdigit()])))
    genes_org1=copy.deepcopy(genes_org1_ordered)
    
    #removing non-coding genes from dictionary of lists
    for each in contig_names_1:
        for every in gene_pos_per_chr_1_ordered[each]:
            if every not in coding_genes1:
                gene_pos_per_chr_1[each].remove(every)
                
    #removing non-coding genes from list of coding genes in organism1
    for one in genes_org1_ordered:
        if one not in coding_genes1:
            genes_org1.remove(one)

    """
    To correlate the genomic positions of the genes from the GFF3 file of organism 2
    """
    #making dictionary of lists for organism2
    gene_pos_per_chr_2_unordered = {}
    genes_org2_unordered=[]
    with open( gff3_input_file_2, "r" ) as f:
        line = f.readline()
        while line:
            if line[0] != "#":
                parts = line.strip().split('\t')
                if len( parts ) > 2:
                    if parts[2] == "gene":
                        partsnew=parts[-1].strip().split(';')
                        for each in partsnew:
                            if 'ID=' in each:
                                partsnew1=str(each).replace("ID=","")
                        genes_org2_unordered.append(partsnew1)
                        try:
                            gene_pos_per_chr_2_unordered[ parts[0] ].append( partsnew1 )
                        except KeyError:
                            gene_pos_per_chr_2_unordered.update( { parts[0]: [ partsnew1 ] } )
            line = f.readline()
    contig_names_2 = gene_pos_per_chr_2_unordered.keys()
    
    #sorting the dictionary of lists to arrange the genes in each sequence
    gene_pos_per_chr_2_ordered={}
    for each in contig_names_2:
        gene_pos_per_chr_2_ordered.update({each:sorted(gene_pos_per_chr_2_unordered[each],key=lambda x: int("".join([i for i in x if i.isdigit()])))})
        
    #creating copy of the ordered dictionary of lists
    gene_pos_per_chr_2=copy.deepcopy(gene_pos_per_chr_2_ordered)
        
    #sorting the gene list in org2 and creating its copy
    genes_org2_ordered=sorted(genes_org2_unordered,key=lambda x: int("".join([i for i in x if i.isdigit()])))
    genes_org2=copy.deepcopy(genes_org2_ordered)
    
    #removing non-coding genes from dictionary of lists
    for each in contig_names_2:
        for every in gene_pos_per_chr_2_ordered[each]:
            if every not in coding_genes2:
                gene_pos_per_chr_2[each].remove(every)
                
    #removing non-coding genes from list of coding genes in organism2
    for one in genes_org2_ordered:
        if one not in coding_genes2:
            genes_org2.remove(one)
    
    #making dictionary of forward blast best hits
    fwdblastbest_dic={}
    fwdblastbest_dic_rev=[]
    fwdblastbest_dic_both=[]
    with open(fwdblast_besthits,'r')as f:
        line=f.readline()
        while line:
            part=line.strip().split('\t')
            fwdblastbest_dic.update({str(part[0]):str(part[1])})
            fwdblastbest_dic_rev.append(str(part[1]))
            fwdblastbest_dic_both.append(str(part[0])+'\t'+str(part[1]))
            line=f.readline()
            
    #creating lists of fwdblast hits, fwdblast_besthits, and dictionary of org1 self-blast file
    rbh=[]
    rbhboth=[]
    blast=[]
    blastboth=[]
    selfblast=[]
    selfblast_dic={}
    with open(fwdblast_besthits, "r" ) as f3:
        line3=f3.readline()
        while line3:
            partblast=line3.strip().split('\t')
            blast.append(partblast[0])
            blastboth.append(str(partblast[0])+'\t'+str(partblast[1]))
            line3=f3.readline()
    with open(fwdblast_out,'r') as f4:
        fwdlist=f4.readlines()
    org1_fwd=[]
    org2_fwd=[]
    with open(fwdblast_out,'r') as f5:
        line5=f5.readline()
        while line5:
            sep=str(line5).strip().split('\t')
            org1_fwd.append(sep[0])
            org2_fwd.append(sep[1])
            line5=f5.readline()
            
        
    tdic=datetime.now()
    print('self blast dictionary making started:'+'\t'+str(tdic))
    with open(self_blast,'r') as f4:
        line4=f4.readline()
        while line4:
            part=line4.strip().split('\t')
            gene_pair=str(part[0])+'\t'+str(part[1])
            if gene_pair in selfblast_dic:
            	pass
            else:
                selfblast_dic.update({gene_pair:str(part[2])})
            line4=f4.readline()
    tdicend=datetime.now()
    print('self blast dictionary made:'+'\t'+str(tdicend))
    
    #finding tandem duplications and multiplications in org1 with respect to org2
    if not os.path.isfile( multiplications):
        multiplications_org1=[]
        if continuity=='normal':
            print('potential tandems being recovered in normal mode!')
            for every in contig_names_1:
                for each in gene_pos_per_chr_1[every]:
                    ind=gene_pos_per_chr_1[every].index(each)
                    if ind!=(len(gene_pos_per_chr_1[every])-1):
                        for eachb in blastboth:
                            parts=str(eachb).strip().split('\t')
                            if each==parts[0]:
                                each_next=gene_pos_per_chr_1[every][ind+1]
                                for everyb in blastboth:
                                    partsb=str(everyb).strip().split('\t')
                                    if each_next==partsb[0]:
                                        if parts[1]==partsb[1]:
                                            multiplications_org1.append(str(each)+'\t'+str(parts[1]))
                                            multiplications_org1.append(str(each_next)+'\t'+str(partsb[1]))
        elif continuity=='fragmented':
            print('potential tandems being recovered in fragmented mode!')
            for each in genes_org1:
                ind=genes_org1.index(each)
                if ind!=(len(genes_org1)-1):
                    for eachb in blastboth:
                        parts=str(eachb).strip().split('\t')
                        if each==parts[0]:
                            each_next=genes_org1[ind+1]
                            for everyb in blastboth:
                                partsb=str(everyb).strip().split('\t')
                                if each_next==partsb[0]:
                                    if parts[1]==partsb[1]:
                                        multiplications_org1.append(str(each)+'\t'+str(parts[1]))
                                        multiplications_org1.append(str(each_next)+'\t'+str(partsb[1]))
                                    
                                    
        with open(multiplications,'w')as out:
            for eachm in multiplications_org1:
                out.write(str(eachm)+'\n')
                
    t8=datetime.now()
    print('temporary file of tandem multiplications with anchors is ready:'+'\t'+str(t8))
    
    
    """                                
    creating a dictionary of lists with the anchor genes from org2 and tandem duplicated genes 
    (multiplied) genes from org1
    """
    multiplied_genes = {}
    with open( multiplications, "r" ) as f:
        line = f.readline()
        while line:
            parts = line.strip().split('\t')
            try:
                multiplied_genes[ parts[1] ].append( parts[0] )
            except KeyError:
                multiplied_genes.update( { parts[1]: [ parts[0] ] } )
            line = f.readline()
    gene_anchors_org2 = multiplied_genes.keys()
    
    """
    calculating the similarity among each group of duplicated genes using self blast results and writing the
    set of confident tandem multiplications into an output file
    """
    all_mean_score={}
    similar_tandems=[]
    for genes in gene_anchors_org2:
        mean_score={}
        dups_list=[]
        for each in multiplied_genes[genes]:
            score=[]
            ind=multiplied_genes[genes].index(each)
            target=multiplied_genes[genes].pop(ind)
            for every in multiplied_genes[genes]:
                bait=str(target)+'\t'+str(every)
                if bait in selfblast_dic:
                    score.append(float(selfblast_dic[bait]))
            if (len(score)!=0 and len(score)==len(multiplied_genes[genes])):
                avg=numpy.mean(score)
                mean_score.update({(str(each)):avg})
                all_mean_score.update({(str(each)):avg})
            multiplied_genes[genes].insert(ind,target)
        gene_dups=mean_score.keys()
        for dups in gene_dups:
            if (float(mean_score[dups])>float(threshold)):
                dups_list.append(dups)
        if len(dups_list)>1:
            gene_contig={}
            indice=[]
            abs_diff=[]
            counter=0
            for tandem_array in dups_list:
                for contigs in contig_names_1:
                    for one in gene_pos_per_chr_1[contigs]:
                        if str(tandem_array)==str(one):
                            ind=gene_pos_per_chr_1[contigs].index(one)
                            gene_contig.update({ one:contigs })
                            indice.append(ind)
            con=list((gene_contig.values()))
            if(con.count(con[0]) == len(con)):
                res='true'
            else:
                res='false'
            differences=numpy.diff(indice)
            for eachd in differences:
                abs_diff.append(abs(eachd))
            if ((abs_diff.count(1)==len(abs_diff)) and res=='true'):
                similar_tandems.append(str(dups_list).replace('[','').replace(']','').replace("'",'')+'\t'+str(genes))
            else:
                pass
                
    with open(tandem_multiplications,'w') as out1:
        for each in similar_tandems:
            out1.write(str(each)+'\n') 
            
    tsim=datetime.now()
    print('Highly similar tandem multiplications file is ready!!!:'+'\t'+str(tsim))  
                     
    tandems_syntelogs=[]
    tandems_potential=[]
    confident_tandems=[]
    mixed_tandems=[]
    for each in similar_tandems:
        
        l=1
        r=1
        l2=1
        r2=1
        m=1
        k=1
        org1_left=[]
        org1_right=[]
        org2_left=[]
        org2_right=[]
        part=str(each).strip().split('\t')
        parts=str(part[0]).strip().split(',')
        start=str(parts[0])
        end=str(parts[-1]).strip()
        anchor=str(part[1])
        
        #finding the left flanking region in org1
        for contig in contig_names_1:
            for gene in gene_pos_per_chr_1[contig]:
                if start==str(gene):
                    ind=gene_pos_per_chr_1[contig].index(gene)
                    posl=ind+1
                    length=len(gene_pos_per_chr_1[contig])
                    p=posl-1
                    
                    if ind==0:
                        break
                    elif ind==(length-1):
                        if length>int(flank_org1):
                            looper=int(flank_org1)
                        else:
                            looper=length
                        while l<looper:
                            org1_left.append(str(gene_pos_per_chr_1[ contig ][ind-l]))
                            l+=1
                    elif ind!=0 and ind!=(length-1):
                        if p>int(flank_org1):
                            looper=int(flank_org1)
                        elif p<=int(flank_org1):
                              looper=p
                        while l<looper:
                            org1_left.append(str(gene_pos_per_chr_1[ contig ][ind-l]))
                            l+=1
                            
        #finding the right flanking region in org1                   
        for contig in contig_names_1:
            for gene in gene_pos_per_chr_1[contig]:
                if end==str(gene):
                    inde=gene_pos_per_chr_1[contig].index(gene)
                    posr=inde+1
                    lengthe=len(gene_pos_per_chr_1[contig])
                    q=lengthe-posr
                    if inde==(lengthe-1):
                        break
                    elif inde==0:
                        if lengthe>int(flank_org1):
                            loopere=int(flank_org1)
                        else:
                            loopere=lengthe
                        while r<loopere:
                            org1_right.append(str(gene_pos_per_chr_1[ contig ][inde+r]))
                            r+=1
                    elif inde!=0 and inde!=(lengthe-1):
                        if q>int(flank_org1):
                            loopere=int(flank_org1)
                        elif q<=int(flank_org1):
                              loopere=q
                        while r<loopere:
                            org1_right.append(str(gene_pos_per_chr_1[ contig ][inde+r]))
                            r+=1
                             
        
        #finding the left and right flanking regions in org2
        for contig in contig_names_2:
            for gene in gene_pos_per_chr_2[contig]:
                if anchor==str(gene):
                    inda=gene_pos_per_chr_2[contig].index(gene)
                    lengtha=len(gene_pos_per_chr_2[contig])
                    posa=inda+1
                    p=posa-1
                    q=lengtha-posa
                    if inda==0:
                        if lengtha>int(flank_org2):
                            loopera=int(flank_org2)
                        else:
                            loopera=lengtha
                        while r2<loopera:
                            org2_right.append(str(gene_pos_per_chr_2[ contig ][inda+r2]))
                            r2+=1
                    elif inda==(lengtha-1):
                        if lengtha>int(flank_org2):
                            loopera=int(flank_org2)
                        else:
                            loopera=lengtha
                        while l2<loopera:
                            org2_left.append(str(gene_pos_per_chr_2[ contig ][inda-l2]))
                            l2+=1
                    elif inda!=0 and inda!=(lengtha-1):
                        if p>int(flank_org2):
                            looperal=int(flank_org2)
                        elif p<=int(flank_org2):
                              looperal=p
                        while l2<looperal:
                            org2_left.append(str(gene_pos_per_chr_2[ contig ][inda-l2]))
                            l2+=1
                        if q>int(flank_org2):
                            looperar=int(flank_org2)
                        elif q<=int(flank_org2):
                              looperar=q
                        while r2<looperar:
                            org2_right.append(str(gene_pos_per_chr_2[ contig ][inda+r2]))
                            r2+=1
                            
        org1_flank=org1_left+org1_right
        org2_flank=org2_left+org2_right
        best=0
        for flank in org1_flank:
            if flank in fwdblastbest_dic:
                if fwdblastbest_dic[flank] in org2_flank:
                    best+=1
        ratio=float(best)/float(2*int(flank_org1))
        
        if float(ratio)>=float(cutoff):
            tandems_syntelogs.append(each)
        else:
            tandems_potential.append(each)           
                   
    with open(non_syntelog_tandems,'w') as out:
        for each in tandems_potential:
            out.write(str(each)+'\n')
                          
    tpot=datetime.now()
    print('Non-syntelog tandem multiplications file is ready!:'+str(tpot))
                          
    with open(syntelog_tandem_file,'w') as out1:
        for each in tandems_syntelogs:
            out1.write(str(each)+'\n')
                
    t9=datetime.now()
    print('Syntelog tandem multiplications file is ready!:'+'\t'+str(t9))
    #plotting histogram of mutual gene similarities in the tandem arrays
    histo=[]
    for each in tandems_syntelogs:
        cutting=str(each).strip().split('\t')
        splitting=str(cutting[0]).strip().split(',')
        for every in splitting:
            tester=str(every).strip()
            histo.append(all_mean_score[tester])
    for each in tandems_potential:
        cutting=str(each).strip().split('\t')
        splitting=str(cutting[0]).strip().split(',')
        for every in splitting:
            tester=str(every).strip()
            histo.append(all_mean_score[tester])
    sns.set()
    _=pyplot.hist(histo)
    _=pyplot.xlabel('mutual similarity scores')
    _=pyplot.ylabel('number of genes showing the mutual similarity')
    pyplot.savefig(histogram, dpi=300)
                
    
                
    tfinish=datetime.now()
    print("All done!:"+'\t'+str(tfinish))            
        
def main( arguments ):
    """! @brief extract transcript and pepetides from FASTA and GFF3 """
    gff3_input_file_1 = arguments[ arguments.index('--gff')+1 ]
    gff3_input_file_2 = arguments[ arguments.index('--gff')+2 ]
    fasta_input_file_1 = arguments[ arguments.index('--fasta')+1 ]
    fasta_input_file_2 = arguments[ arguments.index('--fasta')+2 ]
    output_dir = arguments[ arguments.index('--out')+1 ]
    tmp_dir=arguments[arguments.index('--tmp')+1]
    
    if '--name' in arguments:
        name1 = arguments[ arguments.index('--name')+1 ]
        name2 = arguments[arguments.index('--name')+2]
    else:
        name1 = "org1"
        name2 = "org2"
    if output_dir[ -1 ] != "/":
        output_dir += "/"
    if not os.path.exists( output_dir ):
        os.makedirs( output_dir )
    if tmp_dir[ -1 ] != "/":
        tmp_dir += "/"
    if not os.path.exists( tmp_dir ):
        os.makedirs( tmp_dir )
    if '--cpu' in arguments:
        cpu = arguments[arguments.index('--cpu')+1]
    else:
        cpu = 10
    if '--similarity' in arguments:
        threshold = arguments[arguments.index('--similarity')+1]
    else:
        threshold=80
    if '--flankq' in arguments:
        flank_org1=arguments[arguments.index('--flankq')+1]
    else:
        flank_org1=6
    if '--flankr' in arguments:
        flank_org2=arguments[arguments.index('--flankr')+1]
    else:
        flank_org2=6
    if '--cutoff' in arguments:
        cutoff=arguments[arguments.index('--cutoff')+1]
    else:
        cutoff=0.5
    if '--continuity' in arguments:
        continuity=arguments[arguments.index('--continuity')+1]
    else:
        continuity='normal'
        
    pep_output_file_1 = tmp_dir + name1 + ".pep.fasta"
    CDS_output_file1 = tmp_dir + name1 + ".cds.fasta"
    pep_output_file_2 = tmp_dir + name2 + ".pep.fasta"
    CDS_output_file2 = tmp_dir + name2 + ".cds.fasta"
    
    gff=[gff3_input_file_1, gff3_input_file_2]
    wgs=[fasta_input_file_1, fasta_input_file_2]
    cds=[CDS_output_file1, CDS_output_file2]
    pep=[pep_output_file_1, pep_output_file_2]
    strict_start = False #option not implemented yet
    strict_end = False #option not implemented yet
    
    #calling functions with GFF3 and FASTA file inputs of organism 1
    seqs_1 = load_sequences( fasta_input_file_1)
    transcript_information_1 = load_transcript_information_from_gff3( gff3_input_file_1 )
    construct_CDS_file( transcript_information_1, CDS_output_file1, seqs_1 )
    transeq( CDS_output_file1, pep_output_file_1, strict_start, strict_end )
    
    #calling functions with GFF3 and FASTA file inputs of organism 2
    seqs_2 = load_sequences( fasta_input_file_2)
    transcript_information_2 = load_transcript_information_from_gff3( gff3_input_file_2 )
    construct_CDS_file( transcript_information_2, CDS_output_file2, seqs_2 )
    transeq( CDS_output_file2, pep_output_file_2, strict_start, strict_end )
    
    
    #calling the function to perform reciprocal blast between organism 1 and 2
    gene_multiplications (output_dir, tmp_dir, name1, name2, pep_output_file_1, pep_output_file_2, gff3_input_file_1, gff3_input_file_2, cpu, threshold, flank_org1, flank_org2, cutoff, continuity)
        
if '--gff' in sys.argv and '--fasta' in sys.argv and '--out' in sys.argv and '--tmp' in sys.argv:
    main( sys.argv )
else:
    sys.exit( __usage__ )

