#!/usr/bin/python
# -*- coding: utf-8 -*-
"""

fisher.py


"""
import sys
import os
import re
import pysam
import scipy.special
from scipy.stats import fisher_exact as fisher
import numpy as np
import argparse
import logging
import math
import subprocess

#
# Globals
#
arg = None
target = None
remove_chr = None
filter_quals = None

#
# Class definitions
#

############################################################
class AutoVivification(dict):
    """Implementation of perl's autovivification feature."""
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)

        except KeyError:
            value = self[item] = type(self)()

        return value

#
# Subroutines
#

# data_pair IDs
POS_CHR = 0
POS_COORD = 1
POS_REF = 2
POS_DATA1 = 3
POS_DATA2 = 4
# POS_FISHER_SNV = 5
# POS_FISHER_INS = 6
# POS_FISHER_DEL = 7
POS_COUNT = 5

def math_log_fisher_pvalue(fisher_pvalue):

    val = float(0.0)
    if fisher_pvalue < 10**(-60):
        val = float(60.0)
    elif fisher_pvalue  > 1.0 - 10**(-10) :
        val = float(0.0)
    else:
        val = -math.log( fisher_pvalue, 10 )
                
    return val


############################################################
def Pileup_out( mpileup, min_depth, min_variant_read, compare ):

    #
    # mpileup format
    #
    # chr1 272 T 24  ,.$.....,,.,.,...,,,.,..^+. <<<+;<<<<<<<<<<<=<;<;7<&
    #
    # 0 chromosome,
    # 1 1-based coordinate,
    # 2 reference base,
    # 3 the number of reads covering the site (1)
    # 4 read bases (1)
    # 5 base qualities (1)
    # 6 the number of reads covering the site (2)
    # 7 read bases (2)
    # 8 base qualities (2)
    #
    global target
    global remove_chr
    global filter_quals


    #
    # Prepare mpileup data
    #
    mp_list = str( mpileup.translate( None, '\n' ) ).split( '\t' )
    mp_list_len = len( mp_list )
    ref_base_U = mp_list[ 2 ].upper()
    coordinate = mp_list[ 0:3 ]

    #
    # skip if depth is 0
    #
    if mp_list[ 3 ] == '0' or ( mp_list_len > 7 and mp_list[ 7 ] == '0' ):
    # if int(mp_list[ 3 ]) < min_depth or ( mp_list_len > 6 and int(mp_list[ 7 ]) < min_depth ):
        return None

    ###Fix:0.2.1: not used
    # ref_base_plus  = mp_list[ 4 ].count('.')
    # ref_base_minus = mp_list[ 4 ].count(',')

    ref_base_count = mp_list[ 4 ].count('.') + mp_list[ 4 ].count(',')
    ins_base_count = mp_list[ 4 ].count('+')
    del_base_count = mp_list[ 4 ].count('-')
    ###Fix:0.2.1: skip_base_count (because of deletion)
    skip_base_count = mp_list[ 4 ].count('*') + mp_list[ 4 ].count('#')
    if int(mp_list[ 3 ])  < min_depth or \
       (int(mp_list[ 3 ]) - ref_base_count - skip_base_count + ins_base_count + del_base_count) < min_variant_read:
        return None

    if ref_base_U not in 'ACGTN': return None
    #
    # data_pair IDs
    # POS_CHR = 0
    # POS_COORD = 1
    # POS_REF = 2
    # POS_DATA1 = 3
    # POS_DATA2 = 4
    # # POS_FISHER_SNV = 5
    # # POS_FISHER_INS = 6
    # # POS_FISHER_DEL = 7
    # POS_COUNT = 5
    #
    data_pair = [ mp_list[ 0 ],
                  int( mp_list[ 1 ] ),
                  mp_list[ 2 ],
                  { 'mis_base': ref_base_U, 'mis_rate': 0, 'proper_read_depth': 0, 'proper_read_depth_plus': 0, 'proper_read_depth_minus': 0, 'proper_read_depth_indel': 0, 'proper_read_depth_indel_plus': 0, 'proper_read_depth_indel_minus': 0,'indel': AutoVivification(), 'snv': AutoVivification() },
                  { 'mis_base': ref_base_U, 'mis_rate': 0, 'proper_read_depth': 0, 'proper_read_depth_plus': 0, 'proper_read_depth_minus': 0, 'proper_read_depth_indel': 0, 'proper_read_depth_indel_plus': 0, 'proper_read_depth_indel_minus': 0,'indel': AutoVivification(), 'snv': AutoVivification() },
                  0 ]


    #
    # Loop for 2 bam file case
    #
    if compare:
        data_pair[ POS_COUNT ] = 2
        ###Fix:0.2.1: for qname including
        input_list = [ ( POS_DATA1, mp_list[ 3 ], mp_list[ 4 ], mp_list[ 5 ], mp_list[ 6 ] ),
                       ( POS_DATA2, mp_list[ 7 ], mp_list[ 8 ], mp_list[ 9 ], mp_list[ 10 ] ) ]
    else:
        data_pair[ POS_COUNT ] = 1
        input_list = [ ( POS_DATA1, mp_list[ 3 ], mp_list[ 4 ], mp_list[ 5 ], mp_list[ 6 ] ) ]

    #
    # position id,
    # mpileup output 4th row(number of read covering the site),
    # 5th row(read bases),
    # 6th row(base quality)
    #
    for data_id, depth, read_bases, qual_list, qnames in input_list:
                
        #
        # Remove '^.' and '$'
        #
        read_bases = remove_chr.sub( '', read_bases )
        read_bases = read_bases.translate( None, '$' ) 
        
        qname_list = qnames.split(',')
                
        indel = AutoVivification()

        #
        # Look for deletion/insertion and save info in 'indel' dictionary
        #
        #   ([\+\-])[0-9]+[ACGTNacgtn]+
        #
        # m.group(1): + or - (deletion/insertion)
        # m.group(2): number of deletion/insertion
        # m.group(3): nucleotides
        #
        deleted = 0
        iter = target.finditer( read_bases )
        for m in iter:
            site = m.start()
            indeltype = m.group( 1 )
            num = m.group( 2 )
            bases = m.group( 3 )[ 0:int( num ) ]
            if bases.islower():
                strand = ( '-', '+' )
            else:
                strand = ( '+', '-' )
            ###Fix:0.2.1: In previous method, ovelapped indel was double-counted.
            ###           For single count, here only qnames are added.
            k_base_before = site - deleted - 1
            baseqname = qname_list[ k_base_before ]
    
            key = '\t'.join( coordinate + [ bases.upper() ] )
            if indeltype in indel and key in indel[ indeltype ]:
                # indel[ indeltype ][ key ][ strand[ 0 ] ] += 1
                indel[ indeltype ][ key ][ strand[ 0 ] ].append(baseqname)
            else:
                # indel[ indeltype ][ key ][ strand[ 0 ] ] = 1
                # indel[ indeltype ][ key ][ strand[ 1 ] ] = 0
                indel[ indeltype ][ key ][ strand[ 0 ] ] = [baseqname]
                indel[ indeltype ][ key ][ strand[ 1 ] ] = []
    
            read_bases = read_bases[ 0:site - deleted ] + read_bases[ site + int(num) + len( num ) + 1 - deleted: ]
            deleted += 1 + len( num ) + int( num )
                
    

        #
        # Error check
        #
        if len( read_bases ) != len( qual_list ) or len( read_bases ) != len(qname_list):
            print("mpileup data is not good: {0}, {1}, {2}".format( mpileup, read_bases, qnames ))
            return None
            logging.error( "mpileup data is not good: {0}, {1}".format( mpileup, read_bases ) )
            return None
        # Note:
        # Even in the skip ('* / #') base, base quality and qname are described: same as the last base.

        #
        # Count mismatch
        #
        # mis_base_U = None
        
        bases_plus  = ('A','C','G','T','N','*')
        bases_minus = ('a','c','g','t','n','#')
        
         # if int( depth ) >= min_depth:

        read_bases = read_bases.replace( '.', ref_base_U )
        read_bases = read_bases.replace( ',', ref_base_U.lower() )

        base_num = {
            "total_A": 0,
            "total_C": 0,
            "total_G": 0,
            "total_T": 0,
            "A": 0,
            "C": 0,
            "G": 0,
            "T": 0,
            "a": 0,
            "c": 0,
            "g": 0,
            "t": 0,
        }

        #
        # Set data
        #
        data_pair[ data_id ][ 'bases' ] = read_bases # not used, but useful for debug
        data_pair[ data_id ][ 'depth' ] = int( depth )

        #
        # Count number
        #
        
        ###Fix:0.2.1: In previous method, the depth and variant read of indel is calculated using both overlapped reads.
        # In this version, variant is called when at least one read (of overlap) has variant base (or indel)
        read_bases_ar = np.array(list(read_bases))
        qual_list_ar  = np.array(list(qual_list))
        qname_list_ar = np.array(list(qname_list))
        
        isplus  = np.isin(read_bases_ar, bases_plus)
        isminus = np.isin(read_bases_ar, bases_minus)
        
        isqual  = ~np.isin(qual_list_ar, filter_quals)
        
        #np.isinより少し速い
        isA     = (read_bases_ar == 'A') | (read_bases_ar == 'a')
        isC     = (read_bases_ar == 'C') | (read_bases_ar == 'c')
        isG     = (read_bases_ar == 'G') | (read_bases_ar == 'g')
        isT     = (read_bases_ar == 'T') | (read_bases_ar == 't')
        isN     = (read_bases_ar == 'N') | (read_bases_ar == 'n')
        isskip  = (read_bases_ar == '*') | (read_bases_ar == '#')
        
        qname_list_plus     = qname_list_ar[ isplus  ]
        qname_list_minus    = qname_list_ar[ isminus ]
        
        qname_list_qual       = qname_list_ar[ isqual & ~isN & ~isskip ]
        qname_list_qual_plus  = qname_list_ar[ isqual & ~isN & ~isskip & isplus  ]
        qname_list_qual_minus = qname_list_ar[ isqual & ~isN & ~isskip & isminus ]
        
        # depth_indel is not necessarily equal to depth_indel_plus + depth_indel_minus (because of overlap)
        #　np.uniqueを使うよりもsetを使ったほうが速い
        data_pair[ data_id ][ 'proper_read_depth_indel' ]       = len(set(qname_list))
        data_pair[ data_id ][ 'proper_read_depth_indel_plus' ]  = len(set(qname_list_plus))
        data_pair[ data_id ][ 'proper_read_depth_indel_minus' ] = len(set(qname_list_minus))
        
        data_pair[ data_id ][ 'proper_read_depth' ]       = len(set(qname_list_qual))
        data_pair[ data_id ][ 'proper_read_depth_plus' ]  = len(set(qname_list_qual_plus))
        data_pair[ data_id ][ 'proper_read_depth_minus' ] = len(set(qname_list_qual_minus))
        
        # total_A is not necessarily (A + a)
        # It is guaranteed that base_num is less than 'proper_read_depth'
        for nuc, isB in zip(('A','C','G','T'), (isA, isC, isG, isT)):
            base_num[ nuc            ] = len(set(qname_list_ar[ isqual & isB & isplus  ]))
            base_num[ nuc.lower()    ] = len(set(qname_list_ar[ isqual & isB & isminus ]))
            base_num[ 'total_' + nuc ] = len(set(qname_list_ar[ isqual & isB ]))
        
        # for nuc, qual in zip( read_bases, qual_list ):
            # if nuc in 'ATGCNacgtn':
                # data_pair[ data_id ][ 'proper_read_depth_indel' ] += 1 
            # if nuc in 'ATGCN':
                # data_pair[ data_id ][ 'proper_read_depth_indel_plus' ] += 1 
            # if nuc in 'acgtn':
                # data_pair[ data_id ][ 'proper_read_depth_indel_minus' ] += 1 
            # if nuc in 'ATGCNacgtn' and not ( qual in filter_quals) :
                # base_num[ nuc ] += 1
                # base_num[ 'total_' + nuc.upper() ] += 1
            # if nuc in 'ATGCatgc' and not ( qual in filter_quals):
                # data_pair[ data_id ][ 'proper_read_depth' ] += 1 
            # if nuc in 'ATGC' and not ( qual in filter_quals):
                # data_pair[ data_id ][ 'proper_read_depth_plus' ] += 1 
            # if nuc in 'atgc' and not ( qual in filter_quals):
                # data_pair[ data_id ][ 'proper_read_depth_minus' ] += 1 

        #
        # InsDel
        #
        for indeltype in ( '+', '-' ):
            if indeltype in indel:
                for key in indel[ indeltype ].keys():
                    bases = key.split( '\t' )[ 3 ]
                    # both is not necessarily '+' + '-'
                    indel_number_plus  = len(set(indel[ indeltype ][ key ][ '+' ]))
                    indel_number_minus = len(set(indel[ indeltype ][ key ][ '-' ]))
                    indel_number       = len(set( indel[ indeltype ][ key ][ '-' ] + indel[ indeltype ][ key ][ '+' ] ))
                    
                    # if data_id == POS_DATA1 and indel_number < min_variant_read: continue
                    # if ( data_id == POS_DATA2 and 
                         # (indeltype not in data_pair[ POS_DATA1 ]['indel'] or
                          # bases     not in data_pair[ POS_DATA1 ]['indel'][indeltype]
                         # )
                       # ): continue
                    
                    data_pair[ data_id ][ 'indel' ][ indeltype ][ bases ][ '+' ]    = indel_number_plus
                    data_pair[ data_id ][ 'indel' ][ indeltype ][ bases ][ '-' ]    = indel_number_minus
                    data_pair[ data_id ][ 'indel' ][ indeltype ][ bases ][ 'both' ] = indel_number
                        
                                                                                   

                    ###Fix:0.2.1: not calculate when pair analysis
                    if data_pair[ POS_COUNT ] == 1:
                        data_pair[ data_id ][ 'indel' ][ indeltype ][ bases ][ '0.1' ] = \
                            scipy.special.btdtri( indel_number + 1, float( data_pair[ data_id ][ 'proper_read_depth_indel' ] ) - indel_number + 1, 0.1 )
                        data_pair[ data_id ][ 'indel' ][ indeltype ][ bases ][ 'mid' ] = \
                            ( indel_number + 1 ) / ( float( data_pair[ data_id ][ 'proper_read_depth_indel' ] ) + 2 )
                        data_pair[ data_id ][ 'indel' ][ indeltype ][ bases ][ '0.9' ] = \
                            scipy.special.btdtri( indel_number + 1, int( data_pair[ data_id ][ 'proper_read_depth_indel' ] ) - indel_number + 1, 0.9 )
                        
                    data_pair[ data_id ][ 'indel' ][ indeltype ][ bases ][ 's_ratio' ] = \
                        float( indel_number_plus ) / (indel_number_plus + indel_number_minus) # not equal to indel_number
        #
        # SNV
        # skip if reference is 'N'
        #
        if ref_base_U != 'N':
            # ref_num = base_num[ 'total_' + ref_base_U ]
            # mis_num = 0
            
            data_pair[ data_id ].update(base_num)
            
            # for nuc in ( 'A', 'C', 'G', 'T' ):
                # data_pair[ data_id ][ nuc ] = base_num[ nuc ]
                # tmp = nuc.lower()
                # data_pair[ data_id ][ tmp ] = base_num[ tmp ]
                # tmp = 'total_' + nuc
                # data_pair[ data_id ][ tmp ] = base_num[ tmp ]

                # if nuc != ref_base_U:
                    # if base_num[ tmp ] > mis_num:
                        # mis_num = base_num[ tmp ]
                        # mis_base_U = nuc

            if data_id == POS_DATA1: 
                mis_bases = [nuc for nuc in ( 'A', 'C', 'G', 'T' )
                                if base_num['total_' + nuc] >= min_variant_read and nuc != ref_base_U]
            
            # Note that mutation changed to reference from SNP cannot be detected.
            elif data_id == POS_DATA2 and len(data_pair[ POS_DATA1 ][ 'snv' ]) > 0:
                # mis_num = base_num[ 'total_' + data_pair[ POS_DATA1 ][ 'mis_base' ] ]
                mis_bases = data_pair[ POS_DATA1 ][ 'snv' ].keys() # In python2, keys() returns list

            #
            # Calculate ratio
            #
            ###Fix:0.2.1: fix strandratio
            ###           change structure similar to indel to handle multiple mutations at same point.
            proper_depth = data_pair[ data_id ][ 'proper_read_depth' ]
            for k, mis_base in enumerate(mis_bases):
                mis_num      = base_num['total_' + mis_base]
                ref_num      = proper_depth - mis_num
                mis_rate     = mis_num / float( proper_depth ) if proper_depth != 0 else 0
                sum_mis_num = base_num[ mis_base ] + base_num[ mis_base.lower() ] 
                if sum_mis_num > 0:
                    s_ratio  = base_num[ mis_base ] / float( sum_mis_num )
                else:
                    s_ratio  = None
                
                data_pair[ data_id ][ 'snv' ][ mis_base ][ 'mis_rate' ] = mis_rate
                data_pair[ data_id ][ 'snv' ][ mis_base ][ 'mis_num' ]  = mis_num
                data_pair[ data_id ][ 'snv' ][ mis_base ][ 's_ratio' ]  = s_ratio
                
            
            # data_pair[ data_id ][ 'mis_rate' ] = mis_num / float( data_pair[ data_id ][ 'proper_read_depth' ] )
            # data_pair[ data_id ][ 'mis_base' ] = mis_base_U
            # if mis_base_U and ( base_num[ mis_base_U ] + base_num[ mis_base_U.lower() ] ) > 0:
                # data_pair[ data_id ][ 's_ratio' ]  = float( base_num[ mis_base_U ] ) / ( base_num[ mis_base_U ] + base_num[ mis_base_U.lower() ] )
            # else:
            #    data_pair[ data_id ][ 's_ratio' ]  = float(0)

            #
            # Beta distribution for SNV
            #
                ###Fix:0.2.1: not calculate when pair analysis
                if data_pair[ POS_COUNT ] == 1:
                    data_pair[ data_id ][ 'snv' ][ mis_base ][ '0.1' ] = scipy.special.btdtri( mis_num + 1, ref_num + 1, 0.1 )
                    data_pair[ data_id ][ 'snv' ][ mis_base ][ 'mid' ] = ( mis_num + 1 ) / float( ref_num + mis_num + 2 )
                    data_pair[ data_id ][ 'snv' ][ mis_base ][ '0.9' ] = scipy.special.btdtri( mis_num + 1, ref_num + 1, 0.9 )
                
    #
    # Fisher
    #
    # SNV
    #
    # Note that fisher p-value can be small if alt_tumor is significantly "less" than alt_normal
    if ( data_pair[ POS_COUNT ] == 2 and
         ref_base_U != 'N' and
         len(data_pair[ POS_DATA1 ][ 'snv' ]) > 0
       ):
        for mis_base in data_pair[ POS_DATA1 ][ 'snv' ].keys():
            alt_tumor  = data_pair[ POS_DATA1 ][ 'snv' ][ mis_base ][ 'mis_num' ]
            alt_normal = data_pair[ POS_DATA2 ][ 'snv' ][ mis_base ][ 'mis_num' ]
            
            ref_tumor  = data_pair[ POS_DATA1 ][ 'proper_read_depth' ] - alt_tumor
            ref_normal = data_pair[ POS_DATA2 ][ 'proper_read_depth' ] - alt_normal
            
            _, fisher_pvalue = fisher(
                        ( ( ref_tumor, ref_normal ), ( alt_tumor, alt_normal ) ),
                        alternative='two-sided'
                        )

            data_pair[ POS_DATA1 ][ 'snv' ][ mis_base ][ 'fisher' ] = math_log_fisher_pvalue(fisher_pvalue)

    #
    # INDEL
    #
    if ( data_pair[ POS_COUNT ] == 2 and len( data_pair[ POS_DATA1 ][ 'indel' ] ) > 0 ):
        fisher_pvalue = None
        for indeltype in data_pair[ POS_DATA1 ][ 'indel' ]:
            for bases in data_pair[ POS_DATA1 ][ 'indel' ][ indeltype ].keys():
              
                if not isinstance( data_pair[ POS_DATA2 ][ 'indel' ][ indeltype ][ bases ][ 'both' ], int ):
                    data_pair[ POS_DATA2 ][ 'indel' ][ indeltype ][ bases ][ 'both' ] = 0
                    data_pair[ POS_DATA2 ][ 'indel' ][ indeltype ][ bases ][ '+' ] = 0
                    data_pair[ POS_DATA2 ][ 'indel' ][ indeltype ][ bases ][ '-' ] = 0
                    data_pair[ POS_DATA2 ][ 'indel' ][ indeltype ][ bases ][ 's_ratio' ] = None

                #Maybe this could never happen.
                # if (data_pair[ POS_DATA2 ][ 'proper_read_depth_indel' ] >= data_pair[ POS_DATA2 ][ 'indel' ][ indeltype ][ bases ][ 'both' ] and
                    # data_pair[ POS_DATA1 ][ 'proper_read_depth_indel' ] >= data_pair[ POS_DATA1 ][ 'indel' ][ indeltype ][ bases ][ 'both' ]
                    # ):
                    
                alt_tumor  = data_pair[ POS_DATA1 ][ 'indel' ][ indeltype ][ bases ][ 'both' ]
                alt_normal = data_pair[ POS_DATA2 ][ 'indel' ][ indeltype ][ bases ][ 'both' ]
                ref_tumor  = data_pair[ POS_DATA1 ][ 'proper_read_depth_indel' ] - alt_tumor
                ref_normal = data_pair[ POS_DATA2 ][ 'proper_read_depth_indel' ] - alt_normal
                
                _, fisher_pvalue = fisher(
                            ( ( ref_tumor, ref_normal ), ( alt_tumor, alt_normal ) ),
                            alternative='two-sided'
                            )
                            
    
                if fisher_pvalue != None:
                    data_pair[ POS_DATA1 ][ 'indel' ][ indeltype ][ bases ][ 'fisher' ] = math_log_fisher_pvalue(fisher_pvalue)
                else:
                    data_pair[ POS_DATA1 ][ 'indel' ][ indeltype ][ bases ][ 'fisher' ] = None
                    # if indeltype == '+':
                        # data_id = POS_FISHER_INS
                    # elif indeltype == '-':
                        # data_id = POS_FISHER_DEL

                    # if data_pair[ data_id ] == 'N:1.0':
                        # data_pair[ data_id ] = bases + ':' + str( math_log_fisher_pvalue(fisher_pvalue) )
                    # else:
                        # data_pair[ data_id ] += ',' + bases + ':' + str( math_log_fisher_pvalue(fisher_pvalue) )


    return data_pair


############################################################
def print_data( data, w, min_depth, mismatch_rate_disease, mismatch_rate_normal, posterior_10_quantile, fisher_threshold, min_variant_read ):
    str_indel_dict = AutoVivification()
    f_print_indel = False

    if data[ POS_COUNT ] == 1:
        #
        # barcode SNV output
        #
        data_tumor  = data[ POS_DATA1 ]
        #
        # Fisher
        #
        for mis_base in sorted(data_tumor[ 'snv' ].keys()):
            data_snv_tumor  = data_tumor[ 'snv' ][ mis_base ]
            
            if not (
                data_snv_tumor[ 'mis_rate' ]            >  mismatch_rate_disease    and
                data_tumor[  'proper_read_depth' ]      >= min_depth                and
                data_snv_tumor['0.1']                   >  posterior_10_quantile    and
                data_snv_tumor[ 'mis_num' ]             >= min_variant_read
               ): continue
               
            # Genomon output for barcode
            # chr \t start \t end \t ref \t obs \tdepth \t A,C,G,T \t mis \t s_ratio \t 0.1 \t ratio \t 0.9
            # outstr = '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7},{8},{9},{10}\t{11:.3f}\t{12:.3f}\t{13:.3f}\t{14:.3f}\t{15:.3f}\n'.format(
            outstr = '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7},{8},{9},{10}\t{11},{12},{13},{14}\t{15:.3f}\t{16:.3f}\t{17:.3f}\t{18:.3f}\t{19:.3f}\n'.format(
                        data[ POS_CHR ],
                        data[ POS_COORD ],
                        data[ POS_COORD ],
                        data[ POS_REF ],
                        mis_base,
                        data_tumor[ 'proper_read_depth' ],
                        data_snv_tumor[ 'mis_num' ],
                        data_tumor[ 'proper_read_depth_plus' ],
                        data_tumor[ mis_base ],
                        data_tumor[ 'proper_read_depth_minus' ],
                        data_tumor[ mis_base.lower() ],
                        data_tumor[ 'total_A' ],
                        data_tumor[ 'total_C' ],
                        data_tumor[ 'total_G' ],
                        data_tumor[ 'total_T' ],
                        data_snv_tumor[ 'mis_rate' ],
                        data_snv_tumor[ 's_ratio' ],
                        data_snv_tumor[ '0.1' ],
                        data_snv_tumor[ 'mid' ],
                        data_snv_tumor[ '0.9' ],
                        )
                        
            w.write( outstr)

        #
        # InDel output
        #
        for data_type in ( '-', '+' ):
            if data_type not in data_tumor[ 'indel' ]: continue
            
            for indel_bases in sorted( data_tumor[ 'indel' ][ data_type ] ):
                data_indel_tumor  = data_tumor[ 'indel' ][ data_type ][ indel_bases ]
                misrate_tumor     = data_indel_tumor[ 'both' ] / float(data_tumor[ 'proper_read_depth_indel' ])
            
                if not (
                    misrate_tumor  > mismatch_rate_disease                  and
                    data_tumor[  'proper_read_depth_indel' ] >= min_depth   and
                    data_indel_tumor[ 'both' ]   >= min_variant_read        and
                    data_indel_tumor[ '0.1'  ]   >= posterior_10_quantile
                   ): continue

                outstr = '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7},{8},{9},{10}\t{11}\t{12:.3f}\t{13:.3f}\t{14:.3f}\t{15:.3f}\t{16:.3f}\n'.format(
                            data[ POS_CHR ],
                            data[ POS_COORD ] + 1 if data_type == '-' else data[ POS_COORD ],
                            data[ POS_COORD ] + len( indel_bases ) if data_type == '-' else data[ POS_COORD ],
                            indel_bases if data_type == '-' else '-',
                            indel_bases if data_type == '+' else '-',
                            data_tumor[ 'proper_read_depth_indel' ],
                            data_indel_tumor[ 'both' ],
                            data_tumor[ 'proper_read_depth_indel_plus' ],
                            data_indel_tumor[ '+' ],
                            data_tumor[ 'proper_read_depth_indel_minus' ],
                            data_indel_tumor[ '-' ],
                            '---',
                            misrate_tumor,
                            data_indel_tumor['s_ratio'],
                            data_indel_tumor[ '0.1' ],
                            data_indel_tumor[ 'mid' ],
                            data_indel_tumor[ '0.9' ],
                            )
                            
                w.write(outstr)

    elif data[ POS_COUNT ] == 2:

        fisher_threshold_log = -math.log( fisher_threshold, 10 )
        
        data_tumor  = data[ POS_DATA1 ]
        data_normal = data[ POS_DATA2 ]
        #
        # Fisher
        #
        for mis_base in sorted(data_tumor[ 'snv' ].keys()):
            data_snv_tumor  = data_tumor[  'snv' ][ mis_base ]
            data_snv_normal = data_normal[ 'snv' ][ mis_base ]
            
            if not (
                data_snv_normal['mis_rate' ]            <  mismatch_rate_normal     and
                data_snv_tumor[ 'mis_rate' ]            >  mismatch_rate_disease    and
                data_normal[ 'proper_read_depth' ]      >= min_depth                and
                data_tumor[  'proper_read_depth' ]      >= min_depth                and
                data_snv_tumor['fisher']                is not None                 and
                data_snv_tumor['fisher']                >  fisher_threshold_log     and
                data_snv_tumor[ 'mis_num'  ]            >= min_variant_read
               ): continue
               
            #
            # Genomon output for fisher by comparing nomral and tumor
            # chr \t start \t end \t ref \t Alt
            # depth_tumor variantNum_tumor depth_normal variantNum_normal bases_tumor bases_normal
            # A_C_G_T_tumor A_C_G_T_normal misRate_tumor strandRatio_tumor misRate_normal strandRatio_normal 
            # P-value(fisher)

            # data IDs
            # POS_CHR = 0
            # POS_COORD = 1
            # POS_REF = 2
            # POS_DATA1 = 3
            # POS_DATA2 = 4
            # POS_COUNT = 5
            ###
            outstr = ('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}'.format(
                        data[ POS_CHR ],
                        data[ POS_COORD ],
                        data[ POS_COORD ],
                        data[ POS_REF ],
                        mis_base,
                        data_tumor[ 'proper_read_depth' ],
                        data_snv_tumor[ 'mis_num' ],
                        data_normal[ 'proper_read_depth' ],
                        data_snv_normal[ 'mis_num' ],
                        )
                     + '\t{0},{1},{2},{3}'.format(
                        data_tumor[ 'proper_read_depth_plus' ],
                        data_tumor[ mis_base ],
                        data_tumor[ 'proper_read_depth_minus' ],
                        data_tumor[ mis_base.lower() ],
                        )
                     + '\t{0},{1},{2},{3}'.format(
                        data_normal[ 'proper_read_depth_plus' ],
                        data_normal[ mis_base ],
                        data_normal[ 'proper_read_depth_minus' ],
                        data_normal[ mis_base.lower() ],
                        )
                     + '\t{0},{1},{2},{3}'.format(
                        data_tumor[ 'total_A' ],
                        data_tumor[ 'total_C' ],
                        data_tumor[ 'total_G' ],
                        data_tumor[ 'total_T' ],
                        )
                     + '\t{0},{1},{2},{3}'.format(
                        data_normal[ 'total_A' ],
                        data_normal[ 'total_C' ],
                        data_normal[ 'total_G' ],
                        data_normal[ 'total_T' ],
                        )
                     + '\t{0:.3f}'.format(data_snv_tumor[ 'mis_rate' ])
                     + '\t{0:.3f}'.format(data_snv_tumor[ 's_ratio' ])
                     + '\t{0:.3f}'.format(data_snv_normal[ 'mis_rate' ])
                     )
            if data_snv_normal['s_ratio'] is not None:
                outstr += '\t{0:.3f}'.format(data_snv_normal[ 's_ratio' ])
            else:
                outstr += '\t---'
            outstr += '\t{0:.3f}'.format(data_snv_tumor[ 'fisher' ])
                     
            w.write( outstr +"\n")

        for data_type in ( '-', '+' ):
            if data_type not in data_tumor[ 'indel' ]: continue
            
            for indel_bases in sorted( data_tumor[ 'indel' ][ data_type ] ):
                data_indel_tumor  = data_tumor[  'indel' ][ data_type ][ indel_bases ]
                data_indel_normal = data_normal[ 'indel' ][ data_type ][ indel_bases ]
                
                if data_normal[ 'proper_read_depth_indel' ] == 0: continue
                
                misrate_tumor  = data_indel_tumor[  'both' ] / float(data_tumor[  'proper_read_depth_indel' ])
                misrate_normal = data_indel_normal[ 'both' ] / float(data_normal[ 'proper_read_depth_indel' ])
            
                if not (
                    misrate_normal < mismatch_rate_normal                   and
                    misrate_tumor  > mismatch_rate_disease                  and
                    data_normal[ 'proper_read_depth_indel' ] >= min_depth   and
                    data_tumor[  'proper_read_depth_indel' ] >= min_depth   and
                    data_indel_tumor[ 'both' ]   >= min_variant_read        and
                    data_indel_tumor[ 'fisher' ] is not None                and
                    data_indel_tumor[ 'fisher' ] >= fisher_threshold_log
                   ): continue
                
                #
                # Genomon output for fisher by comparing nomral and tumor
                # chr \t start \t end \t ref \t Alt
                # depth_tumor variantNum_tumor depth_normal variantNum_normal bases_tumor bases_normal
                # A_C_G_T_tumor A_C_G_T_normal misRate_tumor strandRatio_tumor misRate_normal strandRatio_normal 
                # P-value(fisher)
                outstr = ('{0}\t{1}\t{2}\t{3}\t{4}'.format(
                            data[ POS_CHR ],
                            data[ POS_COORD ] + 1 if data_type == '-' else data[ POS_COORD ],
                            data[ POS_COORD ] + len( indel_bases ) if data_type == '-' else data[ POS_COORD ],
                            indel_bases if data_type == '-' else '-',
                            indel_bases if data_type == '+' else '-'
                            )
                         + '\t{0}\t{1}\t{2}\t{3}'.format(
                            data_tumor[ 'proper_read_depth_indel' ],
                            data_indel_tumor[ 'both' ],
                            data_normal[ 'proper_read_depth_indel' ],
                            data_indel_normal[ 'both' ],
                            )
                         + '\t{0},{1},{2},{3}'.format(
                            data_tumor[ 'proper_read_depth_indel_plus' ],
                            data_indel_tumor[ '+' ],
                            data_tumor[ 'proper_read_depth_indel_minus' ],
                            data_indel_tumor[ '-' ],
                            )
                         + '\t{0},{1},{2},{3}'.format(
                            data_normal[ 'proper_read_depth_indel_plus' ],
                            data_indel_normal[ '+' ],
                            data_normal[ 'proper_read_depth_indel_minus' ],
                            data_indel_normal[ '-' ],
                            )
                         + '\t{0}\t{1}'.format(
                            '---',
                            '---',
                            )
                         + '\t{0:.3f}'.format(misrate_tumor)
                         + '\t{0:.3f}'.format(data_indel_tumor['s_ratio'])
                         + '\t{0:.3f}'.format(misrate_normal)
                         )
                if data_indel_normal['s_ratio'] is not None:
                    outstr += '\t{0:.3f}'.format(data_indel_normal['s_ratio'])
                else:
                    outstr += '\t---'
                outstr += '\t{0:.3f}'.format(data_indel_tumor['fisher'])
                
                w.write( outstr + "\n" )


############################################################
def Pileup_and_count(
        in_bam1,
        in_bam2,
        out_file,
        ref_fa,
        baseq_thres,
        mismatch_rate_disease,
        mismatch_rate_normal,
        post_10_q,
        fisher_threshold,
        min_depth,
        print_header,
        min_variant_read,
        samtools,
        samtools_params,
        region
        ):

    global target
    global remove_chr
    global filter_quals
    
    #
    # Initalize filter quality values
    #
    filter_quals = ''
    for qual in range( 33, 33 + baseq_thres ):
        filter_quals += str( unichr( qual ) )
    filter_quals = list(filter_quals)

    #
    # Setup regular expression
    # ([\+\-])[0-9]+[ACGTNacgtn]+
    #
    target = re.compile( '([\+\-])([0-9]+)([ACGTNRMacgtnrm]+)' )
    remove_chr = re.compile( '\^.' )

    samtools_params_list = samtools_params.split(" ")

    #
    # Open output file and write header
    #
    w = open( out_file, 'w' )
    FNULL = open(os.devnull, 'w')
    #
    # Print header only for testing.
    #

    if in_bam1 and in_bam2:
        if print_header:
            header_str = "#chr\tstart\tend\tref\talt\tdepth_tumor\tvariantNum_tumor\tdepth_normal\tvariantNum_normal\tbases_tumor\tbases_normal\tA,C,G,T_tumor\tA,C,G,T_normal\tmisRate_tumor\tstrandRatio_tumor\tmisRate_normal\tstrandRatio_normal\tP-value(fisher)\n"
            w.write( header_str )
        ###Fix:0.2.1: use --output-QNAME, --reverse-del and -x
        # To use this option, samtools version is need to be >= 1.10
        # Note that when samtools is conducted without -x,
        # ignored base in ovelapped region is opposite in samtools 1.2 and 1.13
        cmd_list = [samtools,'mpileup','--output-QNAME','--reverse-del','-x','-f',ref_fa]
        cmd_list.extend(samtools_params_list)
        cmd_list.extend([in_bam1, in_bam2])
        if region:
            cmd_list.insert(2, '-r')
            cmd_list.insert(3, region)
        
        # multithreading in 1 core does not improve performance...
        
        #check whether samtools has --output-QNAME option
        pileup_check = subprocess.Popen(cmd_list[:cmd_list.index('--output-QNAME')+1],
                                             stdout=subprocess.PIPE, stderr = subprocess.PIPE)
        if 'unrecognized option' in pileup_check.stderr.read():
            logging.error( "The samtools version is old. Use >= 1.10")
            raise AssertionError
                
        pileup = subprocess.Popen(cmd_list, stdout=subprocess.PIPE, stderr = FNULL)
            
        end_of_pipe = pileup.stdout
        for mpileup in end_of_pipe:
            data = Pileup_out( mpileup, min_depth, min_variant_read, True)
            if data:
                print_data( data, w, min_depth, mismatch_rate_disease, mismatch_rate_normal, post_10_q, fisher_threshold, min_variant_read )

    elif in_bam1 or in_bam2:
        if print_header:
            header_str = "#chr\tstart\tend\tref\talt\tdepth\tvariantNum\tbases\tA,C,G,T\tmisRate\tstrandRatio\t10%_posterior_quantile\tposterior_mean\t90%_posterior_quantile\n"
            w.write( header_str )
        in_bam = in_bam1 if in_bam1 else in_bam2
        ###Fix:0.2.1: use --output-QNAME, --reverse-del and -x
        # To use this option, samtools version is need to be >= 1.10
        # Note that when samtools is conducted without -x,
        # ignored base in ovelapped region is opposite in samtools 1.2 and 1.13
        cmd_list = [samtools,'mpileup','--output-QNAME','--reverse-del','-x','-f',ref_fa]
        cmd_list.extend(samtools_params_list)
        cmd_list.extend([in_bam])
        if region:
            cmd_list.insert(2, '-r')
            cmd_list.insert(3, region)
            
        #check whether samtools has --output-QNAME option
        pileup_check = subprocess.Popen(cmd_list[:cmd_list.index('--output-QNAME')+1],
                                                      stdout=subprocess.PIPE, stderr = subprocess.PIPE)
        if 'unrecognized option' in pileup_check.stderr.read():
            logging.error( "The samtools version is old. Use >= 1.10")
            raise AssertionError
                
        pileup = subprocess.Popen(cmd_list, stdout=subprocess.PIPE, stderr = FNULL)
            
        end_of_pipe = pileup.stdout
        for mpileup in end_of_pipe:
            data = Pileup_out( mpileup, min_depth, min_variant_read, False )
            if data:
                print_data( data, w, min_depth, mismatch_rate_disease, mismatch_rate_normal, post_10_q, fisher_threshold, min_variant_read )

    else:
        logging.error( "Input file: {file} not found.".format( file = args.in_bam1 +" "+ args.in_bam2 ) )
        raise
    
    FNULL.close()
    w.close()

