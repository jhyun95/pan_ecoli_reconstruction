#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 27 4:27:55 2019

@author: jhyun95
"""

import numpy as np

def process_header(header):
    ''' From an NCBI sequence header, extracts the raw name, annotations, and label
        of the form <strain>|<locus_tag>|<protein_id> '''
    raw_name = header.split()[0][1:]
    strain_name = header.split('_prot_')[0][5:]
    parameters = {}
    for annotation in header.split()[1:]:
        if '=' in annotation:
            annot_split = annotation.strip()[1:-1].split('=')
            if len(annot_split) == 2:
                param, value = annot_split
                parameters[param] = value
    locus_tag = parameters['locus_tag'] if 'locus_tag' in parameters else 'unknown'
    protein_id = parameters['protein_id'] if 'protein_id' in parameters else 'unknown'
    if locus_tag == 'unknown' and 'gene' in parameters and strain_name == 'BA000007.2':
        locus_tag = parameters['gene']
    label = strain_name + '|' + locus_tag + '|' + protein_id
    return raw_name, parameters, label


def get_sequences_as_dict(fasta_file, select_fxn=None, apply_fxn=None):
    ''' Adapted from scripts/utils.py. Loads sequences as a dict {header:sequence},
        can optionally filter by header or transform sequences.

        Parameters
        ----------
        fasta_file : str
            Path to fasta file (amino acid or nucleotide)
        select_fxn : function
            If not None, keeps only sequences where the corresponding header 
            satsifies the function (default None)
        apply_fxn : function
            If not None, applies the function each sequence, such as computing
            length, GC content, reverse complement, etc. (default None)

        Returns
        -------
        seqs : dict
            Dictionary mapping headers (including initial ">") to sequences
            or as transformed by apply_fxn when provided.
        '''

    header = ''; seq = ''; seqs = {}
    with open(fasta_file, 'r') as f:
        for line in f:
            if line[0] == '>': # header line
                satisfies_fxn = select_fxn is None or select_fxn(header) # check select_fxn
                if satisfies_fxn and len(header) > 0:
                    seq = apply_fxn(seq) if apply_fxn else seq # check apply_fxn
                    seqs[header] = seq
                header = line.strip(); seq = ''
            else: # sequence line
                seq += line.strip()

    ''' Repeat for last entry in file '''
    satisfies_fxn = select_fxn is None or select_fxn(header) # check select_fxn
    if satisfies_fxn and len(header) > 0: 
        seq = apply_fxn(seq) if apply_fxn else seq
        seqs[header] = seq
    return seqs


def edit_distance(s1, s2, substitution=1, indel=1, limit=-1):
    ''' 
    Adapted from https://stackoverflow.com/questions/2460177/edit-distance-in-python 
    Simple dynamic programming approach for computing edit distance. Optionally, can 
    terminate the computation early if edit distance is found to be greater than some limit.
    
    Args:
        s1 : str
            First sequence in edit distance pair
        s2 : str
            Second sequence in edit distance pair
        substitution : int
            Penalty for substitution error (default 1)
        indel : int
            Penalty for insertion or deletion (default 1)
        limit : int
            If positive and edit distance exceeds limit, then
            terminates early and returns limit (default -1)
            
    Returns:
        edit_dist : int
            Computed edit distance between the two sequences
    '''
    m = len(s1) + 1; n = len(s2) + 1
    tbl = np.zeros((m,n))
    
    for i in np.arange(m): tbl[i,0] = i * indel
    for j in np.arange(n): tbl[0,j] = j * indel
    for i in np.arange(1, m):
        for j in np.arange(1, n):
            cost = 0 if s1[i-1] == s2[j-1] else substitution
            tbl[i,j] = min(tbl[i, j-1] + indel, tbl[i-1, j] + indel, tbl[i-1, j-1] + cost)
        if limit > 0 and np.min(tbl[i,:]) > limit:
            return limit
    edit_dist = tbl[i,j]
    return edit_dist