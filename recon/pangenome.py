#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 10 17:27:00 2019

@author: jhyun95

Various tools for analyzing pangenomes.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as sch

from utils import get_sequences_as_dict

def get_duplicates(gff_file, verbose=True):
    ''' Identify genes annotated to the same location '''
    features = {}; duplicate_features = {}
    with open(gff_file, 'r') as f:
        for line in f:
            if line[:2] != '##':
                entries = line.split('\t')
                if entries[2] == 'CDS':
                    start, stop = entries[3:5]
                    strand = entries[6]
                    gene_location = (start,stop,strand)
                    for annot in entries[-1].split(';'):
                        param, value = annot.split('=')
                        if param == 'locus_tag':
                            if not gene_location in features:
                                features[gene_location] = [value]
                            else:
                                features[gene_location].append(value)
                                duplicate_features[gene_location] = features[gene_location]
                                if verbose:
                                    print 'CONFLICT:', features[gene_location], gene_location
            elif '##FASTA' in line:
                break
    return duplicate_features


def build_gene_table(clstr_file, faa_files, cluster_log=1000):
    ''' Builds binary gene table from a CD-Hit cluster file and 
        the original FAA files that were provided for clustering.
        Infers strain names from FAA filenames '''

    ''' Assign headers to strains '''
    print 'Extract feature names for strains...'
    feature_to_strain = {}; strain_order = []
    for faa_file in faa_files:
        strain = faa_file.split('/')[-1].split('.')[0]
        strain_order.append(strain)
        with open(faa_file, 'r') as f:
            for line in f:
                if '>' == line[0]:
                    feature_name = line.split()[0]
                    feature_to_strain[feature_name] = strain

    ''' Parse cluster file ''' 
    print 'Parsing cluster file...'
    gene_data = []
    current_cluster = np.zeros(len(strain_order), dtype=int)
    cluster_names = []
    with open(clstr_file, 'r') as f:
        for line in f:
            if line[0] == '>': # starting new cluster
                if len(cluster_names) > len(gene_data):
                    gene_data.append(current_cluster)
                    if len(gene_data) % cluster_log == 0:
                        print 'Clusters processed', len(gene_data)
                cluster_name = line[1:].replace(' ','_').strip()
                cluster_names.append(cluster_name)
                current_cluster = np.zeros(len(strain_order), dtype=int)
            else: # reading current cluster
                feature_name = line.split()[2].replace('...','').strip()
                try:
                    feature_strain = feature_to_strain[feature_name]
                    current_cluster[strain_order.index(feature_strain)] = 1
                except KeyError:
                    print 'Orphan feature:', feature_name
    if len(cluster_names) > len(gene_data): # process last cluster 
        gene_data.append(current_cluster)

    df_genes = pd.DataFrame(index=cluster_names, columns=strain_order, data=np.array(gene_data))
    df_genes.replace(0, np.nan, inplace=True) # replace with NaN to reduce filesize
    return df_genes

def build_allele_table(clstr_file, faa_files, seq_out, cluster_log=1000):
    ''' Builds binary allele table from a CD-Hit cluster file and 
        the original FAA files that were provided for clustering.
        Also exports non-redundant unique alleles.
        Infers strain names from FAA filenames '''
    print 'Loading sequences...'
    all_seqs = {} # feature_name:(seq,strain)
    for faa_file in faa_files:
        strain = faa_file.split('/')[-1].split('.')[0]
        strain_seqs = get_sequences_as_dict(faa_file)
        strain_seqs = {k.split()[0]:(v,strain) for k,v in strain_seqs.items()}
        all_seqs.update(strain_seqs)

    print 'Parsing cluster file...'
    cluster_features = []; cluster_name = ''
    allele_table = {}; clusters = 0
    f_seq = open(seq_out, 'w+')
    with open(clstr_file, 'r') as f:
        for line in f:
            if line[0] == '>': # starting new cluster
                ''' Processing finished cluster '''
                if len(cluster_features) > 0: 
                    alleles = {}; clusters += 1
                    for feature_name in cluster_features: 
                        ''' Identify and name unique alleles '''
                        seq, strain = all_seqs[feature_name]
                        if not seq in alleles: # new allele encountered
                            allele_name = cluster_name + 'A' + str(len(alleles))
                            f_seq.write('>' + allele_name + '\n')
                            f_seq.write(seq + '\n')
                            alleles[seq] = allele_name

                        ''' Updated allele table'''
                        allele_name = alleles[seq]
                        if not strain in allele_table:
                            allele_table[strain] = {}
                        allele_table[strain][allele_name] = 1

                ''' Logging '''
                if clusters % cluster_log == 0:
                    print 'Clusters processed', clusters

                ''' Initializing new cluster '''
                cluster_name = line[1:].replace('Cluster ','C').strip()
                cluster_features = []
            else: # reading current cluster
                feature_name = line.split()[2].replace('...','').strip()
                cluster_features.append(feature_name)

    ''' Processing last cluster '''
    if len(cluster_features) > 0: 
        alleles = {}; clusters += 1
        for feature_name in cluster_features: 
            ''' Identify and name unique alleles '''
            seq, strain = all_seqs[feature_name]
            if not seq in alleles:
                allele_name = cluster_name + 'A' + str(len(alleles))
                alleles[seq] = allele_name

            ''' Updated allele table'''
            allele_name = alleles[seq]
            if not strain in allele_table:
                allele_table[strain] = {}
            allele_table[strain][allele_name] = 1

    df_alleles = pd.DataFrame.from_dict(allele_table)
    return df_alleles


def plot_gene_abundance(df_genes, plot_both=True):
    ''' Generate a plot of gene abundance vs. count '''

    num_strains = df_genes.shape[1]
    gene_counts = df_genes.sum(axis=1).values

    ''' Get gene frequencies '''
    frequencies = np.arange(num_strains)
    freq_counts = np.zeros(num_strains) # index i = # of genes present in (i+1) strains
    for i in frequencies:
        freq_counts[i] = (gene_counts == i+1).sum()
    cumulative_counts_unique_first = np.cumsum(freq_counts)
    cumulative_counts_core_first = np.cumsum(np.flip(freq_counts, axis=0))

    ''' Plot cumulative counts up to frequency '''
    print 'Core:', freq_counts[-1]
    print 'Unique:', freq_counts[0]
    plt.plot(cumulative_counts_unique_first, frequencies)
    if plot_both:
        plt.plot(cumulative_counts_core_first, np.flip(frequencies, axis=0))
    plt.xlabel('# genes')
    plt.ylabel('# strains')


def plot_strain_hierarchy(df_genes, df_labels=None, df_colors=None, 
    linkage='average', figsize=(8,6), polar=False):
    ''' Cluster strains by gene content, using Jaccard distance '''
    X = df_genes.fillna(0).values.T # flip to strain x gene cluster
    print X.shape
    num_strains, num_genes = X.shape

    ''' Compute pairwise Jaccard distances '''
    distances = np.zeros((num_strains,num_strains))
    for i in range(num_strains):
        genes1 = X[i,:]
        for j in range(i):
            genes2 = X[j,:]
            shared = np.logical_and(genes1,genes2).sum()
            union = np.logical_or(genes1,genes2).sum()
            dist = 1.0 - (shared / np.float(union))
            distances[i,j] = dist
            distances[j,i] = dist
            #print df_genes.columns[i], df_genes.columns[j], dist

    ''' Generate dendrogram '''
    dist_condensed = distances[np.triu_indices(num_strains,1)]
    Z = sch.linkage(dist_condensed, linkage)

    if df_labels is None:
        labels = df_genes.columns
    else:
        for strain in df_genes.columns:
            if not strain in df_labels.index:
                print strain, 'missing'
        labels = df_labels.reindex(df_genes.columns)
    
    if polar: 
        dend = sch.dendrogram(Z, labels=labels, no_plot=True)
        plot_polar_dendogram(dend, figsize=figsize, df_colors=df_colors)
    else:
        # TODO: Colored text labels for non-polar dendrogram
        fig, ax = plt.subplots(1,1, figsize=figsize)
        dend = sch.dendrogram(Z, labels=labels, ax=ax)
    df_distances = pd.DataFrame(index=df_genes.columns, columns=df_genes.columns, data=distances)
    return dend, df_distances


def plot_polar_dendogram(dend, figsize=(8,8), df_colors=None):
    ''' Adapted from https://stackoverflow.com/questions/51936574/how-to-plot-scipy-hierarchy-dendrogram-using-polar-coordinates'''
    icoord = np.array(dend['icoord'])
    dcoord = np.array(dend['dcoord'])
    
    ''' Transformations for polar coordinates + formatting '''
    def smooth_segment(seg, Nsmooth=100):
        return np.concatenate([[seg[0]], np.linspace(seg[1], seg[2], Nsmooth), [seg[3]]])

    gap = 0.1
    dcoord = -dcoord
    #dcoord = -np.log(dcoord+1) # log transform, makes polar dendrogram look better
    imax = icoord.max()
    imin = icoord.min()
    icoord = ((icoord - imin)/(imax - imin)*(1-gap) + gap/2)*2*np.pi
    
    ''' Plotting polar dendrogram '''
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, polar=True)
    for xs,ys in zip(icoord, dcoord):
        xs = smooth_segment(xs)
        ys = smooth_segment(ys)
        ax.plot(xs, ys, color="black")
    
    ''' Adjust polar tick labels '''
    ax.spines['polar'].set_visible(False)
    ax.set_rlabel_position(0)
    num_xticks = dcoord.shape[0]+1
    angles = np.linspace(gap/2, 1-gap/2, num_xticks) * np.pi * 2
    ax.set_xticks(angles) #*np.pi*2)
    ax.set_xticklabels(dend['ivl'])

    plt.gcf().canvas.draw()
    for label, angle in zip(ax.get_xticklabels(), angles):
        x,y = label.get_position()
        label_text = label.get_text()
        shift = 0.007 * len(label_text)
        if not df_colors is None:
            color = df_colors.loc[label_text][0]
            #print label_text, color
        else:
            color = 'black'
        lab = ax.text(x,y-shift, label.get_text(), transform=label.get_transform(),
            ha=label.get_ha(), va=label.get_va(), color=color)
        
        if angle > 0.5*np.pi and angle < 1.5*np.pi:
            angle += np.pi
        lab.set_rotation(angle * 180.0 / np.pi)

    ax.set_xticklabels([])

def assembly_stats(fna_file, percent=50):
    ''' Reports # of contigs (L), sum of contig length (N), N50, and L50.
        Can compute other N* and L* statistics by changing the percent arg. '''

    contig_lengths = get_sequences_as_dict(fna_file, apply_fxn=len)
    L = len(contig_lengths) # num contigs
    N = sum(contig_lengths.values()) # sum of contig lengths
    N_run = 0; N_threshold = float(percent / 100.0) * N
    for L_run, contig_length in enumerate(sorted(contig_lengths.values(), reverse=True)):
        N_run += contig_length
        if N_run > N_threshold:
            break
    Lstat = L_run + 1; Nstat = contig_length
    return (L,Lstat,N,Nstat)
