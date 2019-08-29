#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 27 2:26:01 2019

@author: jhyun95
"""

from utils import process_header, get_sequences_as_dict
import os
import cobra
import numpy as np
import pandas as pd

FLUX_MINIMUM = 1e-6
FREE_METABOLITES = ['o2', 'nh4', 'h', 'h2o', 'pi', 'so4', 'co2',
                    'cl', 'ca2', 'cu2', 'cobalt2', 'fe2', 'fe3',
                    'k', 'mg2', 'mn2', 'mobd', 'ni2', 'zn2']

def allele_analysis(query_dir, ref_dir='reference/', rare_limit=1, low_mem=True):
    '''
    Same as allele_analysis_differential_carbon_source, but 
    analyzes all genes instead of those involved in differential reactions
    '''
    allele_analysis_differential_carbon_source(query_dir,
        expt_source=None, ref_dir=ref_dir, rare_limit=rare_limit, low_mem=low_mem) 


def allele_analysis_differential_carbon_source(query_dir, expt_source, 
    control_source='glc__D', ref_dir='reference/', rare_limit=-1, low_mem=True):
    ''' 
    Runs the following analysis for a given model 
    1) Simulate on an experimental carbon source and a control carbon source
    2) Identify reactions active only in the experimental conditions 
    3) Identify genes associated with the differential reactions 
    4) Compare upstream and coding sequences of those genes to reference strains

    Expects there to be four files in query-dir: <strain>.faa, <strain>.json,
    <strain>_cdhit_merged.faa.clstr, and <strain>_upstream.fna.

    If expt_source=None, does not do a differential analysis but instead 
    analyzes all genes as in steps 3/4.
    '''

    ''' Load files relevant to the query  '''
    for filename in os.listdir(query_dir):
        if filename[-4:] == '.faa':
            strain = filename[:-4]
    protein_file = query_dir + '/' + strain + '.faa'
    model_file = query_dir + '/' + strain +  '.json'
    upstream_file = query_dir + '/' + strain + '_upstream.fna'
    cluster_file = query_dir + '/' + strain + '_cdhit_merged.faa.clstr'
    for filename in [protein_file, model_file, upstream_file, cluster_file]:
        if not os.path.exists(filename):
            print 'FILE IS MISSING, ABORTING:', filename
            return  

    ''' Extract all co-clustered reference genes '''
    def query_fxn(feature_name):
        return feature_name[:5] != '>lcl|'
    query_to_cluster = get_co_clustered(cluster_file, query_fxn)
    print 'Loaded query clusters:', len(query_to_cluster)

    ''' Simulate differential growth '''
    model = cobra.io.load_json_model(model_file)
    if expt_source != None: # analyzing differential genes
        expt_diff_genes = get_differential_reactions(model, expt_source, control_source)
        print 'Found differentially active reactions:', len(expt_diff_genes)
    else: # analyzing all genes
        expt_diff_genes = get_gene_mapping_for_reactions(model, map(lambda x: x.id, model.reactions))
        print 'Analyzing all reactions:', len(expt_diff_genes)

    ''' Extract locus tags for relevant reference genes '''
    header_to_label = {} # map raw headers (without ">") to labels
    with open(ref_dir + '/ref_labels.tsv', 'r') as f:
        for line in f:
            name, label = line.split()
            header_to_label[name] = label

    matched_diff_queries = {} # query_to_cluster reduced to just differentially active genes
    diff_clustered_ref_names = set() # set of all reference genes co-clustered with a differential gene 
    for rxnID in expt_diff_genes:
        for query in expt_diff_genes[rxnID].values():
            matched_diff_queries[query] = query_to_cluster['>'+query]
            diff_clustered_ref_names = diff_clustered_ref_names.union(query_to_cluster['>'+query])
    del query_to_cluster 
    diff_clustered_ref_loci = map(lambda x: header_to_label[x[1:]].split('|')[1], diff_clustered_ref_names)
    print 'Found co-clustered reference genes:', len(diff_clustered_ref_loci)
    
    ''' Extract relevant sequences from reference and query files '''
    def filter_ref_seq(header): # identifying reference protein sequences 
        if len(header.strip()) == 0: # catch empty lines
            return False
        return header.split()[0] in diff_clustered_ref_names

    def filter_ref_upstream(header): # identifying reference upstream sequences
        if len(header.strip()) == 0: # catch empty lines
            return False
        return header.split('|')[1] in diff_clustered_ref_loci

    def filter_query_seq(header): # identifying query protein sequences
        return header[1:] in matched_diff_queries

    def filter_query_upstream(header): # identifying query upstream sequences
        label = header[1:].split('|')
        label = label[0] + '|' + label[-1]
        return label in matched_diff_queries

    ref_seqs = {} # maps ref headers ">lcl|..." to protein sequences
    ref_upstreams = {} # maps ref headers "<strain>|<tag>|..." to upstream sequences
    query_seqs = {} # maps query headers "<local gene>|<cluster gene>" to protein sequences
    query_upstreams = {} # maps query headers "<local>|<up>|<cluster gene>" to upstream sequences

    ref_seq_dir = (ref_dir + '/ref_genomes/').replace('//','/')
    for ref_seq_file in os.listdir(ref_seq_dir):
        if ref_seq_file[-4:] == '.faa': # expect amino acid fasta
            ref_seq_path = ref_seq_dir + ref_seq_file
            if low_mem: # Load only relevant sequences
                ref_seqs.update(get_sequences_as_dict(ref_seq_path, select_fxn=filter_ref_seq))
            else: # Load all sequences, faster since no filtering
                ref_seqs.update(get_sequences_as_dict(ref_seq_path))
    ref_seqs = {k.split()[0]: v for k,v in ref_seqs.iteritems()}
    print 'Loaded reference sequences for genes:', len(ref_seqs)
    
    ref_upstream_dir = (ref_dir + '/ref_upstream/').replace('//','/')
    for ref_upstream_file in os.listdir(ref_upstream_dir):
        if ref_upstream_file[-4:] == '.fna': # expect nucleotide fasta
            ref_upstream_path = ref_upstream_dir + ref_upstream_file
            if low_mem: # Load only relevant sequences
                ref_upstreams.update(get_sequences_as_dict(ref_upstream_path, select_fxn=filter_ref_upstream))
            else: # Load all sequences, faster since no filtering
                ref_upstreams.update(get_sequences_as_dict(ref_upstream_path))
    ref_upstreams = {k.split('|')[1]: v for k,v in ref_upstreams.iteritems()}
    print 'Loaded reference upstream sequences for genes:', len(ref_upstreams)

    ''' Extract relevant query sequences '''
    query_seqs = get_sequences_as_dict(protein_file, select_fxn=filter_query_seq)
    query_seqs = {k[1:]: v for k,v in query_seqs.iteritems()}
    print 'Loaded query sequences for genes:', len(query_seqs)

    query_upstreams = get_sequences_as_dict(upstream_file, select_fxn=filter_query_upstream)
    format_header = lambda x: x.split('|')[0] + '|' + x.split('|')[-1]
    query_upstreams = {format_header(k[1:]): v for k,v in query_upstreams.iteritems()}
    print 'Loaded query upstream sequences for genes:', len(query_upstreams)

    ''' Report allele analysis of differential genes '''
    for rxnID in sorted(expt_diff_genes.keys()):
        reported_reaction = False
        rxn_genes = expt_diff_genes[rxnID]
        for match_gene in sorted(rxn_genes.keys()):
            query_gene = rxn_genes[match_gene]
            query_seq = query_seqs[query_gene]
            query_ups = query_upstreams[query_gene][:53]
            
            ''' Get allele distribution of reference gene/upstream sequences '''
            co_clustered = matched_diff_queries[query_gene]
            seq_distr = {}; ups_distr = {}
            for ref_gene in co_clustered:
                ref_tag = header_to_label[ref_gene[1:]].split('|')[1]
                if ref_tag in ref_upstreams and ref_gene in ref_seqs:
                    # exclude rare cases where either piece of information is missing
                    ref_seq = ref_seqs[ref_gene]
                    ref_ups = ref_upstreams[ref_tag][:53]
                    if not ref_seq in seq_distr:
                        seq_distr[ref_seq] = 0
                    if not ref_ups in ups_distr:
                        ups_distr[ref_ups] = 0
                    seq_distr[ref_seq] += 1
                    ups_distr[ref_ups] += 1

            ''' Add in query sequence '''
            if not query_seq in seq_distr:
                seq_distr[query_seq] = 0
            if not query_ups in ups_distr:
                ups_distr[query_ups] = 0
            seq_distr[query_seq] += 1
            ups_distr[query_ups] += 1

            ''' Report allele distribution '''
            query_seq_count = seq_distr[query_seq]
            query_ups_count = ups_distr[query_ups]
            if rare_limit < 0 or (query_seq_count <= rare_limit and query_ups_count <= rare_limit):
                if not reported_reaction:
                    print '\n------------- REACTION:', rxnID, '-------------'
                    reported_reaction = True
                print '\nMATCH:', match_gene , '<->', query_gene
                print 'Coding sequence distribution:'
                print 'Count\tLength\tSeq'
                for seq_allele in sorted(seq_distr.keys()):
                    count = str(seq_distr[seq_allele])
                    if seq_allele == query_seq:
                        count += '*'
                    print count, '\t', len(seq_allele), '\t', seq_allele[:50] + '...'

                print '\nUpstream sequence distribution:'
                print 'Count\tSeq'
                for ups_allele in sorted(ups_distr.keys()):
                    count = str(ups_distr[ups_allele])
                    if ups_allele == query_ups:
                        count += '*'
                    print count, '\t', ups_allele


def get_co_clustered(query_cluster_file, query_fxn):
    ''' 
    Given a CD-Hit cluster file, attempts to identify query sequences 
    (using query_fxn) and all non-query sequences clustered with them.
    Returns a dictionary {query_name:set(co-clustered sequence names)}
    '''
    clusters = {}; current_cluster = set(); current_query = ''
    with open(query_cluster_file, 'r') as f:
        for line in f:
            if line[0] == '>': # new cluster
                if len(current_cluster) > 0 and len(current_query) >0: # if ended last cluster
                    clusters[current_query] = current_cluster
                current_cluster = set(); current_query = ''
            else: # reading current cluster
                feature_name = line.split()[2].replace('...','')
                if query_fxn(feature_name): # if query sequence
                    current_query = feature_name
                else:
                    current_cluster.add(feature_name)
        if len(current_cluster) > 0 and len(current_query) >0: # final cluster
            clusters[current_query] = current_cluster
    return clusters


def get_differential_reactions(model, expt_source, control_source='glc__D'):
    ''' For a given model, simulates growth with an experimental carbon source
        and a reference carbon source (i.e. glucose) using pFBA. Fluxes that
        are active only under the experimental carbon source are reported,
        along with all associated genes.  '''
    model = set_model_carbon_source(model, expt_source)
    expt_soln = cobra.flux_analysis.parsimonious.pfba(model)
    expt_growth = expt_soln.objective_value
    expt_fluxes = expt_soln.fluxes

    if expt_growth < FLUX_MINIMUM: # no growth under experimental conditions
        return None
    else: # growth under experimental conditions
        ''' Repeat under control conditions '''
        model = set_model_carbon_source(model, control_source)
        ctrl_soln = cobra.flux_analysis.parsimonious.pfba(model)
        ctrl_growth = ctrl_soln.objective_value
        ctrl_fluxes = ctrl_soln.fluxes

        ''' Identify reactions active in experimental, inactive in control '''
        expt_active = np.abs(expt_fluxes.values) > FLUX_MINIMUM
        ctrl_inactive = np.abs(ctrl_fluxes.values) < FLUX_MINIMUM
        expt_difference = np.logical_and(expt_active, ctrl_inactive)
        expt_diff_rxns = expt_fluxes.index[expt_difference]

        ''' Get genes for each enzyme-mediated reaction '''
        expt_diff_genes = get_gene_mapping_for_reactions(model, expt_diff_rxns)
        return expt_diff_genes


def get_gene_mapping_for_reactions(model, rxns):
    ''' Get genes associated with a list of reactions for a reconstruction.
        Returns a nested dictionary mapping {rxn.id: geneID: match} '''
    genes_map = {}
    for rxnID in rxns:
        rxn = model.reactions.get_by_id(rxnID)
        gpr = rxn.gene_reaction_rule + ''
        for operator in ['(', ')', ' and ', ' or ']:
            gpr = gpr.replace(operator, '  ')
        rxn_genes = set(gpr.split())
        if len(rxn_genes) > 0 and not 's0001' in rxn_genes:
            genes_map[rxn.id] = {}
            for geneID in rxn_genes:
                gene = model.genes.get_by_id(geneID)
                if 'match' in gene.notes:
                    genes_map[rxn.id][geneID] = gene.notes['match']
    return genes_map


def set_model_carbon_source(model, carbon_source, limit=-10.0):
    ''' For a given model, disables all input exchanges, then enables 
        free exchange of all metabolites in FREE_METABOLITES, and finally
        enables uptakes of the specific carbon source up to a specified limit '''
    for reaction in model.reactions.query('EX_'): # get all exchanges
        if reaction.id.split('_')[1] in FREE_METABOLITES:
            reaction.lower_bound = -1000.0
            reaction.upper_bound = 1000.0
        else: # everything else export-only
            reaction.lower_bound = 0.0
            reaction.upper_bound = 1000.0
    exchangeID = 'EX_' + carbon_source + '_e' if not 'EX_' in carbon_source else carbon_source
    exchange = model.reactions.get_by_id(exchangeID)
    exchange.lower_bound = limit
    return model