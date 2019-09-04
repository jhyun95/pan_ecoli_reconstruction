#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 27 2:26:01 2019

@author: jhyun95

Various tools for finer refinement of draft reconstruction. High level methods include:

add_biolog_reactions() # reactions for arbutin, salicin, and D-tagatose metabolism
allele_analysis() # compute edit distances of aligned genes to their reference clusters
load_model_with_rare_filter() # load model without genes with large distance (from allele_analysis)

"""

from utils import process_header, get_sequences_as_dict, edit_distance, hamming_distance
import os
import cobra
import numpy as np
import pandas as pd

FLUX_MINIMUM = 1e-6
FREE_METABOLITES = ['o2', 'nh4', 'h', 'h2o', 'pi', 'so4', 'co2',
                    'cl', 'ca2', 'cu2', 'cobalt2', 'fe2', 'fe3',
                    'k', 'mg2', 'mn2', 'mobd', 'ni2', 'zn2']

def add_biolog_reactions(model):
    ''' Adds sink reactions for arbutin and salicin, as well as other exchanges
        for biolog-tested substrate. Specifically, adds hydroquinone and 
        salicyl alcohol sinks for arbutin and salicin metabolism, respectively. '''
    if not 'DM_hqn_c' in model.reactions: # hydroquinone sink for arbutin
        if 'hqn_c' in model.metabolites: # check if the model supports arbutin metabolism
            met_hqn = model.metabolites.get_by_id('hqn_c')
            rxn = cobra.core.Reaction('DM_hqn_c')
            rxn.name = 'Sink needed to allow hydroquinone to leave system'
            rxn.subsystem = 'Intracellular demand'
            rxn.lower_bound = 0.0; rxn.upper_bound = 1000.0
            rxn.add_metabolites({met_hqn:-1})
            model.add_reactions([rxn])
    if not 'DM_2hymeph_c' in model.reactions: # salicyl alcohol sink for arbutin
        if '2hymeph_c' in model.metabolites: # check if the model supports salicin metabolism
            met_2hymeph = model.metabolites.get_by_id('2hymeph_c') # aka 2-(Hydroxymethyl)phenol
            rxn = cobra.core.Reaction('DM_2hymeph_c')
            rxn.name = 'Sink needed to allow 2-(Hydroxymethyl)phenol to leave system'
            rxn.subsystem = 'Intracellular demand'
            rxn.lower_bound = 0.0; rxn.upper_bound = 1000.0
            rxn.add_metabolites({met_2hymeph:-1})
            model.add_reactions([rxn])

    biolog_exchanges = { # tested substrates missing in iML1515
        'arab__D': 'D Arabinose exchange',
        'rbt': 'Ribitol exchange',
        'salcn': 'Salicin exchange',
        'raffin': 'Raffinose exchange',
        'tag__D': 'D Tagatose exchange',
        'abt__D': 'D-Arabitol exchange'
    }
    new_reactions = []
    for metID in biolog_exchanges:
        exchange_missing = not (('EX_' + metID + '_e') in model.reactions)
        metabolite_present = (metID + '_e') in model.metabolites
        if exchange_missing and metabolite_present:
            met = model.metabolites.get_by_id(metID + '_e')
            rxn = cobra.core.Reaction('EX_' + metID + '_e')
            rxn.name = biolog_exchanges[metID]
            rxn.subsystem = 'Extracellular exchange'
            rxn.lower_bound = 0.0; rxn.upper_bound = 1000.0
            rxn.add_metabolites({met:-1})
            new_reactions.append(rxn)
    if len(new_reactions) > 0:
        model.add_reactions(new_reactions)
    return model


def load_model_with_rare_filter(query_dir, seq_dist_limit=30, upstream_dist_limit=10):
    ''' Loads the reconstruction model and disables genes that are "extreme", 
        i.e. edit distance to nearest reference gene by protein sequence or 
        upstream sequence exceeds the specified limit '''
    for filename in os.listdir(query_dir):
        if filename[-4:] == '.faa':
            strain = filename[:-4]
    extreme_alleles = get_rare_allele_genes(query_dir, seq_dist_limit=seq_dist_limit,
        upstream_dist_limit=upstream_dist_limit)
    extreme_genes = extreme_alleles.keys()
    model = cobra.io.load_json_model(query_dir + '/' + strain +  '.json')
    if len(extreme_genes) > 0:
        cobra.manipulation.delete.delete_model_genes(model, extreme_genes)
        print 'Disabled', len(extreme_genes), 'extreme genes:'
        for gene in extreme_genes:
            print '\t', gene
    else:
        print 'No extreme genes detected.'
    return model


def get_rare_allele_genes(query_dir, seq_dist_limit=10, upstream_dist_limit=5):
    '''
    Identifies and disables genes that are extreme alleles (edit distance
    to nearest reference gene in terms of coding or upstream sequence is 
    above either of the two thresholds provided). Requires allele_analysis()
    to have been run first.  
    '''
    for filename in os.listdir(query_dir):
        if filename[-4:] == '.faa':
            strain = filename[:-4]
    model_file = query_dir + '/' + strain +  '.json'
    log_file = query_dir + '/' + strain + '_allele_report.txt'

    ''' Identify extreme alleles '''
    extreme_alleles = {} # gene to annotations
    with open(log_file, 'r') as f:
        query_gene = ''; ref_gene = ''; rxns = []; seq_dist = 0; ups_dist = 0
        for line in f:
            if 'GENE:' in line:
                is_extreme = seq_dist >= seq_dist_limit or ups_dist >= upstream_dist_limit
                if len(ref_gene) > 0 and is_extreme:
                    extreme_alleles[ref_gene] = (seq_dist, ups_dist, query_gene, rxns)
                ref_gene = line.split()[3]
                query_gene = line.split()[-2]
                seq_dist = 0; ups_dist = 0; rxns = []
            elif 'Impacted Reactions:' in line:
                rxns = line.split()[2:]
            elif 'Unique sequence, distance to nearest neighbor:' in line:
                seq_dist = int(line.split()[-1])
            elif 'Unique upstream, distance to nearest neighbor:' in line:
                ups_dist = int(line.split()[-1])
    is_extreme = seq_dist >= seq_dist_limit or ups_dist >= upstream_dist_limit
    if len(ref_gene) > 0 and is_extreme:
        extreme_alleles[ref_gene] = (seq_dist, ups_dist, query_gene, rxns)
    return extreme_alleles


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
    log_file = query_dir + '/' + strain + '_allele_report.txt'
    protein_file = query_dir + '/' + strain + '.faa'
    model_file = query_dir + '/' + strain +  '.json'
    upstream_file = query_dir + '/' + strain + '_upstream.fna'
    cluster_file = query_dir + '/' + strain + '_cdhit_merged.faa.clstr'
    for filename in [protein_file, model_file, upstream_file, cluster_file]:
        if not os.path.exists(filename):
            print 'FILE IS MISSING, ABORTING:', filename
            return  

    ''' Prepare output file '''
    log_f = open(log_file, 'w+')
    def log_to_file(*argv):
        line = ' '.join(map(str, argv))
        print line # print to console first
        log_f.write(line + '\n') # write to file

    ''' Extract all co-clustered reference genes '''
    def query_fxn(feature_name):
        return feature_name[:5] != '>lcl|'
    query_to_cluster = get_co_clustered(cluster_file, query_fxn)
    log_to_file('Loaded query clusters:', len(query_to_cluster))

    ''' Simulate differential growth '''
    model = cobra.io.load_json_model(model_file)
    if expt_source != None: # analyzing differential genes
        expt_diff_genes = get_differential_reactions_fva(model, expt_source, control_source)
        # for i in range(iters - 1):
        #     diff_genes = get_differential_reactions(model, expt_source, control_source)
        #     for gene in expt_diff_genes.keys(): # only records genes that are differential across multiple runs
        #         if not gene in diff_genes:
        #             print 'Marginal diff:', gene
        #             del expt_diff_genes[gene]
        log_to_file('Found differentially active genes:', len(expt_diff_genes))
    else: # analyzing all genes
        expt_diff_genes = get_gene_mapping_for_reactions(model, map(lambda x: x.id, model.reactions))
        log_to_file('Analyzing all genes:', len(expt_diff_genes))

    ''' Extract locus tags for relevant reference genes '''
    header_to_label = {} # map raw headers (without ">") to labels
    with open(ref_dir + '/ref_labels.tsv', 'r') as f:
        for line in f:
            name, label = line.split()
            header_to_label[name] = label

    matched_diff_queries = {} # query_to_cluster reduced to just differentially active genes
    diff_clustered_ref_names = set() # set of all reference genes co-clustered with a differential gene 
    for match_gene in expt_diff_genes:
        query_gene = expt_diff_genes[match_gene][0]
        matched_diff_queries[query_gene] = query_to_cluster['>'+query_gene]
        diff_clustered_ref_names = diff_clustered_ref_names.union(query_to_cluster['>'+query_gene])
    del query_to_cluster 
    diff_clustered_ref_loci = map(lambda x: header_to_label[x[1:]].split('|')[1], diff_clustered_ref_names)
    log_to_file('Found co-clustered reference genes:', len(diff_clustered_ref_loci))
    
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
    log_to_file('Loaded reference sequences for genes:', len(ref_seqs))
    
    ref_upstream_dir = (ref_dir + '/ref_upstream/').replace('//','/')
    for ref_upstream_file in os.listdir(ref_upstream_dir):
        if ref_upstream_file[-4:] == '.fna': # expect nucleotide fasta
            ref_upstream_path = ref_upstream_dir + ref_upstream_file
            if low_mem: # Load only relevant sequences
                ref_upstreams.update(get_sequences_as_dict(ref_upstream_path, select_fxn=filter_ref_upstream))
            else: # Load all sequences, faster since no filtering
                ref_upstreams.update(get_sequences_as_dict(ref_upstream_path))
    ref_upstreams = {k.split('|')[1]: v for k,v in ref_upstreams.iteritems()}
    log_to_file('Loaded reference upstream sequences for genes:', len(ref_upstreams))

    ''' Extract relevant query sequences '''
    query_seqs = get_sequences_as_dict(protein_file, select_fxn=filter_query_seq)
    query_seqs = {k[1:]: v for k,v in query_seqs.iteritems()}
    log_to_file('Loaded query sequences for genes:', len(query_seqs))

    query_upstreams = get_sequences_as_dict(upstream_file, select_fxn=filter_query_upstream)
    format_header = lambda x: x.split('|')[0] + '|' + x.split('|')[-1]
    query_upstreams = {format_header(k[1:]): v for k,v in query_upstreams.iteritems()}
    log_to_file('Loaded query upstream sequences for genes:', len(query_upstreams))

    ''' Report allele analysis of differential genes '''
    extreme_cases = 0
    for match_gene in sorted(expt_diff_genes.keys()):
        reported_gene = False
        query_gene = expt_diff_genes[match_gene][0]
        impacted_reactions = expt_diff_genes[match_gene][1]
        if query_gene != None:
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
                extreme_cases += 1
                if not reported_gene:
                    log_to_file('\n-------------', extreme_cases, 'GENE:', match_gene, '<->', query_gene, '-------------')
                    log_to_file('\nImpacted Reactions:', ', '.join(impacted_reactions), '\n')
                    reported_gene = True

                ''' If sequence is unique, compute distances to nearest sequence in cluster '''
                seq_neighbor = ''; ups_neighbor = ''
                if query_seq_count == 1: # query sequence is new
                    min_seq_dist = len(query_seq); 
                    for seq in seq_distr: # start with hamming distance as quick estimate
                        if seq != query_seq:
                            dist = hamming_distance(seq, query_seq)
                            if dist < min_seq_dist:
                                min_seq_dist = dist; seq_neighbor = seq
                    if min_seq_dist > 2: # if hamming distance is large, compute edit distance in full
                        for seq in seq_distr: 
                            if seq != query_seq:
                                dist = edit_distance(seq, query_seq)
                                if dist < min_seq_dist:
                                    min_seq_dist = dist; seq_neighbor = seq
                    min_seq_dist = int(min_seq_dist)
                    log_to_file('Unique sequence, distance to nearest neighbor:', min_seq_dist)

                if query_ups_count == 1: # query upstream sequence is new
                    min_ups_dist = len(query_ups); 
                    for ups in ups_distr:
                        if ups != query_ups:
                            dist = hamming_distance(ups, query_ups)
                            if dist < min_ups_dist:
                                min_ups_dist = dist; ups_neighbor = ups
                    if min_ups_dist > 2: # if hamming distance is large, compute edit distance in full
                        for ups in ups_distr:
                            if ups != query_ups:
                                dist = edit_distance(ups, query_ups)
                                if dist < min_ups_dist:
                                    min_ups_dist = dist; ups_neighbor = ups
                    min_ups_dist = int(min_ups_dist)
                    log_to_file('Unique upstream, distance to nearest neighbor:', min_ups_dist)

                ''' Report previews of co-clustered sequences/upstreams '''
                log_to_file('\nCoding sequence distribution:')
                log_to_file('Count\tLength\tSeq')
                for seq_allele in sorted(seq_distr.keys()):
                    count = str(seq_distr[seq_allele])
                    if seq_allele == query_seq:
                        count += '*'
                    elif seq_allele == seq_neighbor:
                        count += '^'
                    log_to_file(count, '\t', len(seq_allele), '\t', seq_allele[:50] + '...')

                log_to_file('\nUpstream sequence distribution:')
                log_to_file('Count\tSeq')
                for ups_allele in sorted(ups_distr.keys()):
                    count = str(ups_distr[ups_allele])
                    if ups_allele == query_ups:
                        count += '*'
                    elif ups_allele == ups_neighbor:
                        count += '^'
                    log_to_file(count, '\t', ups_allele)


def get_co_clustered(query_cluster_file, query_fxn):
    ''' 
    Given a CD-Hit cluster file, attempts to identify query sequences 
    (using query_fxn) and all non-query sequences clustered with them.
    Returns a dictionary {query_name:set(co-clustered sequence names)}
    '''
    clusters = {}; current_cluster = set(); current_queries = set()
    with open(query_cluster_file, 'r') as f:
        for line in f:
            if line[0] == '>': # new cluster
                if len(current_cluster) > 0 and len(current_queries) > 0: # if ended last cluster
                    for query in current_queries:
                        clusters[query] = current_cluster
                current_cluster = set(); current_queries = set()
            else: # reading current cluster
                feature_name = line.split()[2].replace('...','')
                if query_fxn(feature_name): # if query sequence
                    current_queries.add(feature_name)
                else:
                    current_cluster.add(feature_name)
        if len(current_cluster) > 0 and len(current_queries) > 0: # for last entry
            for query in current_queries:
                clusters[query] = current_cluster
    return clusters

def get_differential_reactions_fva(model, expt_source, control_source='glc__D', opt_frac=0.95):
    ''' For a given model, simulates growth with an experimental carbon source
        and a reference carbon source (i.e. glucose) using FVA. Then, identifies 
        which reactions have non-zero minimum absolute fluxes under only the 
        experimental conditions. '''
    model = set_model_carbon_source(model, expt_source)
    expt_growth = model.slim_optimize()
    if expt_growth < FLUX_MINIMUM: # no growth under experimental conditions
        return None
    else: # growth under experimental conditions
        ''' Compute FVA under both conditions '''
        df_expt_fva = cobra.flux_analysis.variability.flux_variability_analysis(model, 
            fraction_of_optimum=opt_frac)
        model = set_model_carbon_source(model, control_source)
        df_ctrl_fva = cobra.flux_analysis.variability.flux_variability_analysis(model, 
            fraction_of_optimum=opt_frac)
        
        ''' Get obligate reactions '''
        df_expt_fva['Required'] = np.logical_or(
            df_expt_fva['minimum'].values > FLUX_MINIMUM, 
            df_expt_fva['maximum'].values < -FLUX_MINIMUM)
        df_ctrl_fva['Required'] = np.logical_or(
            df_ctrl_fva['minimum'].values > FLUX_MINIMUM, 
            df_ctrl_fva['maximum'].values < -FLUX_MINIMUM)
        df_expt_req = df_expt_fva[df_expt_fva['Required']]
        df_ctrl_req = df_ctrl_fva[df_ctrl_fva['Required']]
        
        ''' Get differential reactions '''
        expt_req_rxns = df_expt_req.index.tolist()
        ctrl_req_rxns = df_ctrl_req.index.tolist()
        expt_diff_rxns = filter(lambda x: not x in ctrl_req_rxns, expt_req_rxns)
        expt_diff_genes = get_gene_mapping_for_reactions(model, expt_diff_rxns)
        return expt_diff_genes


def get_differential_reactions(model, expt_source, control_source='glc__D'):
    ''' For a given model, simulates growth with an experimental carbon source
        and a reference carbon source (i.e. glucose) using pFBA. Fluxes that
        are active only under the experimental carbon source are reported,
        along with all associated genes. '''
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
        Returns a nested dictionary mapping {geneID: (match, affected rxns)} '''
    genes_of_interest = set()
    for rxnID in rxns: # first get all genes involved in these reactions 
        rxn = model.reactions.get_by_id(rxnID)
        gpr = rxn.gene_reaction_rule + ''
        for operator in ['(', ')', ' and ', ' or ']:
            gpr = gpr.replace(operator, '  ')
        rxn_genes = set(gpr.split())
        genes_of_interest = genes_of_interest.union(rxn_genes)

    rxn_to_genes = {}
    for rxn in model.reactions: # next, map each reaction to its genes
        gpr = rxn.gene_reaction_rule + ''
        for operator in ['(', ')', ' and ', ' or ']:
            gpr = gpr.replace(operator, '  ')
        rxn_genes = set(gpr.split())
        if len(rxn_genes) > 0:
            rxn_to_genes[rxn.id] = rxn_genes

    genes_map = {}
    for geneID in genes_of_interest: # finally, get all reactions impacted by genes of interest
        if geneID != 's0001': # non-spontaneous
            gene = model.genes.get_by_id(geneID)
            match = gene.notes['match'] if 'match' in gene.notes else None
            affected_rxns = []
            for rxnID in rxn_to_genes:
                if geneID in rxn_to_genes[rxnID]:
                    affected_rxns.append(rxnID)
            genes_map[geneID] = (match, tuple(affected_rxns))
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