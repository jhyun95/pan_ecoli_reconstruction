#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 18:12:43 2019

@author: jhyun95

Code for downloading and constructing the E. coli reference pan-genome
from the genomes and models described in doi:10.1073/pnas.1307797110.
To setup from scratch, the functions should be run with default parameters.

download_reference_data()
create_label_table() 
create_non_redundant_pangenome()
create_cdhit_reference_database() # if using CD-Hit reconstruction
download_reference_data_dna() # if intending on analyzing upstream DNA
get_upstream_sequences() # if intending on analyzing upstream DNA

"""

from utils import process_header
import os, urllib, subprocess
import pandas as pd
import cobra

module_path = os.path.abspath(os.path.dirname(__file__)) + '/'
REFERENCE_LIST_PATH = module_path + 'data/bigg_ecoli_recons.csv' # path to list of reference models
FIX_REACTIONS_PATH = module_path + 'data/rxn_fixes.tsv' # path to metabolite definition overrides
FIX_METABOLITES_PATH = module_path + 'data/met_fixes.tsv' # path to reaction definition overrides
FIX_METABOLITES_DB_PATH = module_path + 'data/pan-ecoli-met.tsv' # path to E. cli pan-metabolome

WHITELIST_REACTIONS = ['HEPKB2', 'HEPKA2', 'GALR1TRA2'] # reactions to ignore when balancing 
WHITELIST_METABOLITES = ['colipaOA'] # metabolites to ignore when fixing metabolite formulas

def download_reference_data_dna(ref_dir='reference/', overwrite=True):
    '''
    Checks whether reference model DNA has been downloaded (specifically,
    the full nucleotide fasta). If missing, downloads from NCBI. Intended 
    to be used for extracting upstream DNA sequences of corresponding proteins. 

    Parameters
    ----------
    ref_dir : str 
        Directory to store genomes. Will create subdirectory "ref_genomes_dna/"
    overwrite : bool
        If true, re-downloads genomes already present (default False)
    '''

    df_rec = pd.read_csv(REFERENCE_LIST_PATH, index_col=0)

    ''' Check for/create output directories '''
    dna_dir = (ref_dir + '/ref_genomes_dna/').replace('//','/')
    if not os.path.isdir(ref_dir):
        os.mkdir(ref_dir)
    if not os.path.isdir(dna_dir):
        os.mkdir(dna_dir)

    ''' Downloading genome DNA and feature tables '''
    for strain in df_rec.index:
        ncbi = df_rec.loc[strain,'NCBI Full']
        ncbi_dna_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?'
        ncbi_dna_url += 'db=nucleotide&id=' + ncbi + '&rettype=fasta'
        ncbi_dna_out = dna_dir + ncbi + '.fna'
        if not os.path.exists(ncbi_dna_out) or overwrite:
            print 'Downloading full genome DNA', ncbi, 'from', ncbi_dna_url
            urllib.urlretrieve(ncbi_dna_url, filename=ncbi_dna_out)
        else: 
            print 'Full genome DNA', ncbi, 'already downloaded'


def download_reference_data(ref_dir='reference/', overwrite=False, repair_downloaded_models=True):
    ''' 
    Checks whether reference models and their corresponding genomes
    have been downloaded. If missing, downloads from BiGG and NCBI.
    Optionally, implements curated fixes to downloaded models. 
    Assumes models that are already present do not need to be repaired. 

    Parameters
    ----------
    ref_dir : str 
        Directory to store genomes and models. Will create subdirectories 
        "ref_genomes/" and "ref_models/" in this directory. 
    overwrite : bool
        If true, re-downloads models and genomes already present (default False)
    repair_downloaded_models : bool
        If true, repairs newly downloaded models with fix_model_balances() (default True) 
    '''

    df_rec = pd.read_csv(REFERENCE_LIST_PATH, index_col=0)

    ''' Check for/create output directories '''
    genomes_dir = (ref_dir + '/ref_genomes/').replace('//','/')
    models_dir = (ref_dir + '/ref_models/').replace('//','/')
    if not os.path.isdir(ref_dir):
        os.mkdir(ref_dir)
    if not os.path.isdir(genomes_dir):
        os.mkdir(genomes_dir)
    if not os.path.isdir(models_dir):
        os.mkdir(models_dir)

    ''' Downloading models and genomes '''
    downloaded_models = [] # models that have been freshly downloaded
    for strain in df_rec.index:
        ncbi = df_rec.loc[strain,'NCBI Full']
        bigg = df_rec.loc[strain,'Model Name']
        if bigg == 'iJO1366': # update the K12 MG1655 model
            bigg = 'iML1515'

        if bigg != 'iY75_1357': # don't use iY75_1357, no GPRs = cannot call homologs
            ncbi_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?'
            ncbi_url += 'db=nucleotide&id=' + ncbi + '&rettype=fasta_cds_aa' # get fasta of protein sequences
            ncbi_out = genomes_dir + ncbi + '.faa'
            bigg_url = 'http://bigg.ucsd.edu/static/models/' + bigg + '.json'
            bigg_out = models_dir + bigg + '.json'

            ''' Downloading NCBI protein features as FAA '''
            if not os.path.exists(ncbi_out) or overwrite:
                print 'Downloading genome', ncbi, 'from', ncbi_url
                urllib.urlretrieve(ncbi_url, filename=ncbi_out)
            else: 
                print 'Genome', ncbi, 'already downloaded'

            ''' Downloading BiGG models as json '''
            if not os.path.exists(bigg_out) or overwrite:
                print 'Downloading model', bigg, 'from', bigg_url
                urllib.urlretrieve(bigg_url, filename=bigg_out)
                downloaded_models.append(bigg_url)
            else: 
                print 'Model', bigg, 'already downloaded'

    ''' Optionally applying curated fixes '''
    if repair_downloaded_models:
        iML1515 = cobra.io.load_json_model(models_dir + 'iML1515.json')
        for model_path in downloaded_models:
            if not 'iML1515' in model_path: # don't correct iML1515
                ''' Fix metabolite formulas and unbalanced reactions '''
                model = cobra.io.load_json_model(model_path)
                model = fix_model_balances(model, clean_model=iML1515)

                ''' Fix proton symport to be irreversible '''
                for reaction in model.reactions:
                    if 'via proton symport (periplasm)' in reaction.name:
                        reaction.lower_bound = 0.0 
                cobra.io.save_json_model(model, model_path)

    ''' Merge into a non-redundant sequence database '''
    create_non_redundant_pangenome(ref_dir)


def get_upstream_sequences(ref_dir='reference/', limits=(-50,3)):
    ''' 
    Gets non-coding upstream DNA sequences for all protein coding features 
    for the reference sequences.

    Parameters
    ----------
    ref_dir : str 
        Directory to look for protein features and full genome DNA. Will
        also create subdirectory "ref_upstream" with upstream sequences 
    limits : tuple
        Length of upstream region to extract, formatted (-X,Y). Will extract X 
        upstream bases (up to but excluding first base of start codon) and Y coding 
        bases (including first base of start codon), for total length of X+Y bases.
        (default (-50,3))
    '''

    genomes_dir = (ref_dir + '/ref_genomes/').replace('//','/') # protein sequences, .faa files
    dna_dir = (ref_dir + '/ref_genomes_dna/').replace('//','/') # full genomic DNA, .fna files
    upstream_dir = (ref_dir + '/ref_upstream/').replace('//','/') # protein DNA upstream, fna files
    if not os.path.isdir(upstream_dir):
        os.mkdir(upstream_dir)

    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 
                  'W': 'W', 'S': 'S', 'R': 'Y', 'Y': 'R', 
                  'M': 'K', 'K': 'M', 'N': 'N'}
    reverse_complement = lambda s: (''.join([complement[base] for base in list(s)]))[::-1]
    upstream_label = 'upstream' + str(limits).replace(' ','')
    for genome_file in os.listdir(genomes_dir):
        strain = genome_file[:-4]
        dna_path = dna_dir + strain + '.fna'
        genome_path = genomes_dir + genome_file
        if genome_file[-4:] == '.faa' and os.path.exists(dna_path): # if faa/fna pair found
            print 'Extracting upstream sequences for', strain
            ''' Load genome DNA, assuming the NCBI genomes are a single contig '''
            genome = ''
            with open(dna_path, 'r') as f_dna:
                for line in f_dna:
                    if not '>' in line:
                        genome += line.strip()

            ''' Extract upstream sequence for each protein sequence '''
            upstream_out = upstream_dir + strain + '_upstream.fna'
            with open(genome_path, 'r') as f_in:
                with open(upstream_out, 'w+') as f_out:
                    for line in f_in:
                        if line[0] == '>': # header line
                            raw_name, parameters, label = process_header(line)
                            full_label = '>' + label + '|' + upstream_label
                            pseudo = 'pseudo' in parameters and parameters['pseudo'] == 'true'
                            if 'location' in parameters and not pseudo:
                                location = parameters['location']
                                if 'join' in location: # rare, frameshift/slippage case
                                    print '\tJOIN:', location, '->', 
                                    location = location.replace('join(','').replace('))',')')
                                    location = location.replace(',','..')
                                    tmp = location.split('..')
                                    location = '..'.join((tmp[0], tmp[-1]))
                                    location = location.replace(')','') if not 'complement' in location else location
                                    print location
                                if '>' in location or '<' in location: # rare, undefined CDS limits
                                    print '\tWARN (ignoring >/<):', raw_name, location
                                    location = location.replace('>','').replace('<','')
                                if 'complement' in location: # negative strand
                                    location = location.replace('complement(','').replace(')','')
                                    start, end = map(int, location.split('..'))
                                    upstream = genome[end-limits[1]:end-limits[0]]
                                    upstream = reverse_complement(upstream)
                                else: # positive strand
                                    start, end = map(int, location.split('..'))
                                    upstream = genome[start+limits[0]-1:start+limits[1]-1]
                                f_out.write(full_label + '\n')
                                f_out.write(upstream + '\n')


def create_label_table(ref_dir='reference/'):
    ''' Creates a table mapping protein names to labels of the form
        <strain>|<locus_tag>|<protein_id> based on processing all
        NCBI headers of the reference genomes '''
    genomes_dir = (ref_dir + '/ref_genomes/').replace('//','/')
    label_file = (ref_dir + '/ref_labels.tsv').replace('//','/')
    with open(label_file, 'w+') as f_label:
        for genome_file in os.listdir(genomes_dir):
            genome_path = genomes_dir + genome_file
            with open(genome_path, 'r') as f:
                for line in f:
                    if line[0] == '>': # if header
                        raw_name, parameters, label = process_header(line)
                        f_label.write(raw_name + '\t' + label + '\n')


def create_cdhit_reference_database(ref_dir='reference/', recreate=True,
    extra_cdhit_args=['-c', 0.8, '-aL', 0.8, '-T', 10, '-d', 0]):
    '''
    Creates a CD-Hit cluster database from the reference sequences. Requires
    pan-ecoli.faa in ref_dir generated by create_non_redundant_pangenome() with
    delete_redundant=False.

    Parameters
    ----------
    ref_dir : str 
        Directory with pan-ecoli.faa generated by create_non_redundant_pangenome().
        CD-Hit database will also be generated in this folder.
    recreate : bool
        If true, overwrites/recreates from given CD-Hit parameters (default True)
    extra_cdhit_args : list
        List of extra command line arguments to pass to CD-Hit. Default set to 80%
        identity, 90% coverage, 10 threads, print full names 
        (default ['-c', 0.8, '-aL', 0.8, -T', 10, '-d', 0])
    '''

    ''' Prepare directories/file paths'''
    cdhit_dir = (ref_dir + '/cdhit/').replace('//','/')
    cdhit_fasta = (ref_dir + '/cdhit/pan-ecoli-cdhit.faa').replace('//','/')
    ref_raw_fasta = (ref_dir + '/pan-ecoli.faa').replace('//','/')
    if not os.path.isdir(cdhit_dir):
        os.mkdir(cdhit_dir)

    ''' Create CD-Hit reference clusters '''
    cdhit_args = ['cd-hit', '-i', ref_raw_fasta, '-o', cdhit_fasta]
    cdhit_args += extra_cdhit_args
    print 'Running CD-Hit:', cdhit_args
    print subprocess.check_output(map(str, cdhit_args))


def create_non_redundant_pangenome(ref_dir='reference/', delete_redundant=False):
    ''' 
    Reduces downloaded reference genomes located in <ref_dir>/ref_genomes 
    into a non-redundant list (pan-ecoli-nr.faa), and stores mapping of 
    redundancies to a separate table (pan-ecoli-nr-map.tsv). 
    Note: Does not preserve any word-wrapping.

    Parameters
    ----------
    ref_dir : str 
        Directory with reference genomes as in download_reference_data()
    delete_redundant : bool
        If True, deletes the intermediate list of sequences with redundancies
        (pan-ecoli.faa). Must be false for CD-Hit based clustering. (default False)
    '''

    ''' Prepare directories/file paths '''
    genomes_dir = (ref_dir + '/ref_genomes/').replace('//','/')
    redundant_file = (ref_dir + '/pan-ecoli.faa').replace('//','/')
    nonredundant_file = (ref_dir + '/pan-ecoli-nr.faa').replace('//','/')
    locus_table_out = (ref_dir + '/pan-ecoli-nr-map.tsv').replace('//','/')

    ''' Concatenate genomes '''
    #merge_args = ['cat', genomes_dir+'*', '|', 'tr', '-s', '"\n"', '>', redundant_file] # one line command
    with open(redundant_file, 'w+') as fcat:
        for genome_file in os.listdir(genomes_dir):
            genome_path = genomes_dir + genome_file
            with open(genome_path, 'r') as f:
                for line in f:
                    fcat.write(line)

    ''' Record and remove redundancies '''
    unique_seqs = {} # maps sequences to stored headers
    locus_map = {} # maps stored headers to labels <strain>|<locus_tag>|<protein_id>
    entry_order = [] # order that sequences were stored
    counts = [0,0] # redundant/total

    def __store_sequence__(seq_header, sequence):
        ''' Save the current header/sequence pair, assumes no whitespace at ends '''
        raw_name, parameters, label = process_header(seq_header)
        if sequence in unique_seqs: # redundant sequence
            counts[0] += 1
            locus_map[unique_seqs[sequence]] += ';' + label
        else: # non-redundant sequence
            unique_seqs[sequence] = raw_name
            locus_map[raw_name] = label
            entry_order.append(raw_name)
            f_nonredundant.write(seq_header + '\n')
            f_nonredundant.write(sequence + '\n')
        counts[1] += 1

    ''' Iterate through all sequences in raw concatenated file '''
    current_seq = ''; header = ''; 
    with open(nonredundant_file, 'w+') as f_nonredundant:
        with open(redundant_file, 'r+') as f_redundant:
            for line in f_redundant:
                if '>' == line[0]: # header line
                    if len(current_seq) > 0 and len(header) > 0: # not first header / completed last sequence
                        __store_sequence__(header, current_seq)
                    header = line.strip() # get the new header
                    current_seq = '' # reset the current sequence
                else: # sequence line
                    current_seq += line.strip() # record sequence while removing word-wrap

            ''' Last header/sequence pair '''
            if len(current_seq) > 0 and len(header) > 0: # same checks for last sequence
                __store_sequence__(header, current_seq)

    ''' Remove the raw concatenated genomes file '''
    if delete_redundant:
        os.remove(redundant_file)

    ''' Report sequence redundancy and save redundancy table '''
    print 'Sequences:', counts[1]
    print 'Unique:', len(unique_seqs)
    print 'Redundant:', counts[0]
    df_locus = pd.DataFrame.from_dict(data=locus_map, orient='index', columns=['locus_tags'])
    df_locus = df_locus.reindex(entry_order, axis='index')
    df_locus.to_csv(locus_table_out, sep='\t')


def check_sequence_mapping(ref_dir='reference/'):
    '''
    Verifies that sequences were mapped correctly by create_non_redundant_pangenome()
    '''

    ''' Prepare directories/file paths '''
    genomes_dir = (ref_dir + '/ref_genomes/').replace('//','/')
    nonredundant_file = (ref_dir + '/pan-ecoli-nr.faa').replace('//','/')
    locus_table_out = (ref_dir + '/pan-ecoli-nr-map.tsv').replace('//','/')

    def __load_sequences__(seq_file):
        ''' Load sequences, does not preserve word wrap '''
        sequences = {} # maps headers to sequences
        header = ''; seq = ''
        with open(seq_file, 'r') as f:
            for line in f:
                if line[0] == '>': # header
                    if len(header) > 0 and len(seq) > 0:
                        sequences[header] = seq
                    header = line.strip(); seq = ''
                else: # sequence line
                    seq += line.strip()
        return sequences

    ''' Load non-redundant sequences and compare against original sequences '''
    discrepancies = 0
    nr_sequences = __load_sequences__(nonredundant_file)
    for genome_file in os.listdir(genomes_dir):
        genome_path = genomes_dir + genome_file
        original_sequences = __load_sequences__(genome_path)
        for header in original_sequences:
            if header in nr_sequences:
                original_seq = original_sequences[header]
                nr_seq = nr_sequences[header]
                if original_seq != nr_seq:
                    discrepancies += 1
                    print 'DISCREPANCY:', header
    print 'Discrepancies:', discrepancies


def check_model_balances(ref_dir='reference/', force_repair=False):
    ''' 
    Report any charge or mass imbalances in a model. Can optionally force
    running fix_model_balances() before verifying.

    Parameters
    ----------
    ref_dir : str 
        Directory with reference models as in download_reference_data()
    force_repair : bool
        If true, runs fix_model_balances() and overwrites existing model
        for all reference models in ref_dir (default false)
    '''
    models_dir = ref_dir + 'ref_models/'
    for model_file in os.listdir(models_dir):
        model_path = models_dir + model_file
        model = cobra.io.load_json_model(model_path)

        ''' Optionally force repair existing models before verifying '''
        if force_repair:
            iML1515 = cobra.io.load_json_model(models_dir + 'iML1515.json')
            if not 'iML1515' in model_path: # don't correct iML1515
                print 'REPAIRING MODEL', model_path
                ''' Fix metabolite formulas and unbalanced reactions '''
                model = cobra.io.load_json_model(model_path)
                model = fix_model_balances(model, clean_model=iML1515)
                ''' Fix proton symport to be irreversible '''
                for reaction in model.reactions:
                    if 'via proton symport (periplasm)' in reaction.name:
                        reaction.lower_bound = 0.0 
                cobra.io.save_json_model(model, model_path)
        
        print 'CHECKING MODEL', model_path
        unbalance = [0,0]
        ''' Check for metabolites without formulas or charges '''
        for met in model.metabolites:
            if met.formula is None or met.charge is None or len(met.formula) == 0:
                print 'Undefined:', met.id, met.formula, met.charge,
                if 'colipaOA' in met.id or 'LptA' in met.id:
                    print '(whitelisted)'
                else:
                    print ''
                unbalance[0] += 1
        print 'Undefined metabolites:', unbalance[0]
        
        ''' Check for imbalanced reactions '''
        for rxn in model.reactions:
            if not 'EX_' in rxn.id and not 'DM_' in rxn.id and \
                not 'SK_' in rxn.id and not 'BIOMASS_Ec' in rxn.id:
                test = rxn.check_mass_balance()
                if len(test) > 0:
                    if rxn.id in ['HEPKB2', 'HEPKA2', 'GALR1TRA2']:
                        print 'Unbalanced:', rxn.id, '(whitelisted)'
                    else:
                        print 'Unbalanced:', rxn, test
                    unbalance[1] += 1
        print 'Unbalanced reactions:', unbalance[1]
 

def fix_model_balances(model, clean_model=None):
    ''' 
    For the original 55 models from doi:10.1073/pnas.1307797110 as available on BiGG,
    attempts to fix undefined metabolites and unbalanced reactions by:
    1) Implement manually curated fixes to metabolite formulas 
    2) Make remaining undefined metabolite formulas consistent with a clean model (i.e. iML1515)
    3) Make remaining undefined metabolite formulas consistent with the SD01 reference table
    4) Implement manually curated fixes to reactions 
    5) Add hydrogens if imbalance is only due to H and charge 

    Parameters
    ----------
    model : cobra Model
        Reference cobra model to be balanced/fixed
    clean_model : str or cobra Model
        Model or path to Model to use as a reference, such as iML1515 (default None)

    Returns
    -------
    model : cobra Model
        Balanced Model (modified in place)
    '''

    df_mets = pd.read_csv(FIX_METABOLITES_DB_PATH, sep='\t', index_col=0)
    df_fixes_met = pd.read_csv(FIX_METABOLITES_PATH, sep='\t', index_col=0)
    df_fixes_rxn = pd.read_csv(FIX_REACTIONS_PATH, sep='\t', index_col=0)
    imbalanced_reactions = 0
    undefined_metabolites = 0

    def __update_metabolite__(met, new_formula, new_charge, label):
        updated = False
        if (not new_formula is None) and (not pd.isnull(new_formula)) and (met.formula != new_formula):
            print met, label, 'Formula:', met.formula, '->', new_formula
            met.formula = new_formula
            updated = True
        if (not new_charge is None) and (not pd.isnull(new_charge)):
            if met.charge is None or int(met.charge) != int(new_charge):
                print met, label, 'Charge:', met.charge, '->', int(new_charge)
                met.charge = int(new_charge)
                updated = True
        return updated

    ''' Step 1: Hard-coded fixes for metabolites from manual curation '''
    print '---------------- Fixing metabolite formulas ----------------'
    resolved_mets = set(df_fixes_met.index.tolist()) # track corrected metabolites to not adjust multiple times
    for met_name in df_fixes_met.index:
        for compartment in model.compartments:
            metID = met_name + '_' + str(compartment)
            if metID in model.metabolites:
                met = model.metabolites.get_by_id(metID)
                ref_formula = df_fixes_met.loc[met_name, 'Formula']
                ref_charge = df_fixes_met.loc[met_name, 'Charge']
                updated = __update_metabolite__(met, ref_formula, ref_charge, 'Curated')

    ''' Step 2/3: Update formulas/charges to match the clean model, or if not in
            the clean model, match the reference table '''
    for met in model.metabolites:
        met_name = met.id[:-2]
        if not met_name in WHITELIST_METABOLITES and not met_name in resolved_mets:
            if met.id in clean_model.metabolites: # check clean model first
                clean_met = clean_model.metabolites.get_by_id(met.id)
                updated = __update_metabolite__(met, clean_met.formula, clean_met.charge, 'Clean')
            else: # then check the SD01 reference table, df_mets
                met_index1 = met.id.replace('__','-')
                met_index2 = met.id[:-2].replace('_','-') + met_index1[-2:]
                possible_matches = [met.id, met_index1, met_index2]; hit = None
                for match in possible_matches:
                    if match in df_mets.index:
                        hit = match; break
                if hit: # found a match in the reference table
                    ref_formula = df_mets.loc[hit, 'Formula']
                    ref_charge = df_mets.loc[hit, 'Charge']
                    if pd.isnull(ref_formula): # no formula value, infer from name
                        # TODO: Safer checking string is formula
                        ref_name = df_mets.loc[hit, 'Name']
                        infer_formula = ref_name.split('_')[-1]
                        if infer_formula.upper() == infer_formula:
                            ref_formula = infer_formula
                    if ref_charge == 0 and not met.charge is None:
                        # Do not update charge to 0 if there is a pre-existing charge assigned
                        ref_charge = met.charge    
                    updated = __update_metabolite__(met, ref_formula, ref_charge, 'Ref')

            if met.formula is None or met.charge is None:
                undefined_metabolites += 1
                print 'UNDEFINED', undefined_metabolites, met, met.formula, met.charge
            elif updated:
                resolved_mets.add(met.id)

    print '\n---------------- Fixing reaction imbalances ----------------'
    ''' Step 4: Hard-coded fixes for reactions from manual curation '''
    for rxnID in df_fixes_rxn.index:
        if rxnID in model.reactions:
            rxn = model.reactions.get_by_id(rxnID)
            ref_rxn_str = df_fixes_rxn.loc[rxn.id, 'Formula']
            if ref_rxn_str == 'DELETE':
                rxn.remove_from_model()
                print '\t', rxnID, 'Curated Reaction: DELETED'
            else:
                rxn.build_reaction_from_string(ref_rxn_str)
                print '\t', rxnID, 'Curated Reaction:', rxn.build_reaction_string()

    ''' Step 5: Add hydrogens if only imbalances are H and charge
        and the reaction is entirely localized to a single compartment '''
    for rxn in model.reactions:
        need_balance = rxn.id[:3] != 'EX_' and rxn.id[:3] != 'DM_' and rxn.id[:3] != 'SK_'
        need_balance = need_balance and not 'BIOMASS' in rxn.id and not rxn.id in WHITELIST_REACTIONS
        balance = rxn.check_mass_balance()
        if need_balance and len(balance) > 0:
            print rxn
            print '\tInitial balance:', balance 
            if len(balance) == 2 and 'H' in balance and 'charge' in balance:
                if balance['H'] == balance['charge']:
                    met_compartments = set(map(lambda m: m.id[-2:], rxn.metabolites))
                    if len(met_compartments) == 1: # single compartment
                        compartment = list(met_compartments)[0]
                        h_met = model.metabolites.get_by_id('h' + compartment)
                        rxn.add_metabolites({h_met:-balance['H']})
                        balance = rxn.check_mass_balance()
                        print '\tAdjusting hydrogens:', rxn.build_reaction_string()

            print '\tFinal Balance:', balance
            if len(balance) > 0:
                imbalanced_reactions += 1
                print '\tUNRESOLVED', imbalanced_reactions, rxn.gene_reaction_rule

    return model