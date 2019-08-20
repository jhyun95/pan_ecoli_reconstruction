#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 18:12:43 2019

@author: jhyun95
"""

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

def download_reference_data(ref_dir='reference/', overwrite=False, repair_downloaded_models=True):
    ''' 
    Checks whether reference models and their corresponding genomes
    have been downloaded. If missing, downloads from BiGG and NCBI.
    Optionally, implements curated fixes to downloaded models. 
    Assumes models that are already present do not need to be repaired. 

    Parameters
    ----------
    ref_dir : str 
        Directory to store genomes and models. Will create two subdirectories 
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

            if not os.path.exists(ncbi_out) or overwrite:
                print 'Downloading genome', ncbi, 'from', ncbi_url
                urllib.urlretrieve(ncbi_url, filename=ncbi_out)
            else: 
                print 'Genome', ncbi, 'already downloaded'
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


def create_non_redundant_pangenome(ref_dir='reference/'):
    ''' 
    Reduces downloaded reference genomes into a non-redundant list,
    and stores mapping of redundancies to a separate table. 
    Note: Does not preserve any word-wrapping.
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
        ''' Save the current header/sequence pair, assumes no whitespace at ends
        1) Extract information from header for storing redundancies '''
        raw_name = seq_header.split()[0][1:]
        strain_name = seq_header.split('_prot_')[0][5:]
        parameters = {}
        for annotation in seq_header.split()[1:]:
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

        ''' 2) Check redundancy to write to file, save labels '''
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
    os.remove(redundant_file)

    ''' Report sequence redundancy and save redundancy table '''
    print 'Sequences:', counts[1]
    print 'Unique:', len(unique_seqs)
    print 'Redundant:', counts[0]
    df_locus = pd.DataFrame.from_dict(data=locus_map, orient='index', columns=['locus_tags'])
    df_locus = df_locus.reindex(entry_order, axis='index')
    df_locus.to_csv(locus_table_out, sep='\t')
    return df_locus


def check_model_balances(ref_dir='reference/', force_repair=False):
    ''' 
    Report any charge or mass imbalances in a model.
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