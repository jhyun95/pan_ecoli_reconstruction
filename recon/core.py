#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 14:36:29 2019

@author: jhyun95
"""

import os, subprocess, simplejson, time
import numpy as np
import numexpr as ne
import pandas as pd
import cobra

BLACKLIST_PROTEINS = ['AYC08180.1','YP_009518752.1'] # ignore alignments to these proteins
BASELINE_GENES = ['s0001'] # genes to include in all models
BASELINE_STRAIN = ('NC000913.3','iML1515') # start with this genome/model during reconstruction
NONENZYMATIC_TERMS = ('BIOMASS', 'ATPM', 'EX_', 'DM_', 'SK_')


def reconstruct_with_cdhit(seq_fasta, work_dir=None, ref_dir='reference/', 
    ref_cdhit='cdhit/', extra_cdhit_args=['-c', 0.8, '-aL', 0.8, '-T', 10, '-d', 0],
    accept_gprless=['DNMPPA']):
    '''
    Attempts a reconstruction with CD-Hit given a set of coding sequences 
    '''

    ''' Timing benchmarks '''
    time_cdhit_run = 0 # time for CD-Hit run 
    time_cdhit_parse = 0 # time for CD-Hit co-cluster parsing
    time_base = 0 # time for base reference model copying/initialization
    time_expand = 0 # time for integration of other reference models

    ''' Prepare file paths and output directory '''
    time_cdhit_run = time.time()
    label_ref_file = (ref_dir + '/ref_labels.tsv').replace('//','/')
    cdhit_ref_dir = (ref_dir + '/' + ref_cdhit).replace('//','/')
    cdhit_ref_db = cdhit_ref_dir + 'pan-ecoli-cdhit.faa'
    head, tail = os.path.split(seq_fasta)
    query_name = '.'.join(tail.split('.')[:-1])
    if work_dir is None:    
        work_dir = query_name + '_recon/'
    if not os.path.isdir(work_dir):
        os.mkdir(work_dir)
    work_dir = (work_dir + '/').replace('//','/')
    cdhit_temp = work_dir + query_name + '_cdhit_temp.faa'
    cdhit_final_out = work_dir + query_name + '_cdhit_merged.faa'

    ''' Run CD-Hit incremental clustering against reference clusters, adapted from
        https://github.com/weizhongli/cdhit/wiki/3.-User's-Guide#CDHIT
        (skipping some steps, don't need preserve clusters not involving reference genes) '''

    # 1) Incremental clustering against reference -> cdhit_temp
    cdhit_args = ['cd-hit-2d', '-i', cdhit_ref_db, '-i2', seq_fasta, '-o', cdhit_temp]
    cdhit_args += map(str,extra_cdhit_args)
    print 'RUNNING:', cdhit_args
    print subprocess.check_output(cdhit_args)

    # 2) Re-clustering -> cdhit_out2 (SKIP: Don't need non-co-clustered genes)
    # cdhit_args = ['cd-hit', '-i', cdhit_out1, '-o', cdhit_out2]
    # cdhit_args += map(str,extra_cdhit_args)
    # print 'RUNNING:', cdhit_args
    # print subprocess.check_output(cdhit_args)

    # 3) Merge clusters -> cdhit_out3 / cdhit_final_out
    cdhit_args = ['clstr_merge.pl', cdhit_ref_db + '.clstr', cdhit_temp + '.clstr']
    print 'RUNNING:', cdhit_args
    with open(cdhit_final_out + '.clstr', 'w+') as f:
        subprocess.call(cdhit_args, stdout=f)

    # 4) Build final cluster set -> cdhit_final_out (SKIP: Don't need non-co-clustered genes)
    # cdhit_args = ['cat', cdhit_out3 + '.clstr', cdhit_out2 + '.clstr']
    # print 'RUNNING:', cdhit_args
    # with open(cdhit_final_out + '.clstr', 'w+') as f:
    #     subprocess.call(cdhit_args, stdout=f) 

    # 5) Delete temporary files
    subprocess.call(['rm', cdhit_temp])
    subprocess.call(['rm', cdhit_temp + '.clstr'])
    time_cdhit_run = time.time() - time_cdhit_run

    ''' Load label table '''
    time_cdhit_parse = time.time()
    header_to_locus = {}
    with open(label_ref_file, 'r') as f:
        for line in f:
            name, label = line.split()
            header_to_locus[name] = label.split('|')[1]

    ''' Identify co-clustered reference genes '''
    #query_to_reference = {} # match query to sets of labels for co-clustered reference genes
    reference_to_query = {} # match reference gene to co-clustered nearest query gene
    current_cluster = set(); current_query = ''; current_pid = 0.0
    with open(cdhit_final_out + '.clstr', 'r') as f:
        for line in f:
            if line[:8] == '>Cluster': # starting a new cluster
                if len(current_query) > 0: # cluster contains a query sequence
                    for ref_gene in current_cluster:
                        reference_to_query[ref_gene] = current_query
                    # query_to_reference[current_query] = current_cluster
                current_cluster = set(); current_query = '' # reset cluster
            else: # continuing existing cluster 
                entries = line.split()
                length = int(entries[1].split('aa')[0])
                name = entries[2][1:-3] # drop '>' and '...'
                if name in header_to_locus: # if reference gene
                    current_cluster.add(header_to_locus[name])
                else: # if query gene
                    pid = float(entries[-1][:-1]) # perecent identity to cluster representative sequence
                    if current_query != '': # if multiple queries in same cluster, choose higher PID
                        if current_pid < pid: # if new hit has higher percent, overwrite
                            current_query = name; current_pid = pid
                    else: # first query gene in this cluster
                        current_query = name; current_pid = pid
    print '------------------------------------------------------------\n'
    print '# of query sequences co-clustered with reference genes:', len(set(reference_to_query.values()))
    time_cdhit_parse = time.time() - time_cdhit_parse

    ''' Initialize the model as a subset of the base strain '''
    time_base = time.time()
    models_dir = (ref_dir + '/ref_models/').replace('//','/')
    ref_model_index = 1
    total_ref_models = len(filter(lambda f: f[-5:] == '.json', os.listdir(models_dir)))
    label = str(ref_model_index) + '/' + str(total_ref_models)

    print '\n' + label + ': Initializing with hits to', BASELINE_STRAIN[1]
    base_model = models_dir + BASELINE_STRAIN[1] + '.json'
    model = cobra.io.load_json_model(base_model)
    model.name = query_name
    base_model_geneIDs = map(lambda g: g.id, model.genes)
    missing_genes = filter(lambda g: not g in reference_to_query, base_model_geneIDs)
    for gene in BASELINE_GENES: # don't delete baseline genes
        missing_genes.remove(gene)
    cobra.manipulation.delete.remove_genes(model, missing_genes) # hard-delete missing genes
    for gene in model.genes: # add note about matched gene
        if gene.id in reference_to_query:
            gene.notes['match'] = reference_to_query[gene.id]
    gprless_enzymatic = [] # enzymatic reactions without GPRs to remove
    for reaction in model.reactions:
        if len(reaction.gene_reaction_rule) == 0:
            rxnID = reaction.id
            is_nonenzymatic = map(lambda x: x in rxnID, NONENZYMATIC_TERMS)
            is_nonenzymatic = reduce(lambda x,y: x or y, is_nonenzymatic)
            if (not is_nonenzymatic) and (not rxnID in accept_gprless): # if enzymatic but not GPR, remove
                gprless_enzymatic.append(rxnID)
    if len(gprless_enzymatic) > 0:
        model.remove_reactions(gprless_enzymatic)
    print '\tGenes:', len(model.genes)
    print '\tReactions:', len(model.reactions)
    print '\tMetabolites:', len(model.metabolites)
    time_base = time.time() - time_base

    ''' Update model with hits to other models '''
    time_expand = time.time()
    for model_file in os.listdir(models_dir):
        ''' For each model that is not the baseline strain '''
        if model_file[-5:] == '.json' and not BASELINE_STRAIN[1] in model_file:
            ''' Load with simplejson not cobrapy for better performance '''
            model_name = model_file[:-5]
            model_path = models_dir + model_file
            ref_model_index += 1
            label = str(ref_model_index) + '/' + str(total_ref_models)
            print label + ': Updating with hits to', model_name
            with open(model_path, 'r') as model_json_file:
                model_json = simplejson.load(model_json_file)

            ''' Build gene presence/absence table '''
            gene_presence = {} # maps locus tags to True/False based on presence in query 
            for gene in model_json['genes']:
                gene_presence[gene['id']] = gene['id'] in reference_to_query
            for gene in BASELINE_GENES:
                gene_presence[gene] = True

            ''' Map gene/metabolite information to IDs '''
            model_gene_dumps = {}; model_met_dumps = {}
            for gene in model_json['genes']:
                model_gene_dumps[gene['id']] = gene
            for met in model_json['metabolites']:
                model_met_dumps[met['id']] = met

            ''' Check reactions not in model with satisfied GPRs '''
            possible_reactions = 0; added_reactions = 0
            defunct_genes = [] # genes in a GPR of a new reaction but not present in strain
            for reaction in model_json['reactions']:
                not_present = not (reaction['id'] in model.reactions) # only new reactions
                not_biomass = not ('BIOMASS' in reaction['id']) # do not overwrite iML1515 biomass
                if not_present and not_biomass: 
                    possible_reactions += 1
                    gpr = reaction['gene_reaction_rule']
                    if len(gpr) == 0: # if no GPR, only allow non-enzymatic reactions
                        if accept_gprless == True: # special option for taking all gprless
                            gpr_state = True
                        else: # otherwise, only allow non-enzymatic reactions without GPRs
                            rxnID = reaction['id']
                            is_nonenzymatic = map(lambda x: x in rxnID, NONENZYMATIC_TERMS)
                            is_nonenzymatic = reduce(lambda x,y: x or y, is_nonenzymatic)
                            gpr_state = is_nonenzymatic or (rxnID in accept_gprless)
                    else: # yes GPR, attempt to evaluate statement
                        ''' Format GPR for numexpr and extract relevant genes '''
                        gpr_eval = gpr.replace(' or ',' | ').replace(' and ',' & ')
                        gpr_genes = gpr_eval + ''
                        for op in ['(', ')', '|', '&']:
                            gpr_genes = gpr_eval.replace(op ,' ')
                        gpr_genes = gpr_genes.split()

                        ''' If GPR contains any unknown genes (i.e. not mapped to locus tags)
                            assume they are not present in the strain '''
                        for gpr_gene in gpr_genes:
                            if not gpr_gene in gene_presence:
                                gene_presence[gpr_gene] = False

                        ''' Evaluate the GPR statement '''
                        try:
                            gpr_state = ne.evaluate(gpr_eval, local_dict=gene_presence)
                        except KeyError:
                            print 'WARNING: Could not resolve GPR', gpr_eval, gpr_genes
                    added_reactions += int(gpr_state)

                    ''' Add a present reaction from raw json data '''
                    if gpr_state: # if reaction is viable, add necessary genes, metabolites, and reaction
                        ''' Add missing metabolites from new reaction '''
                        new_metabolites = []
                        for metID in reaction['metabolites']:
                            met = model_met_dumps[metID]
                            if not met['id'] in model.metabolites:
                                new_met = cobra.core.Metabolite(
                                    id=met['id'],
                                    formula=met['formula'] if 'formula' in met else None,
                                    name=met['name'] if 'name' in met else '',
                                    charge=int(met['charge']) if 'charge' in met else None,
                                    compartment=met['compartment'] if 'compartment' in met else None)
                                new_met.notes = met['notes']
                                new_metabolites.append(new_met)
                        if len(new_metabolites) > 0:
                            model.add_metabolites(new_metabolites)

                        ''' Add the new reaction to the model '''
                        new_rxn = cobra.core.Reaction(
                            id=reaction['id'],
                            name=reaction['name'] if 'name' in reaction else '',
                            subsystem=reaction['subsystem'] if 'subsystem' in reaction else '',
                            lower_bound=float(reaction['lower_bound']) if 'lower_bound' in reaction else -1000.0,
                            upper_bound=float(reaction['upper_bound']) if 'upper_bound' in reaction else 1000.0)
                        new_rxn.notes = reaction['notes']
                        new_rxn.gene_reaction_rule = reaction['gene_reaction_rule']
                        model.add_reactions([new_rxn])
                        new_rxn.add_metabolites(reaction['metabolites'])

                        ''' Annotate newly added genes '''
                        gpr_genes = gpr + ''
                        for c in ['(',')','or','and']: # remove GPR operators to get just genes
                            gpr_genes = gpr_genes.replace(c,' ')
                        gpr_genes = gpr_genes.split() # list of genes in GPR
                        for geneID in gpr_genes:
                            gene = model.genes.get_by_id(geneID)
                            gene_dump = model_gene_dumps[geneID]
                            gene.notes = gene_dump['notes']
                            gene.name = gene_dump['name'] if 'name' in gene_dump else ''
                            if geneID in reference_to_query: # gene is present in strain
                                gene.notes['match'] = reference_to_query[geneID]
                            elif not geneID in BASELINE_GENES: # gene is not present
                                defunct_genes.append(geneID)

            cobra.manipulation.delete.remove_genes(model, defunct_genes) # hard-delete missing geness
            print model_name + ': Added', added_reactions, 'of', possible_reactions,
            print 'possible new reactions.'
            ''' TODO: Try adding models and reactions all at once? '''

    time_expand = time.time() - time_expand
    print 'Final Model:'
    print '\tGenes:', len(model.genes)
    print '\tReactions:', len(model.reactions)
    print '\tMetabolites:', len(model.metabolites)

    print '\nRunning times:'
    print '\tCD-Hit run:', round(time_cdhit_run,3)
    print '\tCD-Hit parse:', round(time_cdhit_parse,3)
    print '\tCopying iML1515 hits:', round(time_base,3)
    print '\tIntegrating other hits:', round(time_expand,3)

    cobra.io.save_json_model(model, work_dir + query_name + '.json')
    return model


def reconstruct_with_blast(seq_fasta):
    ''' 
    Attempts a reconstruction with bi-directional blast given a set of coding sequences.
    '''
    pass

def bidirectional_blast(seq1_fasta, seq2_fasta, report_dir, 
        blast_params={'evalue':0.0000000001, 'num_threads':1}, 
        seq1_blast_db=None, seq2_blast_db=None, verbose=True):
    ''' 
    Takes FASTA files in the format <strain>.faa or <strain>.fna
    and applies a bidirectional blast against a set of reference sequences.

    Parameters
    ----------
    seq1_fasta : str
        Path to first set of sequences as FASTA
    seq2_fasta : str
        Path to second set of sequences as FASTA
    report_dir : str
        Directory to output blast files
    blast_params : dict
        Additional blastp parameters (default: {'evalue':0.0000000001, 'num_threads':1})
    seq1_blast_db : str
        Directory to check for/generate blast databases for first set. 
        If None, will use the same directory and seq1_fasta (default None)
    seq2_blast_db : str 
        Directory to check for/generate blast databases for second set. 
        If None, will use the same directory and seq2_fasta (default None)
    verbose : bool
        Whether or not to output blast command stdout (default True)

    Returns
    -------
    forward_report : str
        Path to TSV blast report for seq1 as query -> seq2 as database
    reverse_report : str
        Path to TSV blast report for seq2 as query -> seq1 as database
    '''

    dtype = 'prot' if seq1_fasta[-4:] == '.faa' else 'nucl'
    add_slash = lambda x: x if x[-1] =='/' else x + '/'
    get_file_dir = lambda x: '/'.join(x.split('/')[:-1])+'/' if '/' in x else ''
    get_filename = lambda x: x.split('/')[-1]
    name1 = '.'.join(get_filename(seq1_fasta).split('.')[:-1])
    name2 = '.'.join(get_filename(seq2_fasta).split('.')[:-1])

    ''' Prepare outputs and check for/generate blast databases '''
    report_dir = add_slash(report_dir) # make sure / is at the end of directory to output blast reports
    forward_report = report_dir + name1 + '_to_' + name2 + '.tsv'
    reverse_report = report_dir + name2 + '_to_' + name1 + '.tsv'
    db1_dir = add_slash(seq1_blast_db) if seq1_blast_db else get_file_dir(seq1_fasta) # directory with query blast files
    db2_dir = add_slash(seq2_blast_db) if seq2_blast_db else get_file_dir(seq2_fasta) # directory with ref blast files

    for directory in [report_dir, db1_dir, db2_dir]:
        if not os.path.isdir(directory):
            os.mkdir(directory)

    def __make_blast_db__(seq_path, db_out_dir=None):
        ''' Generates a blast database for a protein coding sequence fasta in 
            the same directory by default or specified directory. '''
        seq_filename = get_filename(seq_path)
        db_dir = add_slash(db_out_dir) if db_out_dir else get_file_dir(seq_path)
        db_name = db_dir + seq_filename
        db_files = map(lambda ext: db_name + ext, ['.phr', '.pin', '.psq'])
        db_exists = True
        for db_file in db_files: # check that the three blast db files <name>.phr/pin/psq exist
            db_exists = db_exists and os.path.exists(db_file)
        if not db_exists: # if at least one db file doesn't exist, regenerate db
            init_db_args = ['makeblastdb', '-in', seq_path, '-input_type', 'fasta', 
                            '-dbtype', dtype, '-out', db_name]
            __run_commands__(init_db_args, log_to_console=verbose)
        return db_name

    ''' Run bidirectional blast '''
    db1_name = __make_blast_db__(seq1_fasta, db_out_dir=db1_dir) # path to query blast db
    db2_name = __make_blast_db__(seq2_fasta, db_out_dir=db2_dir) # path to reference blast db
    extra_blast_args = []
    for param, value in blast_params.items():
        extra_blast_args += ['-' + str(param), str(value)]
    blasttype = 'blastp' if dtype == 'prot' else 'blastn'
    forward_args = [blasttype, '-db', db2_name, '-query', seq1_fasta,
                    '-out', forward_report, '-outfmt', '6'] + extra_blast_args
    reverse_args = [blasttype, '-db', db1_name, '-query', seq2_fasta,
                    '-out', reverse_report, '-outfmt', '6'] + extra_blast_args
    __run_commands__(forward_args, log_to_console=verbose)
    __run_commands__(reverse_args, log_to_console=verbose)

    return forward_report, reverse_report


def process_bidirectional_blast(blast1, blast2, output_file=None, identity_cutoff=70.0, 
        check_sort=False, verbose=True, deconvolution=None):
    '''
    Processes results from bi-diirection blast, extracting only pairwise top hits. 
    Can optionally limit to hits above an identity cutoff. Can optionally check if
    blast records are sorted. 

    Parameters
    ----------
    blast1 : str
        Path to 1st TSV blast report, i.e. forward_report from bidirectional_blast()
    blast2 : str
        Path to 2nd TSV blast report, i.e. reverse_report from bidirectional_blast()
    output_file : str
        Path to save top bidirectional hits if not None, must be TSV (default None)
    identity_cutoff : float
        Percent identity at which to accept an alignment. Default of 70.0 based on
        reconstruction workflow in Monk et al. 10.1073/pnas.1307797110 (default 70.0)
    check_sort : bool
        Whether or not to check if records are already sorted by bit score, which
        is usually the case. Runs much faster if False (default False)
    verbose : bool
        Whether or not to print the # of forward/reverse/bi-directional hits (default True)
    deconvolution : str
        Path to locus tag deconvolution table from remove_redundant_sequences(), i,e.
        blast1 blasts 1->2, and 2 was had merged redundant sequences (default None)
    '''

    def __get_relevant_hits__(df_blast, strain_name):
        ''' Extract records where strain is present in the hit name '''
        return df_blast[strain_name in df_blast[1]]

    def __get_best_hits_unsorted__(df_blast):
        ''' Get the best hit for each query by bit-score, slow but no assumptions on sorting '''
        best_hits = {}
        for query in df_blast[0].unique():
            df_hits = df_blast1[df_blast1[0]==query] # all hits for this query
            df_hits = df_hits.sort_values(by=11, ascending=False) # sort by bit-score
            if identity >= identity_cutoff:
                best_hits[query] = tuple(df_hits.iloc[0,:].to_list())
                #best_hits[query] = (best_hit, identity, evalue)
        return best_hits

    def __get_best_hits_presorted__(df_blast):
        ''' Get the best hit for each query by bit-score, slow but no assumptions on sorting '''
        best_hits = {}; query = ''
        for row in df_blast.itertuples(name=None):
            query = row[1]
            if not query in best_hits:
                best_hits[query] = row[2:]
        return best_hits

    def __get_bidirectional_best__(forward_best, reverse_best, map_hits=None):
        ''' Extract bi-directional best hits from single-direction best hits '''
        bidirectional_best = []
        for fwd_query in sorted(forward_best.keys()):
            fwd_hit = forward_best[fwd_query][0]
            if fwd_hit in reverse_best and reverse_best[fwd_hit][0] == fwd_query:
                fwd_entry = forward_best[fwd_query]
                rev_entry = reverse_best[fwd_hit]
                entry = list(fwd_entry) # export entry based on forward hit
                entry[1] = np.min([fwd_entry[1], rev_entry[1]]) # use smaller % identity
                entry[-2] = np.max([fwd_entry[-2], rev_entry[-2]]) # use larger evalue
                if map_hits:
                    #fwd_hit = map_hits(fwd_hit)
                    entry[0] = map_hits(entry[0])
                entry = tuple([fwd_query] + entry)
                bidirectional_best.append(entry)
        return bidirectional_best

    get_best_hits = __get_best_hits_unsorted__ if check_sort else __get_best_hits_presorted__
    df_blast1 = pd.read_csv(blast1, delimiter='\t', header=None)
    df_blast2 = pd.read_csv(blast2, delimiter='\t', header=None)
    df_blast1 = df_blast1[df_blast1[2] > identity_cutoff]
    df_blast2 = df_blast2[df_blast2[2] > identity_cutoff]

    if deconvolution: # pass per reference strain, as inferred from the deconvolution table
        ''' Identify list of reference strains for table '''
        if verbose:
            print 'Processing mapping of hits to strain|locus|protein...'
        strains = set()
        df_decon = pd.read_csv(deconvolution, delimiter='\t', index_col=0)
        decon_map = {}
        for raw_name, synonyms in df_decon.itertuples(name=None):
            decon_map[raw_name] = {}
            for entry in synonyms.split(';'):
                strain, locus, protein_id = entry.split('|')
                strains.add(strain)
                decon_map[raw_name][strain] = entry
                
        ''' Convert reference feature names to locus tag mappings '''
        #decon_map = df_decon.to_dict().values()[0]

        ''' Split reports by strain, and identify reference-specific best hits '''
        if verbose:
            print 'Identifying reference-specific best hits...'
        fwd_count = []; rev_count = []; bidir_count = []
        best_bidirectional_hits = []
        for strain in sorted(list(strains)):
            mask1 = df_blast1[1].apply(lambda x: strain in decon_map[x]) #.str.contains(strain) # in forward, check hits
            mask2 = df_blast2[0].apply(lambda x: strain in decon_map[x]) #.str.contains(strain) # in reverse, check queries
            df_strain1 = df_blast1[mask1]
            df_strain2 = df_blast2[mask2]
            #print df_blast1.shape, df_blast2.shape, '->', df_strain1.shape, df_strain2.shape
            forward_best_hits = get_best_hits(df_strain1)
            reverse_best_hits = get_best_hits(df_strain2)
            specific_bidirectional_hits = __get_bidirectional_best__(forward_best_hits, 
                reverse_best_hits, map_hits=lambda x: decon_map[x][strain])
            best_bidirectional_hits += specific_bidirectional_hits
            #print len(forward_best_hits), len(reverse_best_hits), len(specific_bidirectional_hits)
            fwd_count.append(len(forward_best_hits))
            rev_count.append(len(reverse_best_hits))
            bidir_count.append(len(specific_bidirectional_hits))

        if verbose:
            print 'Found', np.min(fwd_count), '~', np.max(fwd_count),
            print 'forward hits per reference above', str(identity_cutoff)+'%', 'cutoff.'
            print 'Found', np.min(rev_count), '~', np.max(rev_count),
            print 'reverse hits per reference above', str(identity_cutoff)+'%', 'cutoff.'
            print 'Found', np.min(bidir_count), '~', np.max(bidir_count), 
            print 'bi-directional best hits per reference.'

    else: # single pass, don't try to deconvolute different strains
        ''' Identify best hits by bit-score in each direction '''
        forward_best_hits = get_best_hits(df_blast1)
        reverse_best_hits = get_best_hits(df_blast2)
        if verbose:
            print 'Found', len(forward_best_hits), 'hits above', str(identity_cutoff)+'%', 'cutoff.'
            print 'Found', len(reverse_best_hits), 'hits above', str(identity_cutoff)+'%', 'cutoff.'

        ''' Save bidirectional best hits, stored as follows: 
            (blast1 query, blast2 query, lower identity, larger evalue) '''
        best_bidirectional_hits = __get_bidirectional_best__(forward_best_hits, reverse_best_hits)
        if verbose:
            print 'Found', len(best_bidirectional_hits), 'bi-directional best hits.'
        
    ''' Format into DataFrame and optionally save to file '''
    df_best_hits = pd.DataFrame(best_bidirectional_hits)
    if not output_file is None:
        df_best_hits.to_csv(output_file, sep='\t')
    return df_best_hits


def __run_commands__(args, log_to_console=True):
    ''' Runs a shell command and optionally prints stdout/stderr to console '''
    if log_to_console:
        print ' '.join(args)
        print subprocess.check_output(args, stderr=subprocess.STDOUT)
    else:
        subprocess.call(args)