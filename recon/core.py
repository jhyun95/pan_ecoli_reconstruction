#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 14:36:29 2019

@author: jhyun95
"""


def reconstruct_with_cdhit(seq_fasta, work_dir=None, ref_dir='reference/'):
    '''
    Attempts a reconstruction with CD-Hit given a set of coding sequences 
    '''

    ''' Prepare output directory'''
    models_dir = ref_dir + 'ref_models/'
    genomes_dir = ref_dir + 'ref_genomes/'
    if work_dir is None:
        head, tail = os.path.split(seq_fasta)
        name = '.'.join(tail.split('.')[:-1])
        work_dir = name + '_recon/'
    if not os.path.isdir(work_dir):
        os.mkdir(work_dir)

    ''' Combine query fasta and reference fastas '''
    


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