#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 15:37:44 2019

@author: asmith
"""
import os
import sys
import numpy as np
import pandas as pd
from ruffus import *
from cgatcore import pipeline as P

# Put parameter YAML here
P.get_parameters('RNA_seq.yml')

@follows(mkdir('fastqc'))
@transform('*.fastq.gz', regex(r'(.*_.*).fastq.gz'), r'fastqc/\1_fastqc.zip')
def run_fastqc(infile, outfile):
    cmd = 'fastqc -q -t %(threads)s --nogroup %(infile)s --outdir fastqc/'
    P.run(cmd, job_queue=P.PARAMS['queue'], job_threads=P.PARAMS['threads'])



@follows(mkdir('bam'))
@collate('*.fastq.gz', regex(r'(.*)_[1|2].fastq.gz'), r'bam/\1.bam')
def align_fastq_paired(infiles, outfile):
    
    if not P.PARAMS['hisat_options']:
        P.PARAMS['hisat_options'] = ''
        
        fq1, fq2 = infiles
        bam = os.path.join('bam',
                           os.path.splitext(os.path.basename(fq1))[0].split('_')[0]
                           + '.bam')
        
        sorted_bam = bam.replace('.bam', '_sorted.bam')

        
        cmd = '''hisat2 -x %(hisat_index)s -1 %(fq1)s -2 %(fq2)s -p %(threads)s %(hisat_options)s |
                      samtools view -b - > %(bam)s &&
                      samtools sort -@ %(threads)s -m 5G -o %(sorted_bam)s %(bam)s &&
                      mv %(sorted_bam)s %(bam)s'''
        
        P.run(cmd, job_queue=P.PARAMS['queue'], job_threads=P.PARAMS['threads'])
    

@transform('*.fastq.gz', 
           regex(r'(.*)_(?![1|2])(.*).fastq.gz'),
           r'bam/\1_\2.bam')
def align_fastq_single(infile, outfile):
    
    fq = infile
    sorted_bam = outfile.replace('.bam', '_sorted.bam')
    
    cmd = '''hisat2 -x %(hisat_index)s -U %(fq)s  -p %(threads)s %(hisat_options)s |
              samtools view -b - > %(bam)s &&
              samtools sort -@ %(threads)s -m 5G -o %(sorted_bam)s %(bam)s &&
              mv %(sorted_bam)s %(outfile)s'''
    
    P.run(cmd, job_queue=P.PARAMS['queue'], job_threads=P.PARAMS['threads'])
    
    

@transform([align_fastq_paired, align_fastq_single], regex(r'bam\/(.*)'), r'bam/\1.bai')
def index_bam(infile, outfile):
    
    cmd = 'samtools index %(infile)s'
    P.run(cmd, job_queue=P.PARAMS['queue'], job_threads=P.PARAMS['threads'])

@follows(mkdir('picard_stats'))
@follows(index_bam)
@transform([align_fastq_paired, align_fastq_single],
           regex(r'bam/(.*.bam)'), 
           r'picard_stats/\1.picard_out.txt')
def get_picard_alignment_stats(infile, outfile):
    
    cmd = ' '.join(['picard CollectAlignmentSummaryMetrics',
                    'R=%(picard_reference)s',
                    'I=%(infile)s',
                    'O=%(outfile)s',]) 
    
    P.run(cmd, job_queue=P.PARAMS['queue'], job_threads=P.PARAMS['threads'])
    

@follows(mkdir('featureCounts'))
@merge([align_fastq_paired, align_fastq_single], f'featureCounts/{P.PARAMS["featureCounts_output"]}')
def count_reads(infiles, outfile):
    
    fnames = ' '.join(infiles)
    output = os.path.join('featureCounts', P.PARAMS['featureCounts_output'])
    
    cmd = '''featureCounts -a %(featureCounts_gtf)s 
             -o %(output)s -T %(threads)s  
             %(featureCounts_options)s %(fnames)s'''
             
    P.run(cmd, job_queue=P.PARAMS['queue'], job_threads=P.PARAMS['threads']) 

@follows(get_picard_alignment_stats)
@originate( 'multiqc/final_multiqc_report.html')
def run_multiqc(outfile):
    
    cmd = '''multiqc . -o multiqc -n final_multiqc_report.html'''
    P.run(cmd, job_queue=P.PARAMS['queue'], job_threads=1)

@follow(mkdir('bigwigs'))
@transform([align_fastq_single, align_fastq_paired],
           regex(r'bam/(.*.bam)'),
           r'bigwigs/\1_plus.bigWig')
def make_bigwig(infile, outfile):
    
    plus = outfile
    minus = outfile.replace('plus', 'minus')
    
    cmd = '''bamCoverage -b %(infile)s -o %(outfile)s --filterRNAstrand forward -p %(threads)s &&
             bamCoverage -b %(infile)s -o %(outfile)s --filterRNAstrand reverse -p %(threads)s'''
    
    P.run(cmd, job_queue=P.PARAMS['queue'], job_threads=P.PARAMS['threads'])

@follows(mkdir(hub_dir), 
         mkdir(assembly_dir), 
         make_bigwig)
@originate(os.path.join(hub_dir, 'hub.txt'))
def generate_hub_metadata(outfile):

    content = {'hub': P.PARAMS['hub_name'],
               'shortLabel': P.PARAMS['hub_short'] if P.PARAMS['hub_short'] else P.PARAMS['hub_name'],
               'longLabel': P.PARAMS['hub_long'] if P.PARAMS['hub_long'] else P.PARAMS['hub_name'],
               'genomesFile': 'genomes.txt',
               'email': P.PARAMS['hub_email'],
               'descriptionUrl': f'http://userweb.molbiol.ox.ac.uk/{P.PARAMS["hub_publoc"].strip("/")}',
               }

    with open(outfile, 'w') as w:
        for label, info in content.items():
            w.write(f'{label} {info}\n')

@transform(generate_hub_metadata, regex(r'.*.txt'), 'hub_address.txt')
def get_hub_address(infile, outfile):
    with open(outfile, 'w') as w:
        w.write(f'http://userweb.molbiol.ox.ac.uk/{os.path.abspath(infile).lstrip("/")}')

@follows(generate_hub_metadata)
@originate(os.path.join(hub_dir, 'genomes.txt'))
def generate_assembly_metadata(outfile):

    content = {'genome': P.PARAMS['hub_genome'],
               'trackDb': os.path.join(P.PARAMS['hub_genome'], 'trackDb.txt')}

    with open(outfile, 'w') as w:
        for label, info in content.items():
            w.write(f'{label} {info}\n')


@follows(generate_hub_metadata, mkdir(assembly_dir))
@merge(make_bigwig, 
       f'{assembly_dir}/trackDb.txt')
def generate_trackdb_metadata(infiles, outfile):
    def get_track_data(fn):
        return {'track': fn,
                'bigDataUrl': f'http://userweb.molbiol.ox.ac.uk/{(os.path.join(assembly_dir, fn)).lstrip("/")}',
                'shortLabel': fn,
                'longLabel': fn,
                'type': f'{fn.split(".")[-1]}',
                'autoscale': 'on',
                'windowingFunction': 'mean',
                }
     
    # Generate all separate tracks
    bigwig_tracks_all = [get_track_data(os.path.basename(fn)) for fn in infiles]
    
    # Add colours to tracks
    colors = sns.color_palette('husl', len(bigwig_tracks_all))
    for track, color in zip(bigwig_tracks_all, colors):
        track['color'] = ','.join([str(c * 255) for c in color])
    
    
    # Write track data separated
    with open(outfile, 'w') as w:
        for track in bigwig_tracks_all:
            for label, data in track.items():
                w.write(f'{label} {data}\n')
            # Need to separate each track with a new line
            w.write('\n')

@follows(generate_trackdb_metadata)
@transform([make_bigwig],
           regex(r'(peaks|bigwigs)/(.*)'),
           f'{assembly_dir}/' + r'\2')
def link_hub_files(infile, outfile):
        
    infile_fp = os.path.abspath(infile)
    os.symlink(infile_fp, outfile)
    
    
    
    
    
if __name__ == "__main__":
    sys.exit( P.main(sys.argv) )