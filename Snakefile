#!/usr/bin/env python3


def GetNewFa(wildcards):
    return(new_fa[wildcards.spec])


def GetOrg(wildcards):
    return(orgs[wildcards.spec])


def GetOrigFa(wildcards):
    return(orig_fa[wildcards.spec])


def GetOrigGff(wildcards):
    return(orig_gff[wildcards.spec])


def GetPrefix(wildcards):
    return(prefix[wildcards.spec])


def GetTbl(wildcards):
    return(tbls[wildcards.spec])


# globals
new_fa = {
    'Vgerm': 'data/fasta/Vgerm_24_August_2020.assembly.annotated.fna',
    'Vpens': 'data/fasta/Vpens_24_August_2020.assembly.annotated.fna',
    'Vvulg': 'data/fasta/Vvulg_24_August_2020.assembly.annotated.fna'}
orgs = {
    'Vgerm': 'Vespula germanica',
    'Vpens': 'Vespula pensylvanica',
    'Vvulg': 'Vespula vulgaris'}
orig_fa = {
    'Vgerm': 'data/fasta/Vgerm.assembly.fna',
    'Vpens': 'data/fasta/Vpens.assembly.fna',
    'Vvulg': 'data/fasta/Vvulg.assembly.fna'}
orig_gff = {
    'Vgerm': 'data/gff/Vgerm_23_June_2020.gff3',
    'Vpens': 'data/gff/Vpens_23_June_2020.gff3',
    'Vvulg': 'data/gff/Vvulg_23_June_2020.gff3'}
prefix = {
    'Vgerm': 'HZH68',
    'Vpens': 'H0235',
    'Vvulg': 'HZH66'}
tbls = {
    'Vgerm': 'data/template/vg_template.sbt',
    'Vpens': 'data/template/vp_template.sbt',
    'Vvulg': 'data/template/vv_template.sbt'}


rule target:
    input:
        expand('output/munged/{spec}/ncbi/{spec}_ncbi.sqn',
               spec=['Vgerm', 'Vpens', 'Vvulg'])


rule table2asn_munged:
    input:
        fa = GetNewFa,
        tbl = GetTbl,
        gff = 'output/munged/{spec}/{spec}.gff3'
    params:
        org = GetOrg,
        prefix = GetPrefix,
        outdir = 'output/munged/{spec}/ncbi'
    output:
        val = 'output/munged/{spec}/ncbi/{spec}_ncbi.val',
        sqn = 'output/munged/{spec}/ncbi/{spec}_ncbi.sqn'
    shell:
        'bin/linux64.table2asn_GFF '
        '-M n -J -c -w -euk -Z -r '
        '-gaps-min 10 '
        '-j \"[organism={params.org}]\" '
        '-locus-tag-prefix \"{params.prefix}\" '
        '-i {input.fa} '
        '-f {input.gff} '
        '-outdir {params.outdir} '
        '-o {output.sqn} '
        '-augustus-fix '
        '-t {input.tbl} '
        '-l proximity-ligation '
        '-gap-type scaffold '

rule mung_gff:
    input:
        gff = 'output/ms_version/{spec}/{spec}.notrna.gff3',
        val = 'output/ms_version/{spec}/ncbi/ncbi.val',
        fai = lambda wildcards: GetNewFa(wildcards) + '.fai'
    output:
        gff = 'output/munged/{spec}/{spec}.gff3'
    log:
        'output/munged/{spec}/mung.log'
    script:
        'src/rename_gff.R'

rule table2asn_orig:
    input:
        fa = GetNewFa,
        tbl = GetTbl,
        gff = 'output/ms_version/{spec}/{spec}.notrna.gff3'
    params:
        org = GetOrg,
        prefix = GetPrefix,
        outdir = 'output/ms_version/{spec}/ncbi'
    output:
        val = 'output/ms_version/{spec}/ncbi/ncbi.val',
        sqn = 'output/ms_version/{spec}/ncbi/ncbi.sqn'
    shell:
        'bin/linux64.table2asn_GFF '
        '-M n -J -c -w -euk -Z -r '
        '-gaps-min 10 '
        '-j \"[organism={params.org}]\" '
        '-locus-tag-prefix \"{params.prefix}\" '
        '-i {input.fa} '
        '-f {input.gff} '
        '-outdir {params.outdir} '
        '-o {output.sqn} '
        '-augustus-fix '
        '-t {input.tbl} '
        '-l proximity-ligation '
        '-gap-type scaffold '

rule remove_trna:
    input:
        GetOrigGff
    output:
        'output/ms_version/{spec}/{spec}.notrna.gff3'
    shell:
        'grep -v '
        '\"tRNAScan-SE\" '
        '{input} > {output} ' 

rule index_fa:
    input:
        '{path}/{file}.{ext}'
    output:
        '{path}/{file}.{ext}.fai'
    wildcard_constraints:
        ext = 'fasta|fa|fna'
    shell:
        'samtools faidx {input}'
