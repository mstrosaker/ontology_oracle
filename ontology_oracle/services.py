# Copyright (c) 2015 Michael Strosaker
# MIT License
# http://opensource.org/licenses/MIT

import json
from urllib2 import quote
import utils, accession

class uniprot:
    @classmethod
    def search(self, term):
        url = 'http://www.uniprot.org/uniprot/?query=%s&format=tab&sort=score' \
                    '&%s' % (quote(term),
                    quote('columns=id,reviewed,protein names'))
        results = utils.download(url)
        return [x.split('\t')[0] for x in results.split('\n')[1:] if x != '']

    @classmethod
    def accession(self, id):
        return accession.accession('uniprot', id)

class ncbi:
    @classmethod
    def protein_accession(self, id):
        return accession.accession('ncbi:protein', id)

    @classmethod
    def protein_search(self, term):
        # this returns GI numbers rather than accessions
        url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?' \
              '&db=protein&rettype=seqid&sort=relevance&retmax=20&' \
              'retmode=json&term=%s' % quote(term)
        results = utils.download(url)
        j = json.loads(results)
        if 'esearchresult' in j:
            if 'idlist' in j['esearchresult']:
                return j['esearchresult']['idlist']
        return None


def mine_protein(accession, database, gene_name=None, organism=None):
    acc = None
    if database == 'ncbi':
        acc = ncbi.protein_accession(accession)
    elif database == 'uniprot':
        acc = uniprot.accession(accession)

    if acc is None:
        return None

    if gene_name is None:
        # see if the gene name is in the accession
        if len(acc.gene) > 0:
            gene_name = acc.gene

    if gene_name is None:
        return acc

    # see if we can find it in another database
    supplement = None
    if database == 'ncbi':
        term = 'gene:%s organism:%s' % (gene_name, organism)
        results = uniprot.search(term)
        if results is not None and len(results) > 0:
            supplement = uniprot.accession(results[0])
    elif database == 'uniprot':
        term = '%s %s' % (gene_name, organism)
        results = ncbi.protein_search(term)
        if results is not None and len(results) > 0:
            supplement = ncbi.protein_accession(results[0])

    if supplement is not None:
        # only use if protein sequence match is exact
        if acc.sequence.upper() == supplement.sequence.upper() and \
                    len(acc.sequence) > 0:
            acc.enrich_ontology(supplement)

    return acc

#class ontology_table:
#    def __init__(self, progress=True):
#
#    def dump(filename):

