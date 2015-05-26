# Copyright (c) 2015 Michael Strosaker
# MIT License
# http://opensource.org/licenses/MIT

import utils
from BeautifulSoup import BeautifulSoup, Tag

class _go:
    def __init__(self, obo):
        self.id = ''
        self.name = ''
        self.definition = ''
        self.synonyms = []
        self.xrefs = []
        self.isa = []
        self.ancestors = []
        self._children = []
        self.namespace = ''
        self.slims = []
        self.obo = obo.split('\n')
        for line in self.obo:
            if line.startswith('id: '):
                self.id = line[4:]
            elif line.startswith('name: '):
                self.name = line[6:]
            elif line.startswith('def: '):
                self.definition = line[5:].strip('"')
            elif line.startswith('synonym: '):
                self.synonyms.append(line[9:])
            elif line.startswith('xref: '):
                self.xrefs.append(line[6:])
            elif line.startswith('is_a: '):
                entries = line[6:].split('!')
                self.isa.append(entries[0].strip())

        current = []
        current.extend(self.isa)
        while len(current) > 0:
            parent = go(current.pop(0))
            if len(parent.isa) == 0:
                self.namespace = parent.name
            self.ancestors.append(parent.id)
            current.extend(parent.isa)
        self.ancestors = list(set(self.ancestors))

        if len(self.ancestors) == 0:
            self.namespace = self.name

        for current in self.ancestors:
            if current in cache_go_slims:
                self.slims.append(current)
        self.slims = list(set(self.slims))

    def print_ancestry(self):
        print '%s %s' % (self.id, self.name)
        for isa in self.isa:
            parent = go(isa)
            print '  -> %s %s' % (parent.id, parent.name)

        for term in self.ancestors:
            parent = go(term)
            print '%s %s' % (parent.id, parent.name)
            for isa in parent.isa:
                ancestor = go(isa)
                print '  -> %s %s' % (ancestor.id, ancestor.name)

    @property
    def children(self):
        # this is nasty, but the children are not encoded in the OBO
        if len(self._children) == 0:
            file = utils.download('%smini' % (url_go_lookup[:-3] % self.id))
            soup = BeautifulSoup(file)
            tab = soup.first('table')
            for entry in tab.contents:
                if isinstance(entry, Tag):
                    self._children.append(entry.findAll('a')[1].contents[0])

        return self._children

    def children_detail(self):
        l = []
        for child in self.children:
            c = go(child)
            l.append('%s %s' % (c.id, c.name))
        return l
                


url_go_lookup = 'http://www.ebi.ac.uk/QuickGO/GTerm?id=%s&format=obo'
cache_go = {}

url_go_slims = 'http://www.geneontology.org/ontology/subsets/goslim_generic.obo'
cache_go_slims = []

def _populate_go_slims():
    global cache_go_slims

    if len(cache_go_slims) == 0:
        file = utils.download(url_go_slims)

        for line in file.split('\n'):
            if line.startswith('id: GO:'):
                cache_go_slims.append(line[4:])

def go(term):
    global cache_go

    if term not in cache_go:
        _populate_go_slims()
        file = utils.download(url_go_lookup % term)
        cache_go[term] = _go(file)

    return cache_go[term]


url_ec2go = 'http://geneontology.org/external2go/ec2go'
local_ec2go = 'local_cache/ec2go.txt'
cache_ec2go = None

def ec2go(ec):
    global cache_ec2go
    if cache_ec2go is None:
        #print 'Building ontology mapping: EC numbers -> GO terms'
        try:
            file = utils.download(url_ec2go)
        except:
            # open from the local cache
            with open(local_ec2go, 'r') as f:
                file = f.read()

        cache_ec2go = {}
        for line in file.split('\n'):
            if line.startswith('!') or line.strip() == '':
                continue
            entries = line.split(' ')
            e = entries[0].split(':')
            if e[1] not in cache_ec2go:
                cache_ec2go[e[1]] = []
            cache_ec2go[e[1]].append(entries[-1])

    # accept either a string or a list of strings as a parameter
    if isinstance(ec, basestring):
        ec = [ec]

    go_terms = []
    for e in ec:
        if e.startswith('ec:') or e.startswith('EC:'):
            e = e[3:]

        if e in cache_ec2go:
            go_terms.extend(cache_ec2go[e])

    return go_terms

url_pfam2go = 'http://geneontology.org/external2go/pfam2go'
local_pfam2go = 'local_cache/pfam2go.txt'
cache_pfam2go = None

def pfam2go(pfam):
    global cache_pfam2go
    if cache_pfam2go is None:
        #print 'Building ontology mapping: Pfams -> GO terms'
        try:
            file = utils.download(url_pfam2go)
        except:
            # open from the local cache
            with open(local_pfam2go, 'r') as f:
                file = f.read()

        cache_pfam2go = {}
        for line in file.split('\n'):
            if line.startswith('!'):
                continue
            entries = line.split(' ')
            if entries[0] not in cache_pfam2go:
                cache_pfam2go[entries[0]] = []
            cache_pfam2go[entries[0]].append(entries[-1])

    # accept either a string or a list of strings as a parameter
    if isinstance(pfam, basestring):
        pfam = [pfam]

    go_terms = []
    for pf in pfam:
        if pf.startswith('pfam'):
            pf = 'Pfam:PF%s' % pf[4:]
        elif pf.startswith('PFAM'):
            pf = 'Pfam:PF%s' % pf[4:]
        elif pf.startswith('PF'):
            pf = 'Pfam:PF%s' % pf[2:]

        if pf in cache_pfam2go:
            go_terms.extend(cache_pfam2go[pf])

    return go_terms

url_tigrfam2go = 'http://geneontology.org/external2go/tigrfams2go'
local_tigrfam2go = 'local_cache/tigrfam2go.txt'
cache_tigrfam2go = None

def tigrfam2go(tigrfam):
    global cache_tigrfam2go
    if cache_tigrfam2go is None:
        #print 'Building ontology mapping: TIGRFAMs -> GO terms'
        try:
            file = utils.download(url_tigrfam2go)
        except:
            # open from the local cache
            with open(local_tigrfam2go, 'r') as f:
                file = f.read()

        cache_tigrfam2go = {}
        for line in file.split('\n'):
            if line.startswith('!') or line.strip() == '':
                continue
            entries = line.split(' ')
            tf = entries[0].split(':')
            if tf[1] not in cache_tigrfam2go:
                cache_tigrfam2go[tf[1]] = []
            cache_tigrfam2go[tf[1]].append(entries[-1])

    # accept either a string or a list of strings as a parameter
    if isinstance(tigrfam, basestring):
        tigrfam = [tigrfam]

    go_terms = []
    for tf in tigrfam:
        if tf in cache_tigrfam2go:
            go_terms.extend(cache_tigrfam2go[tf])

    return go_terms

url_smart2go = 'http://geneontology.org/external2go/smart2go'
local_smart2go = 'local_cache/smart2go.txt'
cache_smart2go = None

def smart2go(smart):
    global cache_smart2go
    if cache_smart2go is None:
        #print 'Building ontology mapping: SMART entries -> GO terms'
        try:
            file = utils.download(url_smart2go)
        except:
            # open from the local cache
            with open(local_smart2go, 'r') as f:
                file = f.read()

        cache_smart2go = {}
        for line in file.split('\n'):
            if line.startswith('!') or line.strip() == '':
                continue
            entries = line.split(' ')
            s = entries[0].split(':')
            if s[1] not in cache_smart2go:
                cache_smart2go[s[1]] = []
            cache_smart2go[s[1]].append(entries[-1])

    # accept either a string or a list of strings as a parameter
    if isinstance(smart, basestring):
        smart = [smart]

    go_terms = []
    for s in smart:
        if s.startswith('smart'):
            s = 'SM%s' % s[5:]

        if s in cache_smart2go:
            go_terms.extend(cache_smart2go[s])

    return go_terms

url_interpro2go = 'http://geneontology.org/external2go/interpro2go'
local_interpro2go = 'local_cache/interpro2go.txt'
cache_interpro2go = None

def interpro2go(interpro):
    global cache_interpro2go
    if cache_interpro2go is None:
        #print 'Building ontology mapping: InterPro entries -> GO terms'
        try:
            file = utils.download(url_interpro2go)
        except:
            # open from the local cache
            with open(local_interpro2go, 'r') as f:
                file = f.read()

        cache_interpro2go = {}
        for line in file.split('\n'):
            if line.startswith('!') or line.strip() == '':
                continue
            entries = line.split(' ')
            s = entries[0].split(':')
            if s[1] not in cache_interpro2go:
                cache_interpro2go[s[1]] = []
            cache_interpro2go[s[1]].append(entries[-1])

    # accept either a string or a list of strings as a parameter
    if isinstance(interpro, basestring):
        interpro = [interpro]

    go_terms = []
    for i in interpro:
        if i in cache_interpro2go:
            go_terms.extend(cache_interpro2go[i])

    return go_terms

