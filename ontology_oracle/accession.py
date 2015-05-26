# Copyright (c) 2015 Michael Strosaker
# MIT License
# http://opensource.org/licenses/MIT

import urllib2
import utils, ontology

class reference:
    def __init__(self):
        self.authors = ''
        self.title = ''
        self.number = ''
        self.position = ''
        self.comments = ''
        self.cross_references = ''
        self.group = ''
        self.location = ''

    def set(self, var, value):
        self.__dict__[var] = value.strip(' ;')

    def __str__(self):
        ret = []
        if self.title != '':
            ret.append(self.title)
        if self.authors != '':
            ret.append(self.authors)
        if self.location != '':
            ret.append(self.location)
        if self.position != '':
            ret.append('Cited for: %s' % self.position)
        return '  '.join(ret)

class accession:
    def __init__(self, db, id):
        self.db = db
        self.id = id
        self.references = []
        self.gene = ''
        self.sequence = ''
        self.definition = ''

        # ontology data
        self.ec = []
        self.pfam = []
        self.go = []
        self.tigrfam = []
        self.smart = []
        self.interpro = []

        if self.db == 'uniprot':
            self._populate_from_uniprot()
        elif self.db == 'ncbi:protein':
            self._populate_from_entrez('protein')

        # cleanup ontology data
        new_list = []
        for pf in self.pfam:
            if pf.startswith('PF'):
                new_list.append('pfam%s' % pf[2:])
            else:
                new_list.append(pf)
        self.pfam = list(set(new_list))

        self.ec = list(set(self.ec))
        self.tigrfam = list(set(self.tigrfam))
        self.smart = list(set(self.smart))
        self.interpro = list(set(self.interpro))

    def _lookup_go_terms(self):
        # determine if there are any go terms associated with ontology entries
        if len(self.pfam) > 0:
            self.go.extend(ontology.pfam2go(self.pfam))

        full_ecs = [x for x in self.ec if not x.endswith('-')]
        if len(full_ecs) > 0:
            self.go.extend(ontology.ec2go(full_ecs))

        if len(self.tigrfam) > 0:
            self.go.extend(ontology.tigrfam2go(self.tigrfam))

        if len(self.smart) > 0:
            self.go.extend(ontology.smart2go(self.smart))

        if len(self.interpro) > 0:
            self.go.extend(ontology.interpro2go(self.interpro))

        self.go = list(set(self.go))

    def _process_db_xref(self, entry):
        # handles database cross-reference in uniprot accessions
        refs = entry.split(';')
        if refs[0] == 'GO':
            self.go.append(refs[1].strip())
        elif refs[0] == 'InterPro':
            self.interpro.append(refs[1].strip())
        elif refs[0] == 'Pfam':
            self.pfam.append(refs[1].strip())
        elif refs[0] == 'TIGRFAMs':
            self.tigrfam.append(refs[1].strip())
        #elif refs[0] == 'SUPFAM':
        #elif refs[0] == 'EMBL':
        #elif refs[0] == 'Gene3D':
        #elif refs[0] == 'EnsemblBacteria':
        #elif refs[0] == 'SMR':
        #elif refs[0] == 'ProteinModelPortal':
        #elif refs[0] == 'PSIRF':
        #elif refs[0] == 'PANTHER':
        #elif refs[0] == 'KEGG':
        #elif refs[0] == 'PATRIC':

        self._lookup_go_terms()

    def _parse_2char_code_style(self):
        contents = {}
        for line in self.record.split('\n'):
            code = line[0:2]
            text = line[5:]
            if code in contents:
                contents[code].append(text)
            else:
                contents[code] = [text]

        for key, value in contents.iteritems():
            if key == 'ID':
                self.full_id = value[0]
            elif key == 'AC':
                self.accession_numbers = ''.join(value)
            elif key == 'DT':
                self.dates = value
            elif key == 'DE':
                self.description = ' '.join(value)
            elif key == 'GN':
                self.gene_name = ' '.join(value)
                if self.gene_name.startswith('Name='):
                    entries = self.gene_name.split(' ')
                    self.gene = entries[0][5:]
            elif key == 'OS':
                self.organism_species = ' '.join(value)
            elif key == 'OG':
                self.organelle = value[0]
            elif key == 'OC':
                self.organism_classification = ' '.join(value)
            elif key == 'OX':
                self.taxonomy_cross_reference = value[0]
            elif key == 'OH':
                self.organism_host = value[0]
            elif key == 'RN':
                if len(self.references) == 0:
                    self.references.append(reference())
                self.references[0].set('number', ' '.join(value))
            elif key == 'RP':
                if len(self.references) == 0:
                    self.references.append(reference())
                self.references[0].set('position', ' '.join(value))
            elif key == 'RC':
                if len(self.references) == 0:
                    self.references.append(reference())
                self.references[0].set('comments', ' '.join(value))
            elif key == 'RX':
                if len(self.references) == 0:
                    self.references.append(reference())
                self.references[0].set('cross_references', ' '.join(value))
            elif key == 'RG':
                if len(self.references) == 0:
                    self.references.append(reference())
                self.references[0].set('group', ' '.join(value))
            elif key == 'RA':
                if len(self.references) == 0:
                    self.references.append(reference())
                self.references[0].set('authors', ' '.join(value))
            elif key == 'RT':
                if len(self.references) == 0:
                    self.references.append(reference())
                self.references[0].set('title', ' '.join(value))
            elif key == 'RL':
                if len(self.references) == 0:
                    self.references.append(reference())
                self.references[0].set('location', ' '.join(value))
            elif key == 'CC':
                self.comments = ' '.join(value)
            elif key == 'DR':
                self.database_cross_references = ' '.join(value)
                for entry in value:
                    self._process_db_xref(entry)
            elif key == 'PE':
                self.protein_existence = value[0]
            elif key == 'KW':
                self.keywords = ' '.join(value)
            elif key == 'FT':
                self.feature_table_data = ' '.join(value)
            elif key == 'SQ':
                self.sequence_data = value[0]
                if '  ' in contents:
                    self.sequence = ''.join(contents['  '])
                self.sequence = self.sequence.replace(' ', '')

    def _parse_genbank_line(self, line, keylen=12):
        if line == '//':
            return None, None
        return (line[0:keylen].strip(), line[keylen:])

    def _parse_genbank_stanza(self, pos, joinchar=' ', keylen=12):
        lines = self.record.split('\n')
        code, value = self._parse_genbank_line(lines[pos], keylen)
        if code is None:
            return (pos, '')
        l = [value]
        pos += 1
        code, value = self._parse_genbank_line(lines[pos], keylen)
        while code is not None and code == '':
            l.append(value)
            pos += 1
            code, value = self._parse_genbank_line(lines[pos], keylen)
        return (pos, joinchar.join(l))

    def _extract_ontology_data(self):
        # the features section is in genbank-style entries
        if 'features' in self.__dict__:
            features = self.features.split('\n')
            note = None
            for i in range(len(features)):
                line = features[i].strip()
                # some notes span multiple lines
                if note is not None:
                    note = '%s %s' % (note, line)
                    if not note.endswith('"'):
                        continue
                if line.startswith('/note='):
                    note = line
                    if not note.endswith('"'):
                        continue
                if note is not None:
                    # check for pfam, TIGR, smart, etc. entries in the notes
                    note = note[6:].strip('"')
                    vals = note.split(';')
                    for val in vals:
                        val = val.strip()
                        if val.startswith('pfam'):
                            self.pfam.append(val)
                        elif val.startswith('TIGR'):
                            self.tigrfam.append(val)
                        if val.startswith('smart'):
                            self.smart.append(val)
                    note = None
                if line.startswith('/EC_number='):
                    vals = line.split('=')
                    self.ec.append(vals[1].strip('"'))
                if line.startswith('/gene='):
                    vals = line.split('=')
                    self.gene = vals[1].strip('"')

        self._lookup_go_terms()

    def _parse_genbank_style(self):
        lines = self.record.split('\n')
        pos = 0
        while pos < len(lines):
            code, value = self._parse_genbank_line(lines[pos])

            if code == 'LOCUS':
                self.locus = value
                pos += 1
            elif code == 'DEFINITION':
                pos, self.definition = self._parse_genbank_stanza(pos)
            elif code == 'ACCESSION':
                self.accession = value
                pos += 1
            elif code == 'VERSION':
                self.version = value
                pos += 1
            #elif code == 'DBLINK':
            #elif code == 'DBSOURCE':
            #elif code == 'KEYWORDS':
            #elif code == 'SOURCE':
            #elif code == 'REFERENCE':
            elif code == 'COMMENT':
                pos, self.comment = self._parse_genbank_stanza(pos, '\n')
            elif code == 'FEATURES':
                pos += 1
                pos, self.features = self._parse_genbank_stanza(pos, '\n',
                            keylen=1)
            elif code == 'ORIGIN':
                pos += 1
                pos, self.sequence_data = self._parse_genbank_stanza(pos,
                            '\n', keylen=1)
                seq = [x[9:] for x in self.sequence_data.split('\n')]
                self.sequence = ''.join(seq).replace(' ', '')
            else:
                pos += 1

        self._extract_ontology_data()

    def _populate_from_uniprot(self):
        url = 'http://www.uniprot.org/%s' % \
                    urllib2.quote('uniprot/%s.txt' % self.id)
        self.record = utils.download(url)
        self._parse_2char_code_style()

    def _populate_from_entrez(self, entrez_db):
        url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?' \
                    'db=%s&id=%s&rettype=gb&retmode=text' % (entrez_db,
                    urllib2.quote(self.id))
        self.record = utils.download(url)
        self._parse_genbank_style()

    def enrich_ontology(self, acc):
        '''
        Includes the ontology information from another accession into this
        one.  This should only be done with accessions referring to the
        same gene/protein, presumably obtained from different databases.
        '''
        self.ec.extend(acc.ec)
        self.ec = list(set(self.ec))

        self.pfam.extend(acc.pfam)
        self.pfam = list(set(self.pfam))

        self.go.extend(acc.go)
        self.go = list(set(self.go))

        self.tigrfam.extend(acc.tigrfam)
        self.tigrfam = list(set(self.tigrfam))

        self.smart.extend(acc.smart)
        self.smart = list(set(self.smart))

    def print_record(self):
        for line in self.record.split('\n'):
            print line


