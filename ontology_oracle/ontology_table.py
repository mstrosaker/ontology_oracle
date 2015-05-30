# Copyright (c) 2015 Michael Strosaker
# MIT License
# http://opensource.org/licenses/MIT

from utils import dataset, fold_change
from services import mine_protein
from ontology import go

class _feature_row:
    def __init__(self, feature, description, gene, protein_id):
        if isinstance(feature, basestring) and feature != '':
            self.feature = feature
        else:
            self.feature = None

        if isinstance(description, basestring) and description != '':
            self.description = description
        else:
            self.description = None

        if isinstance(gene, basestring) and gene != '':
            self.gene = gene
        else:
            self.gene = None

        if isinstance(protein_id, basestring) and protein_id != '':
            self.protein_id = protein_id
        else:
            self.protein_id = None

        self.expression_labels = []
        self.expression = {}
        self.foldchange_labels = []
        self.foldchanges = {}
        self.annotation_labels = []
        self.annotation = {}

    def retrieve_accession(self, organism, lookup_db):
        self.acc = mine_protein(self.protein_id, lookup_db, self.gene, organism)
        if self.gene is None:
            if self.acc is not None and len(self.acc.gene) > 0:
                self.gene = self.acc.gene
            else:
                self.gene = self.feature

        self.go = []
        self.go_slims = None
        self.ec = []
        self.pfam = []
        self.tigrfam = []
        self.smart = []
        self.interpro = []
        if self.acc is not None:
            self.go = self.acc.go
            self.ec = self.acc.ec
            self.pfam = self.acc.pfam
            self.tigrfam = self.acc.tigrfam
            self.smart = self.acc.smart
            self.interpro = self.acc.interpro

    def add_expression(self, label, value):
        self.expression_labels.append(label)
        self.expression[label] = value

    def calc_foldchange(self, from_label, to_label):
        label = '%s:%s' % (from_label, to_label)
        if from_label in self.expression_labels and \
                    to_label in self.expression_labels:
            self.foldchange_labels.append(label)
            self.foldchanges[label] = fold_change(self.expression[from_label],
                                                  self.expression[to_label])

    def add_annotation(self, label, value):
        self.annotation_labels.append(label)
        self.annotation[label] = value

    def csv(self, all_exprs, all_foldchanges, all_annotations):
        ret = []
        if self.gene:
            ret.append(self.gene)
            ret.append(self.feature)
        else:
            ret.append(self.feature)
            ret.append('')
        ret.append('"' + self.description + '"')
        ret.append(';'.join(self.go))

        if not self.go_slims:
            self.go_slims = []
            for term in self.go:
                self.go_slims.extend(go(term).slims)
            self.go_slims = list(set(self.go_slims))
        ret.append(';'.join(self.go_slims))

        ret.append(';'.join(self.ec))
        ret.append(';'.join(self.pfam))
        ret.append(';'.join(self.tigrfam))
        ret.append(';'.join(self.smart))
        ret.append(';'.join(self.interpro))
        for label in all_exprs:
            if label in self.expression_labels:
                ret.append(str(self.expression[label]))
            else:
                ret.append('')
        for label in all_foldchanges:
            if label in self.foldchange_labels:
                ret.append(str(self.foldchanges[label]))
            else:
                ret.append('')
        for label in all_annotations:
            if label in self.annotation_labels:
                ret.append('"%s"' % str(self.annotation[label]))
            else:
                ret.append('')

        return ','.join(ret)

class ontology_table:
    def __init__(self, organism=None, feature_table=None, locus_col=None,
                 gene_col=None, lookup_db=None, label_col=None,
                 accession_col=None, description_col=None,
                 locus_tag_prefix=None, progress=True, filename = None):

        self.expression_labels = []
        self.foldchanges = []
        self.annotation_labels = []

        if filename is not None:
            # build the table from an existing file
            self.feat_rows = []
            self.feat_tab = dataset(filename, 'csv')
            firstrow = True
            for row in self.feat_tab.rows:
                feat = _feature_row(row['locus'], row['product'],
                                    row['feature'], '')

                if isinstance(row['go-term'], basestring):
                    feat.go = row['go-term'].split(';')
                else:
                    feat.go = []

                if isinstance(row['go-slim'], basestring):
                    feat.go_slims = row['go-slim'].split(';')
                else:
                    feat.go_slims = []

                if isinstance(row['ec'], basestring):
                    feat.ec = row['ec'].split(';')
                else:
                    feat.ec = []

                if isinstance(row['pfam'], basestring):
                    feat.pfam = row['pfam'].split(';')
                else:
                    feat.pfam = []

                if isinstance(row['tigrfam'], basestring):
                    feat.tigrfam = row['tigrfam'].split(';')
                else:
                    feat.tigrfam = []

                if isinstance(row['smart'], basestring):
                    feat.smart = row['smart'].split(';')
                else:
                    feat.smart = []

                if isinstance(row['interpro'], basestring):
                    feat.interpro = row['interpro'].split(';')
                else:
                    feat.interpro = []

                for col in self.feat_tab.colnames:
                    if col.startswith('expr:'):
                        if firstrow:
                            self.expression_labels.append('"' + col + '"')
                        label = ':'.join(col.split(':')[1:])
                        feat.expression_labels.append(label)
                        feat.expression[label] = row[col]
                    elif col.startswith('foldchange:'):
                        if firstrow:
                            self.foldchanges.append('"' + col + '"')
                        label = ':'.join(col.split(':')[1:])
                        feat.foldchange_labels.append(label)
                        feat.foldchanges[label] = row[col]
                    elif col.startswith('annotation:'):
                        if firstrow:
                            self.annotation_labels.append('"' + col + '"')
                        label = ':'.join(col.split(':')[1:])
                        feat.annotation_labels.append(label)
                        feat.annotation[label] = row[col]

                self.feat_rows.append(feat)
                firstrow = False
            self.build_index()
            return

        if not organism or not feature_table or not locus_col or \
                    not gene_col or not lookup_db or not label_col or \
                    not accession_col or not description_col or \
                    not locus_tag_prefix:
            raise Exception('missing required parameter(s)')

        self.feat_tab = dataset(feature_table, 'tab-delimited')
        #index = self.feat_tab.index([locus_col, gene_col], accession_col)

        self.feat_rows = []
        finished = 0
        for row in self.feat_tab.rows:
            if isinstance(row[accession_col], basestring):
                feature = _feature_row(row[label_col], row[description_col],
                                       row[gene_col], row[accession_col])

                feature.retrieve_accession(organism, lookup_db)

                self.feat_rows.append(feature)
            finished += 1
            if progress and (finished % 10) == 0:
                print 'finished %4d records' % finished

        self.build_index()

    def build_index(self):
        self.index = {}
        for feat in self.feat_rows:
            self.index[feat.feature] = feat
            if feat.gene:
                self.index[feat.gene] = feat

    def dump(self, filename):
        outfile = open(filename, 'w')
        cols = ['feature', 'locus', 'product', 'go-term', 'go-slim', 'ec',
                'pfam', 'tigrfam', 'smart', 'interpro']
        cols.extend([('"expr:%s"' % x) for x in self.expression_labels])
        cols.extend([('"foldchange:%s"' % x) for x in self.foldchanges])
        cols.extend([('"annotation:%s"' % x) for x in self.annotation_labels])

        outfile.write('%s\n' % ','.join(cols))

        for f in self.feat_rows:
            outfile.write('%s\n' % f.csv(self.expression_labels,
                                      self.foldchanges, self.annotation_labels))

        outfile.close()

    def add_RNAseq_values(self, filename, label, locus_col, value_col):
        if label in self.expression_labels:
            raise Exception('Expression data label %s already in table' % label)
        expr_data = dataset(filename, 'csv')
        self.expression_labels.append(label)
        for row in expr_data.rows:
            if isinstance(row[locus_col], basestring):
                if row[locus_col] in self.index:
                    self.index[row[locus_col]].add_expression(label,
                                                         row[value_col])

    def calc_foldchange(self, from_label, to_label):
        label = '%s:%s' % (from_label, to_label)
        if label in self.foldchanges:
            raise Exception('Fold change label %s already in table' % label)
        self.foldchanges.append(label)
        for row in self.feat_rows:
            row.calc_foldchange(from_label, to_label)

    def add_annotation(self, filename, label, value_col, locus_col=None,
                       product_col=None):
        if label in self.annotation_labels:
            raise Exception('Annotation label %s already in table' % label)
        annot_data = dataset(filename, 'csv')
        self.annotation_labels.append(label)
        for row in annot_data.rows:
            if locus_col and isinstance(row[locus_col], basestring):
                if row[locus_col] in self.index:
                    self.index[row[locus_col]].add_annotation(label,
                                                         row[value_col])
            elif product_col and isinstance(row[product_col], basestring):
                for tab_row in self.feat_rows:
                    if tab_row.description and \
                                row[product_col] in tab_row.description:
                        row[product_col].add_annotation(label, row[value_col])

