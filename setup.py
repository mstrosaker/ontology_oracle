# Copyright (c) 2015 Michael Strosaker
# MIT License
# http://opensource.org/licenses/MIT

import sys
from distutils.core import setup

try:
    with open('../README.md', 'rt') as readme:
        description = readme.read()
except IOError:
    description = ''

setup(
    name='ontology_oracle',
    version='0.1.0',
    description='Ontology Oracle',
    long_description=description,
    license='MIT License',
    author='Mike Strosaker',
    author_email='mstrosaker@gmail.com',
    maintainer='Mike Strosaker',
    maintainer_email='mstrosaker@gmail.com',
    url='https://github.com/mstrosaker/ontology_oracle',
    download_url='https://github.com/mstrosaker/ontology_oracle',
    platforms='Cross Platform',
    classifiers = [
        'Environment :: Console',
        'Programming Language :: Python :: 2',
        'Operating System :: POSIX',
    ],

    packages = [ 'ontology_oracle' ],
    data_files = [('/opt/ontology_oracle/examples', [
                      'examples/sample_feature_table.txt',
                      'examples/sample_RNAseq.csv',
                  ]),
                  ('/opt/ontology_oracle/plotting', [
                      'plotting/ontology_tornado',
                      'plotting/heatmap',
                  ]),
                  ('/opt/ontology_oracle/plotting/helpers', [
                      'plotting/helpers/gen_tornado',
                      'plotting/helpers/gen_heatmap',
                  ]),
                 ],
)

