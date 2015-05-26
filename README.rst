Introduction: what is the Ontology Oracle?
------------------------------------------

The **Ontology Oracle** is a Python module for collecting data for the
categorization of genes.  It is specifically designed to gather ontological
information (Gene Ontology terms, Enzyme Commission numbers, etc.) to
enhance the analysis of RNA-Seq results.  The package also includes a
number of Python helper scripts and R scripts for generating plots to
summarize the collected data.

Prerequisites
-------------

The **Ontology Oracle** has been tested with Python 2.7.  It requires the
**pandas** library; it is recommended that the
`Anaconda Scientific Python Distribution <https://store.continuum.io/cshop/anaconda/>`_ 
is installed to ensure that all prerequisites are satisfied.


Installing
----------

The source distribution for the most recent version can be obtained from the
`Ontology Oracle project page <https://github.com/mstrosaker/ontology_oracle>`_ 
by clicking on the Download ZIP button.  The module can be installed with::

    > sudo python setup.py install

(sudo is required to create the /opt directory for ancillary scripts and
plotting tools.)

How to use it?
--------------

The **Ontology Oracle** operates on mappings of locus tags to gene
accessions, and mines online databases such as NCBI and Uniprot for
ontological data associated with the locus tags.  This data can then be
annotated with the raw RPKM numbers from an RNA-Seq run.  At that point,
an assessment of the amount of upreguation or downregulation per ontological
category can be performed.

For a detailed usage guide and an example, please consult the
`user's guide <https://github.com/mstrosaker/ontology_oracle/wiki/User's-guide>`_.

License
-------

Copyright (c) 2015 Michael Strosaker.  See the LICENSE file for license
rights and limitations (MIT).

