
TIDAL1.2 Documentation
======================

For code visit `GitHub <https://github.com/laulabbumc/TIDAL1.2>`_


About TIDAL1.2
--------------
TIDAL is pipeline built to identify Transposon Insertion and Depletion in flies. The pipeline uses shell scrips, PERL, C and a combination of other publicly available tools. TIDAL uses a split read approach to identify Transposon Insertion and Depletion sites, and addition steps are implemented to reduce false positives.

TIDAL1.2 includes a few additional features

- Incorporation of IGEs (Immobile Genetic Elements) and tracking there insertion/deletion frequency
- Generate counts for transposon as proxy for copy number (not part of the default TIDAL output)
- update dependencies for the pipeline

Features
--------
- Flexible (TIDAL can be run with SRA files, Paired End and Single End libraries)
- Annotation Rich (TIDAL outputs have detailed annotation)
- Increased Specificity (Less false positives in outputs)
- Calculate coverage ratio to determine heterogenity of sites


Contents
--------

.. toctree::
   :maxdepth: 2

   installation
   run
   outputs

Support
-------
Please use Github issues to bring up any errors that occur with software.

License
-------
The project is licensed under the BSD license.

Acknowledgement
---------------
Thanks to David Tang for his blat parser script (`psl_to_bed_best_score.pl <http://davetang.org/muse/2012/05/15/using-blat/>`_)
