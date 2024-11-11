configuration
=============

the SWIBRID pipeline is configured via entries in the ``config.yaml`` file

.. literalinclude:: ../../swibrid/config.yaml

for a more in-depth control over parameters to individual functions, have a look at the snakemake file ``pipeline.snake``. e.g., resource requirements are specified there, or which exact plots should be produced. 

for extended testing, the ``config.yaml`` draws on a section specifying input parameters for generation of synthetic reads::

        SIMULATION_PARAMS:
          input_clones: 'input/HD_clones.bed'    # bed file with coordinates of switch segments for input clones
          input_pars: 'input/input_pars.par'     # parameters for adding mutations and insertions/deletions
          model: "poisson"                       # distribution of clone sizes
          mix_n10_s1_10k:                        # parameters for specific sample (name should appear in SAMPLES)
            nclones: 10                          # number of input clones
            nreads: 1000                         # total number of input reads
            seed: 1                              # random seed
          mix_n1000_s2_10k:
            nclones: 1000
            nreads: 10
            seed: 2



