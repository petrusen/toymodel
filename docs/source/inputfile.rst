Input file
=================================

This manual explains the parameters the *config.yaml* file which sets
the configuration for the kinetic simulation.


.. code-block:: yaml

   dkie:
     k1_16: 1
     k2_16: 0.1
     k3_16: 100
     k4_16: 0.000001
     kie1p: 1.02
     kie1s: 1.00
     kie2p: 1
     kie2s: 1
     kie3p: 1.02
     kie3s: 0.995
     kie4s: 1
     kie4p: 1
     dh2o: -9.59
     dcpo: 8.20
     dco: 0
     dpo: 0
     conc_h2o: 55.5
     conc_cpo: 0.004
     conc_co: 0
     conc_po: 0
     filtreactions: True
     verbose: False

   plot:
     4plotlayout: mk_mainplot.png
     C_time: mk_C_time.tsv  # tsv istead of csv because compounds contains ','
     dO_time: mk_dO_time.csv
     dO_rxnprogress: mk_dO_rxnprogress.csv
     rxnref: "CPO4"

   timerange:
     type: logspace
     start: -8
     stop: +4
     num: 100

----

Here we break down all the parameters that can be defined in the input file:

**1. dkie**

This section concerns the definition of the reaction rate constants (k1,k2,k3,k4) with primary and secondary kinetic isotopic effect. Furthermore, initial concentrations and \ δ\ :sup:`18`\ O can also defined.

- **k1_16** (``float``): rate constants for the k1 reaction involving a :sup:`16`\ O
- **k2_16** (``float``): rate constants for the k2 reaction involving a :sup:`16`\ O 
- **k3_16** (``float``): rate constants for the k3 reaction involving a :sup:`16`\ O
- **k4_16** (``float``): rate constants for the k4 reaction involving a :sup:`16`\ O
- **kie1p** (``float``): rate constants for the k1 reaction with primary effect involving a :sup:`18`\ O
- **kie1s** (``float``): rate constants for the k2 reaction with secundary effect involving a :sup:`18`\ O
- **kie2p** (``float``): rate constants for the k2 reaction with primary effect involving a :sup:`18`\ O 
- **kie2s** (``float``): rate constants for the k2 reaction with secundary effect involving a :sup:`18`\ O
- **kie3p** (``float``): rate constants for the k3 reaction with primary effect involving a :sup:`18`\ O
- **kie3s** (``float``): rate constants for the k3 reaction with secundary effect involving a :sup:`18`\ O
- **kie4p** (``float``): rate constants for the k4 reaction with primary effect involving a :sup:`18`\ O
- **kie4s** (``float``): rate constants for the k4 reaction with secundary effect involving a :sup:`18`\ O
- **dh2o** (``float``): \ δ\ :sup:`18`\ O for water
- **dcpo** (``float``): \ δ\ :sup:`18`\ O for organophosphate
- **dco** (``float``): \ δ\ :sup:`18`\ O for alcohol
- **dpo** (``float``): \ δ\ :sup:`18`\ O for phosphate
- **conc_h2o** (``float``): initial concentration of water
- **conc_cpo** (``float``): initial concentration of organophosphate
- **conc_co** (``float``): initial concentration of the alcohol
- **conc_po** (``float``): initial concentration of the phosphate
- **filtreactions** (``bool``): reduce the number of isotop combinations by not allowing scrambling in the transition state
- **verbose** (``bool``): print additional information during runtime

.. warning::

   Word of caution: reaction rate constants are based on the phenomenological isotopic findings.

**2. plot**

This section concerns the generation of the output files. The main results of the kinetic model are summarized in the PNG file, but the user can access to the raw data of the three plots in the PGN via the CSV and TSV files depicted below.

- **4plotlayout** (``str``): path where the main plot of the toy model will be saved
- **C_time** (``str``): path where the data of the conc vs time plot will be stored (as a tsv)
- **dO_time** (``str``): path where the data of the dO vs time plot will be stored (as a csv)
- **dO_rxnprogress** (``str``): path where the data of the dO vs rxnprogress will be stored (as a csv)
- **rxnref** (``str``): name of the compound (e.g., O, CO, PO4, CPO4, CPO5) which will be used to calculate the rxnprogress coordinate

**3. timerange**

This section concerns the definition of the time range which is used during the simulation.

- **type** (``str``): grid of the time steps
- **start** (``int``): minimum time in the kinetic simulation
- **stop** (``int``): maximum time in the kinetic simulation
- **num** (``int``): number of time steps in the kinetic simulation

.. warning::

   Word of caution: depending on the value of the reaction rate constants, the time range will need to be adjusted. For example, low rate constants
   will require higher maximum times, whereas high rate constants will require lower maximum times.
