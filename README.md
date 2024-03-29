This repository contains the family of scripts which integrate a yeast SBML model with experimentally constrained fluxes and then relate these fluxes to experimentally-measured metabolomics and proteomics

## Standard Workflow ##

1. Infer fluxes through metabolism
    * Use boundaryDataBlender.R to construct FBA-boundary fluxes - boundaryFluxes.Rdata.
    * Use FBA_run_full_reco.R to run QP-FBA and generate files for FVA.
    * Run FVA using the gurobi python interace (qp_fba_clust.py) on a cluster.
    * Use FBA_run_full_reco.R to combine QP-estimates and FVA flux bounds for each reaction.

2. Infer parametric forms for reactions
    * run reactionStructures.r to generate reaction-oriented relational information - creating lists of reaction forms and all relevent experimental and bioinformation information.
    * match_brenda.R - use scraped brenda files to determine supported regulators and their kinetic properties.

3. Evaluate parametric forms using MCMC and visualize results
    * Run part of FBGA.R to pare lists of reaction forms down to relevant reactions and pass to cluster.
    * Determine a posterior distribution of parameter values for each reaction form on the cluster using FluxOptimClus.R - computationally intensive, but highly parallelizable.
    * Use FBGA.R to assess the significance of proposed parameteric form alterations, and generate summaries of each reaction's kinetic properties, species and control.
    * Interactively comb through reaction summary information laid out in shinyapp in order to search for control principals.


---

## ChemicalSpeciesQuant ##

Relative or absolute quantification of various crucial species:

* boundaryDataBlender.R - Aggregates measurements of cellular composition and media depletion to summarize all of the major metabolic sources and sinks as specific fluxes
    * Specifies the costs associated with the production of each biomass component (such as GTP and ATP used for protein synthesis)
    * Generates boundaryFluxes.Rdata - a file which stand-alone summarizes each chemostats boundary flux expectation and coefficient of variation
* brauer-microarray: RNA abundance across 36 chemostat conditions
* ~~FTIR: attempt to measure composition using FTIR (not used)~~
* BulkComposition: measuring the dry weight and absolute amounts of abundant macromolecules which go into making more yeast
    * Protein
    * RNA
    * Carbohydrates
    * Polyphosphate
* chemostatInformation: Indicating which chemostat samples will be used for subsequent analysis and in the case of a subset of samples where boundary fluxes and composition fluxes were taken from different chemostats, the way in which these samples are combined is indicated
    * ListOfChemostats.txt - manually annotated data frame indicating for each chemostat where the relevant cell density, DR and NMR-based fluxes are found
    * ListOfChemostats_augmented.txt - ListOfChemostats.txt with cell density, volume and flow rate added from the relevant files
    * summarize_chemostat.R -  processes ListOfChemostats.txt into ListOfChemostats_augmented.txt
* Metabolites: Boer relative abundances of metabolites.  Yeast absolute concentration yifan absolute measures of metabolite abundance which need to be associated with Boer data using a paired comparison.
    * boerquant - reanalysis of boer et al. 2010 data in order to determine point estimate, SD and corr.
    * yeast_absolute_concentration_chemo.txt - scaling of metabolite relative abundance to absolute abundances.
       * Comparison of N15-labelled batch to C-limited chemostat in order to use metabolite concentrations from batch culture.
       * N15-labelled phosphate-limited chemostats compared to unlabelled spike-in amino acid standards. 
* Proteomics: Relative abundance of ~1200 proteins across 25 chemostat conditions.
* Lipids: Absolute quantification of abundant fatty acid.  Relative measurements of low abundance fatty acids.  



---

## ~~EcoliYeastMatch (Replaced by direct annotation using component-contribution)~~##

~~In order to incorporate thermodynamic information into the yeast flux-balance model, the gibbs free energy's available for ecoli rxns and species must be associated with the corresponding yeast ids.~~
~~This involves converting yeast ChEBI IDs to KEGG ids and then associating these yeast KEGG ids with ecoli KEGG ids~~

* ~~CTSpy.py (written by Dave Robinson): All conversions between metabolite names & standards were done using the chemical translation service (CTS).~~
* ~~ecoliYeastMatch.py: Because the same chemical can have multiple ChEBI or KEGG ids, the association is made more degenerate by first searching ChEBI ids to metabolite synonyms and then re-searching these against KEGG ids, to get a set of possible KEGG ids.~~
* ~~ecoliNameDict.txt - correspondence between geneID (unknown standard) and gene name. Downloaded from: http://coli.berkeley.edu/ecoli/gene.txt~~

* ~~gene-links.dat - correspondence between geneIDs and gene name, etc.  Downloaded from BioCyc~~

---
**~~gibbsAnalysis.R~~**

~~Applies two measures of thermodynamic reaction reversibility to predict the reversibility of yeast reactions~~

* ~~Match yeast genes -> E.C. and compare to e. coli genes -> E. C. to pass e. coli annotated free energy of reaction to yeast reactions.  This flips free energy to match substrates and products between corresponding yeast and e. coli reactions.~~
* ~~Predict free energy of reaction from free energy of formation.~~
* ~~Write revRxns.tsv, showing the free energy quantiles from these methods as well as the directionality from the Shlomi e. coli model ~~

---

## KEGGrxns ##

Relational files associating KEGG IDs and yeast systematic names
**keggIDparser.py:** Generates a tsv containing yeast systematic name \t SGD ID \t KEGG ID \t KEGG description and EC number

---

## MCA ##

Use reaction elasticities for adjacent reactions to calculate control coefficients for a subset of pathways

** MCA_SRH ** This script takes specified pathway models (currently hard-coded) and constructs 3-dimensional arrays populated by elasticities - following inversion these give pathway and concentration control coefficients.

---


**sbml_FBA_parser.py and SBML_genome_scale.py**

XML parsers which generate four files:
* The species involved in a reaction - proteins and metabolites
* The reaction, its stoichiometry and the compartment where it occurs
* Annotations of the protein and metabolite identity
* Mapping between a compartment designation and common name

## SBML_genomes ##

SBML reconstructions of several organisms from bioModels

---

## ~~Thermo (No longer used)~~##

~~Prepare the absolute concentrations of metabolites, so that these can be used to modify the reaction free energy~~
* ~~yeastNames.py - conversion between names of metabolites and possible ChEBI and KEGG ids~~

---

## Validations ##

Targeted experiments to follow up on the predictions of rOCA


## Yeast_genome_scale ##

This directory includes files for predicting flux using experimental data (FBA) and scripts for determining the extent to which this predicted flux can be associated with proteomics and metabolomics data (FBGA)

* FBA_run_full_reco.R: Integrate experimental data and in silico reconstruction
    * Experimental Data:
      1. Maximum nutrient uptake rates: Boer_nutrients
      2. Currently using an invariant composition: comp_yeat.tsv
    * Reconstruction:
      1. What direction do reactions proceed based on reconstruction annotation.
      2. which reactions can carry flux based on auxotrophy information
* FBA_lib.R: functions which provide more information about a reaction or metabolite and other misc.
* cytoGrapher.R: This program has two somewhat disparate tasks:
  1. Given previously positioned metabolite and reaction nodes, add on additional nodes to generate two positioning files: a node for each metabolite (some of these are split according to met_split.txt) and two nodes for each reaction, one linking substrates and one linking products.  In addition, some metabolite species can be classified as cofactors using the cofactorexception.txt file, so these nodes will be displayed primarily when they are being anabolized rather than when they are interconverted.
  2. Using positioning files and flux information: draw a metabolic network in cytoscape and color reactions according to the flux carried.
* metSty: Stores the positional information of metabolites nodes.
* rxnSty: Stores the positional information of reaction nodes.

**FBGA and FBGA_files:**
* FBGA.R: metabolites + enzyme ?= flux. Determine how well metabolite and enzyme experimental data can be aligned with flux levels.  Optimal fitted Km values (as well as a posterior distribution) are found through Metropolis-Hastings MCMC and non-negative least squares.  Essentially, likelihood is based on the consistency of predicted fluxes and measured fluxes given parameter values based upon a prior distribution.
* **FBGA_files:** 
    * brendaECtroll.py: using a list of E.C. numbers which contain largely measured species (enzymes and metabolites) determine previously measured values of Km, Ki and activators from BRENDA.  This is done both for S. cerevisiae and across all organisms
    * generated files with km, ki, and activators from S. cerevisiae BRENDA data, as well as an expanded set (all) from all organisms.


**companionFiles: Necessary files which are generated using companion scripts**

* all_activators.tsv, all_ki.tsv, all_km.tsv - activators, inhibitors and substrates/products for reactions - E.C. number to compound name
* Boer_nutrients.txt - chemostat media formulation
* BoerMetabolites.txt - raw supplemental file from Boer et al. 2010
* boerMean.tsv, boerSD.tsv, boerCorr.tsv - processed from boer data to yield by-condition point estimates, SD and correlation of residuals.
* cc_dG_python.tsv - free energies of reaction IDs from [https://github.com/eladnoor/component-contribution](Elad Noor's component-contribution)
* Yeast consensus reconstruction information - sbml-generated reaction/metabolite information ...
* customRxns.txt - custom reactions and metabolites which are not in reconstruction
* internalFluxMatches.csv - allows the specification of non-boundary fluxes (not currently used)
* proteinAbundance.tsv - Protein relative abundance
* proteinPrecision.tsv - Protein precision (inverse-variance)
* thermoAnnotate.txt - Manual specification of reaction directionality (in a few choice cases)
* yeast_absolute_concentration_chemo.txt - conversion between relative and absolute concentrations for some metabolites

**flux_cache: These files are generated during the standard analysis pipeline**



---

## Yeast_reconstruction ##

**Sequences:** 

* CHEBI_compound - a list of all of the metabolites referenced in the sbml reconstruction, associated with a standard name & structural information ?
* orf_genomic_1000.fasta - the mRNA sequence of all transcripts, including introns & 1000 bp upstream & downstream (couldn't get just the UTR)
* orf_genomic - from ATG to AUG/UAG/UAA including introns for all high-confidence transcripts
* rna_genomic - non-mRNA RNAs - tRNA, rRNA ?
* METAdict.tsv - correspondence between CHEBI IDs and consensus names
* PROTcomposition.tsv - Protein Name + ID with counts of amino acids
* RNAcomposition.tsv - Gene Name with counts of nucleotides


**generated by sbml_FBA_parser_forYeast.py:**
  
* comp_yeast.tsv: compartment IDs to compartments
* flux_dir_yeast.tsv: sbml-specified reaction reversibility
* rxn_par_yeast.tsv - reaction rID to enzyme complexes and KEGG reaction ID
* rxn_yeast.tsv - reactions and their associated stoichiometry as well as the enzyme (modifier) that catalyzes the rxn.  The internal names of all members are cross-referenced to spec_yeast.tsv and then associated with a meaningful CHEBI / uniprot or SGD designation in species_par_yeast.tsv
* spec_yeast.tsv: association between an internal id and a species id, that is unique for a metabolite or protein.  compartment information (extraceullar versus intracellular) is included.
* species_par_yeast.tsv: associates a species id with references to a meaningful id; either CHEBI/KEGG (for metabolites) or UniProt/SGD (for proteins).



---
