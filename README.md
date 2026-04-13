# deep_phylo
  
This tool is built on components of ***PhyloKit***, a toolkit originally designed for various evolutionary analyses of large protein superfamilies. ***deep_phylo*** implements a workflow for exploring natural protein sequence space in the endeavour of curating representative, evolutionarily contiguous sequence datasets for robust phylogenetic analysis under deep ancestral nodes.

The workflow is iterative and involves successive profile HMM (pHMM)-based searches, phylogenetic reconstruction, and refinement of pHMMs for each subsequent iteration. A novel procedure is introduced for the identification of monophyletic segments of sequence space for which representations (profiles) are internally robust, defined by their stable recovery of held out sequence members over more distantly related homologs. Refinement of profiles under this objective enables progressive, phylogenetically grounded exploration rather than global, unconstrained expansion, supporting the curation of representative, contiguous datasets for downstream inference.
  
> **Note:** The interface made available here is preliminary and does not expose all functionality and customisation. A more comprehensive and fully featured CLI is currently under development, as well as various optimisation of the underlying framework.  
  
  
## Installation  
 
 It is recommended to install within a dedicated environment (e.g. conda).
  
Clone the repository and install:  
  
```bash  
git clone https://github.com/samdavis16/deep_phylo.git
cd deep_phylo  
pip install .
```

## Non-Python dependencies

Various tools are used in the construction and processing of intermediate alignments, phylogenies, and profile HMMs. These are invoked via subprocess and must be callable from the command line. Currently, this requires separate installation.

- **[HMMER](http://hmmer.org/)**
For various operations with profile HMMs.

- **[MMseqs2](https://github.com/soedinglab/mmseqs2)**
For fast pairwise sequence searching and clustering of curation hits.

- **[MAFFT](https://mafft.cbrc.jp/alignment/software/source.html)**
For constructing multiple sequence alignments.

- **[trimAl](https://trimal.readthedocs.io/en/latest/installation.html)**
For simple MSA trimming by column occupancy.

- **[FastTree](https://morgannprice.github.io/fasttree/)** 
For the construction of intermediate phylogenies.

Customisation for some major steps of the workflow (e.g. more choices of aligner, tree inference tool etc.) will soon be made available.

## Usage

This version does not currently expose all features and options.

```bash  
python -m deep_phylo \
	--name <RUN_NAME> \
	--hmm <HMM_FILE_1> <THRESHOLD_1> \
	--hmm <HMM_FILE_2> <THRESHOLD_2> \
	... \
	[--iters N_ITERS] \
	[--conv-prop CONVERGENCE_PROPORTION] \
	[--max-iters MAX_ITERS] \
	[--db-dir SEQ_DB_PATH] \
	[--cpu N_CPUS]
```

Entry to the iterative workflow is currently *via* initial pHMM(s) and associated threshold(s) only. Each should be declared separately as `--hmm <HMM_FILE> <THRESHOLD>`. Alternative entry from a starting tree to be exposed. 

The number of iterations can be fixed, **or** a "convergence proportion" can be specified. This defines the minimum proportion of new sequences curated in the current iteration (relative to previous iterations) before the exploration is considered to have converged. A maximum number of iterations can be specified for stopping if convergence has not occurred (default 10).

The path to a locally stored sequence database (one or several Fasta files) can be specified. If not present locally, the database is downloaded prior to the first iteration. This is generally the most efficient option as the workflow involves several profile searches of the sequence database.
Note: Only searching of the full UniParc database is currently exposed.