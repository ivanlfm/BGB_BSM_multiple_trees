######
README
######

	This script uses BioGeoBEARS to estimate (1) ancestral ranges and (2) timing and nature of biogeographic transitions by using stochastic maps, while accounting for phylogenetic uncertainty (topology and ages) by running these analysis over a sample of trees (e.g. trees from the posterior distribution of a Bayesian analysis).
	
	By Ivan L. F. Magalhaes (magalhaes@macn.gov.ar) and Martín J. Ramírez
	(c) 2021

##################
INPUTS AND OPTIONS
##################

	The inputs necessary are:
	1) A file with distribution data (as required by BioGeoBEARS)
	2) A nexus file with any number of phylogenetic trees (please note the script assumes you have already discarded the initial trees as burn-in; also note that MrBayes has a particular notation for branch lengths and you cannot use the trees directly from the .t files, you need to process them using burntrees first, see https://github.com/nylander/Burntrees, use the -myr option)
	3) Optionally, any modifier dispersal matrices / time stratification files / etc for BioGeoBEARS
	
	The script can be easily adapted to your own dataset by modifying the objects in lines ~39-164.
	Below I outline the options and explain them in detail.
	
	DIRECTORY and FILE names
	
		filename: the name that will be used for BioGeoBEARS and stochastic maps outputs
	
		dir_name: name of the directory where outputs will be saved
	
	RUN OPTIONS
	
		run_BGB: TRUE estimates ancestral ranges for each of the trees, FALSE loads results of previous analysis (note that results must have been saved as dir_name/BGB_results/filename_rep1.Rdata, dir_name/BGB_results/filename_rep2.Rdata...)
	
		run_BSM: TRUE runs biogeographic stochastic maps for each of the ancestral range estimates, FALSE loads results of previous analysis (note that results must have been saved as dir_name/BSM_results/BSM_rep1/RES_ana_events_tables.Rdata, etc.)
	
		n_stochastic_maps: number of stochastic maps to run for each tree (recommended: 50-100)
	
		n_cores: number of cores to be used in BioGeoBEARS likelihood estimations
	
		optimization: method for searching for optimal parameters during ML calculations
		FALSE for optim, TRUE for optimx, "GenSA" for GenSA
	
	TREES
	
		posterior_trees_fn: nexus file with the posterior trees. This file may contain any number of trees, provided it is larger than n_trees
	
		n_trees: the number of posterior trees to be sampled from the tree file and used for ancestral range estimates (recommended: 100)
	
		sample_trees: TRUE samples n_trees from trees file
		FALSE uses trees previously sampled and saved as /trees/tree1.nex
		use FALSE if you want to use exactly the same trees for different analysis (e.g. different models such as DEC or DIVA-like, or different time-stratification scenarios)
	
		fossil_tips: set this to TRUE if your tree includes fossils as terminals; this will correct cases in which branch lengths are too short for BioGeoBEARS
	
	DISTRIBUTIONS
	
		geogfn: file in Phylip format with geographic distributions of organisms, as required by BioGeoBEARS.
		IMPORTANT: the species in this file must be present in the trees, but not all the species in the trees must be in the distributions file. Any species absent in the distributions file, but present in the trees, will be pruned from the trees. This is useful for pruning outgroups, etc. But please check carefully that you are not inadvertedly pruning important species!
	
	OPTIONS FOR COUNTING EVENTS
	
		get_SD_common_events: TRUE calculates the standard deviation of the number of the most common biogeographic events across all trees and all stochastic maps. Somewhat slow!

		n_digits_age: number of decimal digits for age estimates in output tables
	
		time_period_size: size of the time periods to estimate the number of lineages through time by area
		E.g. 1 will count the number of lineages in each area in 1 million year intervals
	
		nodes_of_interest (and names(nodes_of_interest)): use this to track particular clades in the table with biogeographic transitions. Because the trees are not identical, the same clade will have different node numbers across trees. Define clades by listing ALL of the clade members, and give them a name. The rows corresponding to the first node of this clade will be identified in the BSM_transitions.csv output table.
	
		outgroup_taxa: defines a vector of outgroup taxa that will be excluded from the LTT by area plot, and from counts of biogeographic events. This is useful if you want to include an outgroup to properly estimate the range at the root of your ingroup, but want to count events only in the ingroup. This option can be left empty if you have no outgroups, or want to use them to count events ("outgroup_taxa <- c()").
	
	BIOGEOBEARS OPTIONS
	
		The model runs the DEC model by default (e.g. if DIVALIKE, BAYAREALIKE and J_founder_event are all set to FALSE).
	
		DIVALIKE: TRUE runs the DIVA-like model. (do not set DIVALIKE and BAYAREALIKE to TRUE simultaneously!)
		
		BAYAREALIKE: TRUE runs the BAYAREA-like model. (do not set DIVALIKE and BAYAREALIKE to TRUE simultaneously!)
		
		J_founder_event: TRUE adds a founder-event speciation (+J) free parameter to the model.
		
		time_stratified: TRUE runs a time stratified analysis
		
		time_stratified_fn: name of the input file containing the time slices, formatted as required by BioGeoBEARS
	
		dispersal_matrix: TRUE adds a dispersal matrix
		
		dispersal_matrix_fn: name of the input file containing the dispersal matrix, formatted as required by BioGeoBEARS
		
		area_adjacency: TRUE adds a file for area adjacency
		
		area_adjacency_fn: name of the input file containing the area adjacency, formatted as required by BioGeoBEARS
		
		distance_matrix: TRUE adds a distance matrix
		
		distance_matrix_fn: name of the input file containing the distance matrix, formatted as required by BioGeoBEARS
		
		areas_allowed: TRUE adds a file for areas allower
		
		areas_allowed_fn: name of the input file containing the areas allowed, formatted as required by BioGeoBEARS
		
		max_range_size: maximum allowed range size. Must be larger than max_observed_range_size
	
	PLOT OPTIONS
	
		pdf_w: width of the output PDF figures with estimated ancestral ranges
		
		pdf_h: height of the output PDF figures with estimated ancestral ranges
		
		plot_95: TRUE plots 95% intervals in the LTT by area plot
	
		plot_BSM: TRUE plots all the stochastic maps; this is slower and results in MANY files. 
		I recommend setting it to FALSE; you can re-run the script later to obtain the plots if you need them by setting run_BGB and run_BSM to FALSE and plot_BSM to TRUE
	
		plot_LTT: TRUE creates a plot of lineages through time by area ("LTT by area")
	
		n_maps_count_LTT: number of stochastic maps used to make the "LTT by area" plot. 
		Using all the maps is slower and does not make a big difference. I recommend using 10% of the stochastic maps
		(e.g. if n_stochastic_maps is 100, I set n_maps_count_LTT to 10)
	
#######
OUTPUTS
#######

	The script outputs:
	
	(1) R objects with the results of BioGeoBEARS ancestral range estimates (BGB_results)
	and corresponding figures (BGB_figures)
	(2) R objects with the results of the stochastic maps (BSM_results) and corresponding figures (BSM_figures)
	(3) All_transitions_by_node_BSM.csv: a table with biogeographic transitions for all nodes of all trees in all BSM replicates, including the identity of the node, its starting/ending ages, starting/ending geographic ranges, and the nature of the transition. Nodes of interested are highlighted in this table.
	(4) CountStatesByTime_all_reps.csv: the number of lineages in each area, in each time period, for each BSM replicate
	(5) Average_LTT_by_area.csv: the average number of lineages in each area through time
	(6) Average_LTT_by_area.emf: the information in (4) in graphic format
	(7) Common_transitions.csv: a summary of the most common transitions
	(8) Replication_parameters.csv: a summary of the LnL, AICc, and estimated parameters for each of the trees
	(9) RootState.csv: most likely state at the root for each of the replicates.

########
CITATION
########

	This scripts basically wraps functions from BioGeoBEARS (including ancestral range estimation and biogeographic stochastic maps), so make sure you cite the package and the associated papers:
	
		Matzke, Nicholas J. 2018. BioGeoBEARS: BioGeography with Bayesian (and likelihood) Evolutionary Analysis with R Scripts. version 1.1.1, published on GitHub on November 6, 2018. DOI: http://dx.doi.org/10.5281/zenodo.1478250
	
		Matzke N.J. 2013. Probabilistic historical biogeography: new models for founder-event speciation, imperfect detection, and fossils allow improved accuracy and model-testing. Front. Biogeogr. 5.
	
		Matzke N.J. 2014. Model selection in historical biogeography reveals that founder-event speciation is a crucial process in island clades. Syst. Biol. 63:951–970.
		
		Dupin J., Matzke N.J., Särkinen T., Knapp S., Olmstead R.G., Bohs L., Smith S.D. 2017. Bayesian estimation of the global biogeographical history of the Solanaceae. J. Biogeogr. 44:887–899.
	
	The rationale for counting events and plotting lineages through time while accounting for phylogenetic uncertainty has been developed in the papers below:
	
		Ceccarelli F.S., Koch N.M., Soto E.M., Barone M.L., Arnedo M.A., Ramírez M.J. 2019. The grass was greener: repeated evolution of specialized morphologies and habitat shifts in ghost spiders following grassland expansion in South America. Syst. Biol. 68:63–77.
		
		Magalhaes I.L.F., Santos A.J., Ramírez M.J. In preparation. Incorporating topological and age uncertainty into event-based biogeography supports palaeoislands in Galapagos and ancient connections among Neotropical dry forests.
	
	Also, please cite the original studies describing the models you use. Below is a non-exhaustive list:
	
		DEC: Ree et al. 2005 (https://doi.org/10.1111/j.0014-3820.2005.tb00986.x), Ree & Smith 2008 (https://doi.org/10.1080/10635150701883881)
		DIVA-like: Ronquist 1997 (https://doi.org/10.2307/2413643), Matzke 2013 (https://doi.org/10.1093/sysbio/syt040)
		BAYAREA-like: Landis et al. 2013 (https://doi.org/10.1093/sysbio/syt040), Matzke 2013 (https://doi.org/10.1093/sysbio/syt040)
		Founder-event speciation parameter (+J): Matzke 2013 (https://doi.org/10.1093/sysbio/syt040), Matzke 2014 (https://doi.org/10.1093/sysbio/syu056)
	
	
###############
VERSION HISTORY
###############

	v14:
		changed plotting
		added option to count events through time using BSM results -- this is not working too well yet
	v15:
		count states is now working fine with and without BSM
		started working on the summary table with most common events but did not finish it
	v16:
		finished summary table, corrected a few bugs and commented some parts of the code.
	v17:
		added AICc calculation... 
	v18:
		corrected a small bug in AICc calculation
		added option to skip plotting of BSMs
	v20:
		corrected some bugs, especially when reading a single tree
		fixed an important bug in using BSM to count lineages through time
	v23:
		only samples trees if BOTH sample_trees and run_BGB  are true.
		added code to adapt modifier matrices to shorter trees.
	v24:
		added axisPhylo to plots for node numbers.
		added some code to fix short branch lengths in sampled fossils.
	v25:
		fixed a serious bug regarding the creation of check_subset_sympatry table for counting events from the BSM (it only created a new table on the first replicate...)
	v26:
		corrected a small bug that prevented runs without time-slices (if(time_stratified) how_many_time_slices <- min(which(time_slices > current_height))
	v27:
		removed the control MCCT_tree (TRUE, FALSE) and replaced it by if (n_trees == 1)
	v28:
		added plot_LTT option
		removed all the options/code/outputs related to the transition table using only the most probable state (now the script needs BSM)
		commented the outputs
		add transition type to common transitions
	v29:
		calculates mean and SD of each type of transition among BSM, for each tree replicate
	v30:
		calculates SD for each of the most common transitions - still awfully slow!
	v31:
		removed some old objects no longer used (n_digits_prob)
		changed RUN_BSM from "RUN"/"LOAD" to TRUE/FALSE
		added get_SD_common_events
		now prints progress as SD is calculated
		added j and w to ReplicationParameters table
		fixed a small bug when counting the number of events for the ReplicationParameters table
		SD calculation for most common events still slow, but now is completed in a factible time
	v32:
		fixed a bug that happened when nodes_of_interest was left empty by adding a check (if(length(nodes_of_interest)>0))
	v33:
		added fossil_tips option
	v34: 
		fixed a bug that prevented dispersal matrices to be read if time_stratified was set to FALSE
	v35:
		fixed the portion that builds averageLTTbyarea to correct a bug that prevents the LTT plot from being output correctly if time_period_size is different from 1