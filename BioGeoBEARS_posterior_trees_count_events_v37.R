################################################################################################
###                   																		 ###
### Estimates ancestral ranges over a distribution of trees using BioGeoBEARS,               ###
### counts lineages in areas through time,                                                   ###
### identifies most common transitions,														 ###		
### and estimates the timing of transitions                                                  ###
### by using biogeographic stochastic maps.                               			         ###
###                                                                                          ###
### By Ivan L.F. Magalhaes & Martín J. Ramírez                                               ###
###                                                                                          ###
### Please refer to README.txt for instructions and information for citing					 ###
###                   																		 ###
################################################################################################
	start_time <- Sys.time()
	
	#load packages
		library(rexpokit)
		library(cladoRcpp)
		library(BioGeoBEARS)
			#if you don't have BioGeoBEARS, install it from github (you need Rtools and devtools installed):
			#install_github(repo="nmatzke/BioGeoBEARS")
		scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS")) # BioGeoBEARS needs this
		library(GenSA)
		library(FD)
		library(snow)
		library(parallel)
		library(ape)
		library(phytools)
		library(sjmisc)
		library(tictoc)
		library(dplyr)
	
#########################################
# create output directory and log file  #
#########################################

	#sets working directory
		setwd("C:/a/Sicarius/example")
	#analysis name to be appended to output files
		filename <- "Sicarius_DVL_strat_15"
	#name of the output directory	
		dir_name <- "example_Sicarius_DVL_strat_15"
		if(file.exists(dir_name)==FALSE) dir.create(dir_name)
	#logfile
		sink(paste(dir_name,"/",filename, "_logfile.txt", sep=""), append=TRUE, split=TRUE)
	
########################################################
# set inputs and some options for running the analyses #
########################################################
	
	### RUN OPTIONS
		#TRUE runs BioGeoBEARS to estimate ancestral ranges, FALSE loads previous runs
			run_BGB = TRUE
		#biogeographic stochastic mapping
			# TRUE runs n_stochastic_maps biogeographic stochastic maps for each tree below;
			# FALSE loads previous results
				run_BSM =  TRUE
			#number of BSM to run per tree
				n_stochastic_maps = 2
		#number of cores to be used in BioGeoBEARS likelihood estimations	
			n_cores <- 2
		#method for searching for optimal parameters during maximum likelihood calculations
			#FALSE for optim, TRUE for optimx, "GenSA" for GenSA
				optimization = TRUE
		#options for resuming a previous aborted run
			#to resume an analysis, set checkpoint_from_previous_run to TRUE, and indicate the number of replicates previously run
				checkpoint_from_previous_run = FALSE
				number_of_trees_run_previously = 0 #check logfile of previous run for "Completed replicate X/X" to set this correctly
	
	### TREES
		#file with the posterior trees. This script assumes you have already discarded the pre-burn-in trees
			#you can use the output of BEAST's LogCombiner directly
				posterior_trees_fn <- "example_Sicarius_200_trees.nex"
		#number of trees to be sampled from the posterior distribution
			# set this to 1 only if you are providing a file with a single consensus / maximum clade credibility tree,
			# or any other single, target tree
				n_trees = 10
		#sample_trees: TRUE samples n_trees from trees file, FALSE uses trees previously sampled and saved as /trees/tree1.nex, /trees/tree2.nex...
			sample_trees = TRUE
		#set fossil_tips to TRUE if your tree includes fossils as terminals; this will correct cases in which branch lengths are too short for BioGeoBEARS
			fossil_tips = FALSE
	### DISTRIBUTIONS
		# name of the file containing the distributions
			geogfn = "example_Sicarius_ranges.txt"
			tipranges = getranges_from_LagrangePHYLIP(geogfn)
		#gets the species names from the distribution file
		#species without range data will be dropped from the trees later, so check this carefully!
			species <- rownames(tipranges@df)
		
	### COUNT EVENTS OPTIONS
		# TRUE calculates SD for most common events (somewhat slow!)
			get_SD_common_events = TRUE
		# define the number of decimal digits for probabilities and ages in tables
			n_digits_age <- 2
		# define the size of the time slices to count lineages / area / time slice
			time_period_size <- 1
		# nodes of interest
			# if you are interested in tracking the timing and transitions of particular nodes in the tree,
			# define them here, and give names as you wish.
			# all taxa belonging to each clade MUST be listed, because in some trees it might not be monophyletic
			# this will tag the transitions table in 3 rows per node of interest (NoI):
			# ancestor -> NoI, NoI -> daughter1, NoI -> daughter2
			# if your nodes of interested are nested within each other, list them from most inclusive to less inclusive
				nodes_of_interest <- list(
					c("S_peruensis_506", "S_utriformis")
				)
				names(nodes_of_interest) <- c("Peru+Galapagos")
		# defines outgroup; this is important if you want to include outgroups for ancestral range estimation
			# but exclude them from the counts of transitions and lineages through time
			# by default this script assumes your tree contains only the ingroup taxa
			#leave this empty if your tree is represented solely by the ingroup
				outgroup_taxa <- c() 
			   
	### BIOGEOBEARS OPTIONS
		# The script will run the DEC model by default. Set the parameters below to TRUE modify it
		# Please not this script will run a single model each time, so do NOT set both DIVALIKE and BAYAREALIKE to TRUE simultaneously
			DIVALIKE = TRUE
			BAYAREALIKE = FALSE
			J_founder_event = FALSE
		
		#set these to TRUE and indicate the filenames for using dispersal matrices, time slices, etc.
			time_stratified = TRUE
			time_stratified_fn = "example_timeperiods_15.txt"
			dispersal_matrix = FALSE
			dispersal_matrix_fn = ""
			area_adjacency = FALSE
			area_adjacency_fn = ""
			distance_matrix = FALSE
			distance_matrix_fn = ""
			areas_allowed = TRUE
			areas_allowed_fn = "example_areaallowed.txt"
		
		# number of areas and maximum observed range size
			areas = names(tipranges@df)
			num_areas <- length(areas)
			max_observed_range_size <- max(rowSums(dfnums_to_numeric(tipranges@df)))
			
		# max range size (it can be equal to or higher than max_observed_range_size, but not lower)
			max_range_size <- 3
		
		# the total number of possible ranges depends both on the number of areas and max number of areas a species can occupy
		# should be less than 1000 to run analysis in under a day, less than 1500 to run in under a week, less than 2500 to run at all
			num_states <- numstates_from_numareas(numareas=num_areas, maxareas=max_range_size, include_null_range=TRUE)
			
	### PLOTTING OPTIONS
		#dimensions for the pdf output
			pdf_w <- 10
			pdf_h <- 15
		#plot error intervals in graphic of lineages through time by area, or not
			plot_95 = TRUE
		#plot figures for individual stochastic maps; TRUE is slower and results in MANY files! 
				plot_BSM <- TRUE
		# plot lineages through time by area or not
			plot_LTT <- TRUE
				# you can use n_maps_count_LTT to use only some of the maps for counting lineages through time
				# this speeds things considerably
					n_maps_count_LTT = 20

	### saves a file with all the options for your reference
	### perhaps it is a good idea to check this carefully before proceding to the next analyses
		all_options <- c(filename, dir_name, run_BGB, sample_trees, run_BSM, n_stochastic_maps, plot_BSM, plot_LTT, n_maps_count_LTT, n_cores, optimization, posterior_trees_fn, n_trees, paste(outgroup_taxa, collapse=","), geogfn, n_digits_age, time_period_size, DIVALIKE, BAYAREALIKE, J_founder_event, time_stratified, time_stratified_fn, dispersal_matrix, dispersal_matrix_fn, area_adjacency, area_adjacency_fn, distance_matrix, distance_matrix_fn, areas_allowed, areas_allowed_fn, num_areas, max_observed_range_size, max_range_size, num_states, pdf_w, pdf_h)
		options <- matrix(all_options, ncol = 1, nrow = length(all_options))
		rownames(options) = c("filename", "dir_name", "run_BGB", "sample_trees", "run_BSM", "n_stochastic_maps", "plot_BSM", "plot_LTT", "n_maps_count_LTT", "n_cores", "optimx", "posterior_trees_fn", "n_trees", "outgroup_taxa", "geogfn", "n_digits_age", "time_period_size", "DIVALIKE", "BAYAREALIKE", "J_founder_event", "time_stratified", "time_stratified_fn", "dispersal_matrix", "dispersal_matrix_fn", "area_adjacency", "area_adjacency_fn", "distance_matrix", "distance_matrix_fn", "areas_allowed", "areas_allowed_fn", "num_areas", "max_observed_range_size", "max_range_size", "num_states", "pdf_w", "pdf_h")
		write.table(options, file = paste(dir_name,"/analysis_options.txt",sep=""), col.names = F)

################################
# set a BioGeoBEARS run object #
################################
	
	# you should check this in detail but in most cases there is no need to change the parameters below
	# defines an object, inputs trees and ranges, sets various parameters for the program
		BioGeoBEARS_run_object = define_BioGeoBEARS_run() 					# Initialize a default model
		BioGeoBEARS_run_object$geogfn = geogfn 								# Give BioGeoBEARS the location of the geography text file
		BioGeoBEARS_run_object$max_range_size = max_range_size 				# Input the maximum range size
		BioGeoBEARS_run_object$min_branchlength = 0.000001 					# Min to treat tip as a direct ancestor (no speciation event)
		BioGeoBEARS_run_object$include_null_range = TRUE 					# set to FALSE for e.g. DEC* model, DEC*+J, etc.
		BioGeoBEARS_run_object$on_NaN_error = -1e50 						# returns very low lnL if parameters produce NaN error (underflow check)
		BioGeoBEARS_run_object$speedup = TRUE          						# shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
		BioGeoBEARS_run_object$use_optimx = optimization  					# if FALSE, use optim() instead of optimx()
		BioGeoBEARS_run_object$num_cores_to_use <- n_cores 					# use more cores to speed it up; this requires library(parallel) and/or library(snow)
		BioGeoBEARS_run_object$force_sparse = FALSE
		BioGeoBEARS_run_object$return_condlikes_table = TRUE 				# Good default settings to get ancestral states
		BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE # Good default settings to get ancestral states
		BioGeoBEARS_run_object$calc_ancprobs = TRUE    						# get ancestral states from optim run
		jmax = 2.99999														# max value for J parameter in DEC+J runs
		
		# adds files with time slices, dispersal matrices, etc. if needed
		if (time_stratified) {BioGeoBEARS_run_object$timesfn = time_stratified_fn}
		if (areas_allowed) {BioGeoBEARS_run_object$areas_allowed_fn = areas_allowed_fn}
		if (area_adjacency) {BioGeoBEARS_run_object$areas_adjacency_fn = area_adjacency_fn}
		if (distance_matrix) {BioGeoBEARS_run_object$distsfn = distance_matrix_fn}
		if (dispersal_matrix) {BioGeoBEARS_run_object$dispersal_multipliers_fn = dispersal_matrix_fn}
		files_to_read <- (time_stratified + areas_allowed + area_adjacency + distance_matrix + dispersal_matrix)> 0
			
		if (DIVALIKE){
			jmax = 1.99999
			# Set up DIVALIKE model
			# Remove subset-sympatry
				BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
				BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
				BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0
				BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
				BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
				BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
				BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"
			# Allow classic, widespread vicariance; all events equiprobable
				BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
				BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
				BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5
		}
		
		if (BAYAREALIKE){
			jmax = 0.99999
			# Set up BAYAREALIKE model
			# No subset sympatry
				BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
				BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
				BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0
			# No vicariance
				BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
				BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
				BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0
			# Adjust linkage between parameters
				BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
				BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
				BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"
			# Only sympatric/range-copying (y) events allowed, and with exact copying (both descendants always the same size as the ancestor)
				BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
				BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
				BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999
		}
		
		if (J_founder_event){
			BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
			BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.0001
			BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = 0.0001
			BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] = 0.00001
			BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = jmax
		}
		
#####################################################################################
# END OF OPTIONS BLOCK                     		                                    #
# the script should run completely from now on if everything above is set correctly #
#####################################################################################

	#checks the number of stochastic maps to be plotted
		if(n_stochastic_maps<n_maps_count_LTT){
			n_maps_count_LTT <- n_stochastic_maps
			print("correcting n_maps_count_LTT so it is equal to n_stochastic_maps")
		}

######################################
# create directories for the outputs #
######################################

	if(file.exists(paste(dir_name,"/state_by_state_by_event_tables", sep=""))==FALSE) dir.create(paste(dir_name,"/state_by_state_by_event_tables", sep=""))
	if(file.exists(paste(dir_name, "/tables", sep= ""))==FALSE) dir.create(paste(dir_name, "/tables", sep= ""))
	if(sample_trees & (file.exists(paste(dir_name, "/trees", sep= ""))==FALSE)) dir.create(paste(dir_name, "/trees", sep= ""))
	if(run_BGB){
		if(file.exists(paste(dir_name, "/BGB_figures", sep= ""))==FALSE) dir.create(paste(dir_name, "/BGB_figures", sep= ""))
		if(file.exists(paste(dir_name, "/BGB_results", sep= ""))==FALSE) dir.create(paste(dir_name, "/BGB_results", sep= ""))
	}
	if(run_BSM){
		if(file.exists(paste(dir_name, "/BSM_inputs", sep= ""))==FALSE) dir.create(paste(dir_name, "/BSM_inputs", sep= ""))
		if(file.exists(paste(dir_name, "/BSM_results", sep= ""))==FALSE) dir.create(paste(dir_name, "/BSM_results", sep= ""))
	}
	if(plot_BSM) if(file.exists(paste(dir_name, "/BSM_figures", sep= ""))==FALSE) dir.create(paste(dir_name, "/BSM_figures", sep= ""))
	
##################################################
# sample N trees from the posterior distribution #
##################################################
	
	#sometimes keep.tip produces topologies with weird branch lengths that prevent the script from running
	#if you get this error, change pruning_method_keep to FALSE
		pruning_method_keep <- FALSE
	
	# samples the trees from a posterior distribution
	# do NOT sample trees again if you plan to load previous results of BioGeoBEARS estimates!
	if (sample_trees==TRUE & run_BGB == TRUE & checkpoint_from_previous_run != TRUE) {
		tic()
		posterior_trees <- read.nexus(posterior_trees_fn)
		n_polytomous <- 0
			if(n_trees == 1) { # if this is a single target tree...
				posterior_trees <- ladderize(posterior_trees, right=TRUE)
				if(pruning_method_keep){
					pruned_trees <- keep.tip(posterior_trees, tip=species)
				} else{
					species_to_drop <- setdiff(posterior_trees$tip.label, species)
					pruned_trees <- drop.tip(posterior_trees, tip=species_to_drop)
				}
				if(is.binary(pruned_trees)){} else{
					#resolves polytomies and zero-length branches
						pruned_trees <- multi2di(pruned_trees, random = TRUE)
						pruned_trees$edge.length[pruned_trees$edge.length==0]<-1e-8
						pruned_trees[[i]] <- pruned_trees
						n_polytomous <- n_polytomous+1
				}
				write.tree(pruned_trees, file = paste(dir_name, "/trees/tree",1,".nex", sep= ""))
				tMax_allreps = as.integer(max(node.depth.edgelength(pruned_trees))) #gets age of the tree
				
			} else { # if this is a distribution of trees...
				class(posterior_trees) <- "multiPhylo"
				tMax_allreps = 0
				all_tree_ages = vector() #creates vector to store all tree ages
				#samples N trees from the posterior distribution
					posterior_trees <- sample(posterior_trees, size = n_trees)
				posterior_trees <-lapply(posterior_trees, ladderize, right=TRUE) #ladderizes trees
				#write.nexus(posterior_trees, file = paste(n_trees,"_posterior_trees.nex", sep="")) #in case you want the unpruned trees
				if(pruning_method_keep){
					pruned_trees <- lapply(posterior_trees, keep.tip, tip=species) #drops unwanted tips (outgroups, etc)
				} else{
					species_to_drop <- setdiff(posterior_trees[[1]]$tip.label, species)
					pruned_trees <- lapply(posterior_trees, drop.tip, tip=species_to_drop) #drops unwanted tips (outgroups, etc)
				}
				class(pruned_trees) <- "multiPhylo"
				
				for(i in 1:length(pruned_trees)) {
					tree <- pruned_trees[[i]]
						if(is.binary(tree)){} else{
							#resolves polytomies and zero-length branches
								tree <- multi2di(tree, random = TRUE)
								tree$edge.length[tree$edge.length==0]<-1e-8
								pruned_trees[[i]] <- tree
								n_polytomous <- n_polytomous+1
						}
					# imposes minimum branch lengths for BioGeoBEARS
					# this might be necessary in some trees with sampled fossils where the fossil edge lengths are too short
					# in these cases you get the error "FATAL ERROR in check_BioGeoBEARS_run(): the input tree has branchlengths <= 0, at these nodes:..."
					# this might fix this but is not necessary in most cases
						if(fossil_tips){
							tree <- impose_min_brlen(tree)
							height <- round(max(nodeHeights(tree)),1)
							for (a in 1:length(species)){
								if((round(nodeheight(tree, a),1))==height) next #checks if tip is a fossil
								if(tree$edge.length[which(tree$edge[,2]==a)] < 0.1) tree$edge.length[which(tree$edge[,2]==a)] <- 0.1 #if its branch length is smaller than 0.1, correct this for BioGeoBEARS
								print(paste("Corrected branch length of taxon ",a," in tree ",i,sep=""))
							}
						}
					#writes trees to individual files	
						write.tree(tree, file = paste(dir_name, "/trees/tree",i,".nex", sep= ""))
					#gets maximum age of oldest tree
						tMax = as.integer(max(node.depth.edgelength(tree)))
						all_tree_ages <- c(all_tree_ages, tMax)
						if (tMax_allreps < tMax) tMax_allreps <- tMax
				}
			}
			
		print(paste("Polytomous trees forcibly resolved:",n_polytomous)) #checks the number of trees forcibly resolved, usually this is 0
		write.nexus(pruned_trees, file = paste(dir_name, "/trees/",n_trees,"_pruned_trees.nex", sep= "")) #all the trees in a single file
		toc()
	} else {tMax_allreps = 0}
	
##################################
# create tables to count events  #
##################################

	# get the list of state names
		statenames = areas_list_to_states_list_new(getareas_from_tipranges_object(tipranges))
		statelabels = vector(length=num_states)
		# get only labels for possible ranges (disconsiders ranges including more areas than max_range_size)
			for (a in 1:num_states){statelabels[[a]] <- paste(statenames[[a]], collapse="")}
		
	# settings for plotting results of the stochastic maps	
		states_list_0based = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=TRUE)
		colors_list_for_states = get_colors_for_states_list_0based(areanames=areas, states_list_0based=states_list_0based, max_range_size=max_range_size, plot_null_range=TRUE)		

	# gets the number of edges
		trfn <- paste(dir_name, "/trees/tree",1,".nex", sep= "")
		tr <- read.tree(trfn)
		ingroup <- drop.tip(tr, tip = outgroup_taxa)
		n_edges <- length(ingroup$edge[,1])
		if(length(tr$tip.label) > length(ingroup$tip.label)) n_edges <- n_edges+1 # if outgroups are skipped, we still want the ancestral node of the ingroup to be counted
		
	# make a list of vectors to place the transitions by replicate and edge
	# this prints the identify of the ancestor and descendant nodes,
	# the distribution range of each of them,
	# the transition ("A -> AB"), if any,
	# and the most likely time of the transition, weighted by probabilities at each node
		ColumnLenght =n_trees*(n_edges)*n_stochastic_maps
		BSM_transitions = list(
			replicate= c(rep(0, ColumnLenght)),
			BSMreplicate= c(rep(0, ColumnLenght)), 
			nodeStart= c(rep(0, ColumnLenght)), 
			nodeEnd= c(rep(0, ColumnLenght)),
			nodeEndLabel= c(rep("", ColumnLenght)),
			ageStart= c(rep(0, ColumnLenght)), 
			ageEnd= c(rep(0, ColumnLenght)),
			rangeAncestor= c(rep("", ColumnLenght)),
			rangeStart= c(rep("", ColumnLenght)),
			rangeEnd= c(rep("", ColumnLenght)),
			transitionCladog= c(rep("", ColumnLenght)),
			transitionCladogType= c(rep("", ColumnLenght)),
			transitionAnag= c(rep("", ColumnLenght)),
			transitionAnagType= c(rep("", ColumnLenght)),
			transitionTime= c(rep(0, ColumnLenght)),
			nodeOfInterest= c(rep("", ColumnLenght)),
			BGBeventTypeClado= c(rep("", ColumnLenght)),
			BGBeventTypeAnag= c(rep("", ColumnLenght))
		)
		counter = 1 # this counter will be used to populate the BSM_transitions table
	
	# make a matrix to store lineages per area per time slice
		CountStatesByTime_all_reps = matrix(NA, ncol = (length(statelabels)+4+num_areas), nrow = 0)
		names_single_areas <- vector(length = num_areas)
		for(i in 1:num_areas){names_single_areas[[i]] <- paste("SUM_",statelabels[[i+1]],sep="")}
		colnames(CountStatesByTime_all_reps) = c("TimePoint","replication","BSM",statelabels, "LTT", names_single_areas)
		
	# event names
		n_insituesp <- "in_situ_speciation"
		n_foundesp <- "founder_event_speciation"
		n_vicsubsym <- "vicariance_or_subset_sympatry" # we will distinguish these events later by looking at both descendants
		n_widesym_subsym <- "wide_sympatry_or_subset_sympatry" # we will distinguish these events later by looking at both descendants
		n_vic <- "allopatry"
		n_subsym <- "subset_sympatry"
		n_dispersal <- "dispersal"
		n_extinction <- "extinction"
		n_disp_ext <- "dispersal_and_extinction"
		n_impos_event <- "X"
		no_event <- "–"
	
	# make a matrix to store general statistics of the trees
		transition_types <- c(n_insituesp, n_foundesp, n_vic, n_subsym, n_dispersal, n_extinction)
		transition_types_mean <-lapply(transition_types, paste, "_mean", sep="")
		transition_types_SD <-lapply(transition_types, paste, "_SD", sep="")
		ReplicationParameters = matrix(0, ncol = 10+(2*length(transition_types)), nrow = n_trees)
		colnames(ReplicationParameters) <- c("replication", "log-likelihood", "AICc", "k", "n_taxa", "d", "e", "j", "w", "treeHeight", transition_types_mean, transition_types_SD)
	
	# make matrices to count the number of transitions between states
		# template for event tables of each replication
		   state_by_state_transitions_template = matrix(0, ncol = length(statelabels), nrow = length(statelabels))
	       colnames(state_by_state_transitions_template) = statelabels
	       rownames(state_by_state_transitions_template) = statelabels 
		   
	# a function by Liam Revell to get node numbers for descendants of a particular node. Thanks Liam!
		getDescendants<-function(tree,node,curr=NULL){
			if(is.null(curr)) curr<-vector()
			daughters<-tree$edge[which(tree$edge[,1]==node),2]
			curr<-c(curr,daughters)
			w<-which(daughters>=length(tree$tip))
			if(length(w)>0) for(i in 1:length(w)) 
				curr<-getDescendants(tree,daughters[w[i]],curr)
			return(curr)
		}

###################################################################################
# create shorter modifier matrices (with dispersal multipliers, area allowed, etc #
###################################################################################

	#this will be important in case some of the trees have total age small enough to cross boundaries of the time slices
		if (time_stratified) {
			if(file.exists(paste(dir_name, "/modifier_matrices", sep= ""))==FALSE) dir.create(paste(dir_name, "/modifier_matrices", sep= ""))
			#gets boundaries for time slices
				time_slices <- scan(time_stratified_fn)
				time_slice_control <- max(time_slices)
			#get matrices to be modified
				modifier_matrices <- c(dispersal_matrix_fn, area_adjacency_fn, distance_matrix_fn, areas_allowed_fn)
				modifier_matrices <- modifier_matrices[modifier_matrices != ""]
			#size of each matrix per time slice
				n_lines_matrix <- 2+num_areas
			#prepares modified matrices by successively removing the last time slice
				for (t in 0:(length(time_slices)-1)){	
					for (i in 1:length(modifier_matrices)) {
						write(c(readLines(modifier_matrices[i], n=(1+t)*n_lines_matrix),"END"),paste(dir_name, "/modifier_matrices/",t+1,"_",modifier_matrices[i], sep= ""))
					}
					reduced_time_slices <- head(time_slices, t+1)
					cat(reduced_time_slices, file=paste(dir_name, "/modifier_matrices/",t+1,"_",time_stratified_fn,sep=""), sep = "\n")
				}
		} else {
			time_slice_control <- tMax_allreps+1	
		}
		
#####################################################################
# create tables to classify transition events among possible ranges #
#####################################################################	
		
	# creates reference tables to classify types of transitions
		state_by_state_Cladogenetic_events <- state_by_state_transitions_template
		state_by_state_Anagenetic_events <- state_by_state_transitions_template
		
	# fills reference table: gets names of source and sink ranges, and classifies the events
	# possible events: in situ speciation, founder event speciation, vicaricance/subset sympatry, dispersal, extinction, dispersal+extinction
	# no event (-), unknown event (?), should-be-impossible event (X)
		#loops through all rows
		for (o in 1:length(statelabels)){
			# gets source range
				range_source <- rownames(state_by_state_Cladogenetic_events)[[o]]
			#loops through all columns
			for (k in 1:length(statelabels)){
				#gets new range
				range_sink <- colnames(state_by_state_Cladogenetic_events)[[k]]
				
				#no change in range:
					if(range_source == range_sink) {
						if (nchar(range_source)==1){
							state_by_state_Cladogenetic_events[range_source, range_sink] <- n_insituesp
						} else {
							state_by_state_Cladogenetic_events[range_source, range_sink] <- n_widesym_subsym #this can be either in-situ speciation in a range with more than 1 area, or the larger range in subset sympatry
						}
						state_by_state_Anagenetic_events[range_source, range_sink] <- no_event
					} 
					
				# if new range is a larger than source range:	
					if (nchar(range_sink) > nchar(range_source)) {
						state_by_state_Cladogenetic_events[range_source, range_sink] <- n_impos_event
						state_by_state_Anagenetic_events[range_source, range_sink] <- n_dispersal
						# checks if source range contains an area not present in new range (at least 1 extinction)
							for (a in 1:nchar(range_source)){
								if(str_contains(range_sink,(substring(range_source,a,a)))==FALSE){state_by_state_Anagenetic_events[range_source, range_sink] <- n_disp_ext}
							}
					}
				
				# if new range is a smaller than source range:	
					if (nchar(range_sink) < nchar(range_source)) {
						state_by_state_Cladogenetic_events[range_source, range_sink] <- n_vicsubsym #this can be either vicariance or subset sympatry
						state_by_state_Anagenetic_events[range_source, range_sink] <- n_extinction
						# checks if new range contains an area not present in source range (at least one dispersal)
							for (a in 1:nchar(range_sink)){
								if(str_contains(range_source,(substring(range_sink,a,a)))==FALSE){
									state_by_state_Cladogenetic_events[range_source, range_sink] <- n_impos_event
									state_by_state_Anagenetic_events[range_source, range_sink] <- n_disp_ext
								}
							}
					}				

				# if new range is of same size, but different areas than source:
					if (nchar(range_sink) == nchar(range_source) & range_source != range_sink) {
						state_by_state_Cladogenetic_events[range_source, range_sink] <- n_impos_event
						state_by_state_Anagenetic_events[range_source, range_sink] <- n_disp_ext
					}
							
				# if new range is a single area not contained in source area:
					if (nchar(range_sink) == 1 & str_contains(range_source,range_sink)==FALSE) {
						state_by_state_Cladogenetic_events[range_source, range_sink] <- n_foundesp
					}				
				}				
			}
		
		# writes tables
		# perhaps a good idea to check them!
			write.csv(state_by_state_Cladogenetic_events, file = paste(dir_name,"/tables/Reference_table_cladogenetic_events.csv", sep=""))
			write.csv(state_by_state_Anagenetic_events, paste(dir_name,"/tables/Reference_table_anagenetic_events.csv", sep=""))
	
##########################################################################
# loop over sampled trees to estimate ancestral ranges and count events  #
##########################################################################
	
	# only proceeds with analysis if everything is OK
		if (DIVALIKE == TRUE & BAYAREALIKE ==TRUE) {
			stop("Set either DIVALIKE or BAYAREALIKE to TRUE, but not both at the same time!")
		} else {
		if (time_slice_control<tMax_allreps) {
			stop("At least one of your trees is older than your oldest time slice; please check your input files.")
		} else{
		#starts the loop
		for (r in 1:n_trees) {
			tic()
			# reads one of the trees of the posterior distribution and loads it to BGB object
				trfn <- paste(dir_name, "/trees/tree",r,".nex", sep= "")
				tr <- read.tree(trfn)
				ingroup <- drop.tip(tr, tip = outgroup_taxa)
			#get the height of the tree; defines how many time slices it spans
				current_height <- max(nodeHeights(tr))
				if(time_stratified) how_many_time_slices <- min(which(time_slices > current_height))
			
			# make a list of nodes of your ingroup
			# use this when you want to include outgroups for estimating ancestral ranges,
			# but not for counting transitions and lineages through time
			# check carefully that ingroupNode is actually the node containing all of your ingroup!
				root_node_ingroup <- length(ingroup$tip.label)+1
				ingroupNode = matchNodes(ingroup, tr, method=c("descendants"))[which(matchNodes(ingroup, tr, method=c("descendants"))[,1]==root_node_ingroup),2]
				ListIngroupNodes = getDescendants(tr, node=ingroupNode,curr=NULL)
				# get the node immediately ancestor to ingroupNode -- to get biogeographic transitions in this node as well
					first_node <- tr$edge[1,which(tr$edge[,2]==min(ListIngroupNodes[ListIngroupNodes > length(ingroup$tip.label)]))]
					if (length(tr$tip.label) > length(ingroup$tip.label)) {
						ListIngroupNodes <- c(ListIngroupNodes, first_node)
					}
							
		#############################################################
		# estimate ancestral ranges for each of the posterior trees #
		#############################################################

			# checkpoint for resuming analysis
				if(checkpoint_from_previous_run) {
					run_BGB = FALSE
					run_BSM = FALSE
					if(number_of_trees_run_previously < r) {
						run_BGB = TRUE
						run_BSM = TRUE
					}	
				}
			
			# estimates ancestral ranges using BioGeoBEARS and plots results for each tree
			if (run_BGB) { #runs analysis
				#loads and sections the tree
					BioGeoBEARS_run_object$trfn = trfn # Give BioGeoBEARS the location of the phylogeny Newick file
			# adds files with time slices, dispersal matrices, etc. if needed, respecting the tree height of current tree
					if (time_stratified) {
						if(how_many_time_slices>1) BioGeoBEARS_run_object$timesfn = paste(dir_name, "/modifier_matrices/",how_many_time_slices,"_",time_stratified_fn,sep="")
						if (areas_allowed) BioGeoBEARS_run_object$areas_allowed_fn = paste(dir_name, "/modifier_matrices/",how_many_time_slices,"_",areas_allowed_fn,sep="")
						if (area_adjacency) BioGeoBEARS_run_object$areas_adjacency_fn = paste(dir_name, "/modifier_matrices/",how_many_time_slices,"_",area_adjacency_fn,sep="")
						if (distance_matrix) BioGeoBEARS_run_object$distsfn = paste(dir_name, "/modifier_matrices/",how_many_time_slices,"_",distance_matrix_fn,sep="")
						if (dispersal_matrix) BioGeoBEARS_run_object$dispersal_multipliers_fn = paste(dir_name, "/modifier_matrices/",how_many_time_slices,"_",dispersal_matrix_fn,sep="")
						BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
						if(how_many_time_slices>1) BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE, cut_fossils=FALSE)
					} else if(files_to_read) BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
				#estimates ancestral ranges
					BGB_res = bears_optim_run(BioGeoBEARS_run_object)
				# saves results object
					save(BGB_res, file=paste(dir_name, "/BGB_results/",filename,"_rep",r,".Rdata", sep =""))
				# plots results as pie/text in pdf/jpeg formats
					while (!is.null(dev.list())) {dev.off()} #cleans plotting area
					pdf(paste(dir_name, "/BGB_figures/text_",filename,"_rep", r,".pdf", sep =""), width=pdf_w, height=pdf_h)
						plot_BioGeoBEARS_results(BGB_res, filename, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
					while (!is.null(dev.list())) {dev.off()}
					pdf(paste(dir_name, "/BGB_figures/pie_",filename,"_rep", r,".pdf", sep =""), width=pdf_w, height=pdf_h)
						plot_BioGeoBEARS_results(BGB_res, filename, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
					while (!is.null(dev.list())) {dev.off()}
					jpeg(paste(dir_name, "/BGB_figures/text_",filename,"_rep", r,".jpg", sep =""), w= pdf_w*100, h= pdf_h*100, pointsize = 18, quality = 90)
						plot_BioGeoBEARS_results(BGB_res, filename, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
					while (!is.null(dev.list())) {dev.off()}
					jpeg(paste(dir_name, "/BGB_figures/pie_",filename,"_rep", r,".jpg", sep =""), w= pdf_w*100, h= pdf_h*100, pointsize = 18, quality = 90)
						plot_BioGeoBEARS_results(BGB_res, filename, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
					while (!is.null(dev.list())) {dev.off()}
				# plots node and tip numbers just in case we have to check something
					pdf(paste(dir_name, "/BGB_figures/nodenumbers_",filename,"_rep", r,".pdf", sep =""), width=pdf_w, height=pdf_h)
						plot(tr)
						nodelabels()
						tiplabels()
						axisPhylo()
					while (!is.null(dev.list())) {dev.off()}	
			} else {  #or loads previous results
				load(paste(dir_name, "/BGB_results/", filename,"_rep",r,".Rdata", sep =""))
			}
			
			# fetch the results for counting events later
				# take the node probabilities from results of analyses
				#probabilities at the end of the branch
					nodeprobs_top_table <- BGB_res$ML_marginal_prob_each_state_at_branch_top_AT_node
				#probabilities at the bottom of the branch, right after cladogenesis
					nodeprobs_bottom_table <- BGB_res$ML_marginal_prob_each_state_at_branch_bottom_below_node
									
			# compose a list of edges, their nodes and ages
			# tr$edge comes as NodeAncestor NodeDescendent (or nodeStart nodeEnd of every edge)
				nodeAges = round((max(node.depth.edgelength(tr)) - node.depth.edgelength(tr)),digits = 6)
				edgeAges = list(start=nodeAges[tr$edge[,1]], end=nodeAges[tr$edge[,2]])
				
		########################################
		# run biogeographic stochastic mapping #
		########################################
		
			if (run_BSM) {
				#creates directory for storing results
					dir.create(paste(dir_name, "/BSM_results/BSM_rep",r, sep= ""))
				#get inputs and runs BSM
					stochastic_mapping_inputs_list = get_inputs_for_stochastic_mapping(res=BGB_res)
					BSM_inputs_fn = paste(dir_name,"/BSM_inputs/BSM_inputs_file",r,".Rdata", sep="")
					save(stochastic_mapping_inputs_list, file=BSM_inputs_fn)
					BSM_output = runBSM(BGB_res, stochastic_mapping_inputs_list=stochastic_mapping_inputs_list, maxnum_maps_to_try=n_stochastic_maps*10, nummaps_goal=n_stochastic_maps, maxtries_per_branch=40000, save_after_every_try=TRUE, savedir=paste(dir_name, "/BSM_results/BSM_rep",r, sep= ""), seedval=12345, wait_before_save=0.01)
					RES_clado_events_tables = BSM_output$RES_clado_events_tables
					RES_ana_events_tables = BSM_output$RES_ana_events_tables
			} else {
				#loads previous results
					load(file=paste(dir_name, "/BSM_results/BSM_rep",r,"/RES_clado_events_tables.Rdata", sep= ""))
					load(file=paste(dir_name, "/BSM_results/BSM_rep",r,"/RES_ana_events_tables.Rdata", sep= ""))
			}
			
			#plots stochastic maps
				if(plot_BSM){
					#loops through all BSM and plots results
						for (i in 1:n_stochastic_maps){
							map_fn=paste(dir_name, "/BSM_figures/BSM_rep",r,"_map",i,sep="")
							if(file.exists(map_fn)) next
							master_table_cladogenetic_events = RES_clado_events_tables[[i]]
							resmod = stochastic_map_states_into_res(res=BGB_res, master_table_cladogenetic_events=master_table_cladogenetic_events, stratified=time_stratified)
							while (!is.null(dev.list())) {dev.off()} #cleans plotting area
							pdf(file=paste(map_fn, ".pdf", sep=""), width=pdf_w, height=pdf_h)
								plot_BioGeoBEARS_results(results_object=resmod, analysis_titletxt="Stochastic map", addl_params=list("j"), label.offset=0.5, plotwhat="text", cornercoords_loc=scriptdir, root.edge=TRUE, colors_list_for_states=colors_list_for_states, skiptree=FALSE, show.tip.label=TRUE)
								paint_stochastic_map_branches(res=resmod, master_table_cladogenetic_events=master_table_cladogenetic_events, colors_list_for_states=colors_list_for_states, lwd=5, lty=par("lty"), root.edge=TRUE, stratified=time_stratified)
								plot_BioGeoBEARS_results(results_object=resmod, analysis_titletxt="Stochastic map", addl_params=list("j"), plotwhat="text", cornercoords_loc=scriptdir, root.edge=TRUE, colors_list_for_states=colors_list_for_states, skiptree=TRUE, show.tip.label=TRUE)
							while (!is.null(dev.list())) {dev.off()}
							jpeg(file=paste(map_fn, ".jpg", sep=""), w= pdf_w*100, h= pdf_h*100, pointsize = 18, quality = 90)
								plot_BioGeoBEARS_results(results_object=resmod, analysis_titletxt="Stochastic map", addl_params=list("j"), label.offset=0.5, plotwhat="text", cornercoords_loc=scriptdir, root.edge=TRUE, colors_list_for_states=colors_list_for_states, skiptree=FALSE, show.tip.label=TRUE)
								paint_stochastic_map_branches(res=resmod, master_table_cladogenetic_events=master_table_cladogenetic_events, colors_list_for_states=colors_list_for_states, lwd=5, lty=par("lty"), root.edge=TRUE, stratified=time_stratified)
								plot_BioGeoBEARS_results(results_object=resmod, analysis_titletxt="Stochastic map", addl_params=list("j"), plotwhat="text", cornercoords_loc=scriptdir, root.edge=TRUE, colors_list_for_states=colors_list_for_states, skiptree=TRUE, show.tip.label=TRUE)
							while (!is.null(dev.list())) {dev.off()}
						}
				}
					
		####################################
		# get the most probable root state #
		####################################
			
			# make a matrix to place the probabilities at the root for each replication
			# this can be easily modified to get probabilities for other nodes of interest by using getMRCA
			# e.g. get_probs_for_node <- getMRCA(tr, tip = c("species1", "species2")
				get_probs_for_node <- length(species)+1
				if(r==1){
					RootState = matrix(NA, ncol = (length(nodeprobs_top_table[1,])+2), nrow = 0)
					colnames(RootState) = c("replication",head(statelabels, length(nodeprobs_top_table[1,])),"most_likely")
				}
				most_likely_state <- statelabels[which.max(nodeprobs_top_table[get_probs_for_node,])]
				RootState <-rbind(RootState,c(r, nodeprobs_top_table[get_probs_for_node,], most_likely_state))
							
		##############################################################
		# get the transitions and their ages from stochastic mapping #
		##############################################################
		
			#gets the node number for clades of interest for this particular tree
				node_interest_n <- sapply(nodes_of_interest, getMRCA, phy = tr)
		
			#loops through edges
				for(b in 1:n_stochastic_maps){
					RES_clado_events_tables[[b]] <- as.data.frame(RES_clado_events_tables[[b]])
					parsed_nodes <- c()
					for(e in 1:length(RES_clado_events_tables[[b]]$node)){
					  if((RES_clado_events_tables[[b]]$node[e] %in% ListIngroupNodes)==FALSE) next
					  unique_node <- RES_clado_events_tables[[b]]$node[e]
					  if((unique_node %in% parsed_nodes)==FALSE){ #some (not all!) BSM tables have duplicate rows for each node where their branches cross boundaries between time slices
						parsed_nodes <- c(parsed_nodes, unique_node) # only gets the data for each node once
						node_label <- RES_clado_events_tables[[b]]$label[e] # gets the name of the current node
						if (time_stratified) {oldest_stratum <- max(RES_clado_events_tables[[b]]$stratum[which(RES_clado_events_tables[[b]]$label==node_label)])} #checks the deepest stratum to where the branch leading to the current node extends
						#current tree
							BSM_transitions$replicate[counter] <- r
						#current stochastic map
							BSM_transitions$BSMreplicate[counter] <- b
						#ancestral of current node
							BSM_transitions$nodeStart[counter] <- RES_clado_events_tables[[b]]$ancestor[e]
						#current node
							BSM_transitions$nodeEnd[counter] <- RES_clado_events_tables[[b]]$node[e]
						#name of current node	
							BSM_transitions$nodeEndLabel[counter] <- node_label
						#age of current node
							BSM_transitions$ageEnd[counter] <- round(RES_clado_events_tables[[b]]$time_bp[e],n_digits_age)
						#range of current node
							BSM_transitions$rangeEnd[counter] <- statelabels[RES_clado_events_tables[[b]]$sampled_states_AT_nodes[e]]
						#skips root when getting ancestral states and start age of nodes
							if(RES_clado_events_tables[[b]]$node[e]==length(species)+1){ #if this is the root node...
								BSM_transitions$rangeAncestor[counter] <- "–"
								BSM_transitions$rangeStart[counter] <- "–"
								BSM_transitions$ageStart[counter] <- "–"
								BSM_transitions$transitionCladog[counter] <- "–"
								BSM_transitions$transitionCladogType[counter] <- "–"
								BSM_transitions$transitionAnag[counter] <- "–"
								BSM_transitions$transitionAnagType[counter] <- "–"
								BSM_transitions$transitionTime[counter] <- "–"
							} else {
								#age of ancestral node
									BSM_transitions$ageStart[counter] <- round(RES_clado_events_tables[[b]]$time_bp[which(RES_clado_events_tables[[b]]$node==RES_clado_events_tables[[b]]$ancestor[e])][1],n_digits_age)
								#range of ancestral node
									BSM_transitions$rangeAncestor[counter] <- statelabels[RES_clado_events_tables[[b]]$sampled_states_AT_nodes[which(RES_clado_events_tables[[b]]$node==RES_clado_events_tables[[b]]$ancestor[[e]])]][1]
								#range of descendant immediately after cladogenesis
									if (time_stratified) {
										BSM_transitions$rangeStart[counter] <- statelabels[RES_clado_events_tables[[b]]$sampled_states_AT_brbots[which(RES_clado_events_tables[[b]]$label==node_label & RES_clado_events_tables[[b]]$stratum==oldest_stratum)]] #gets the starting range of the true bottom of the branch (in the oldest stratum where it is)
									} else {
										BSM_transitions$rangeStart[counter] <- statelabels[RES_clado_events_tables[[b]]$sampled_states_AT_brbots[e]]
									}
								if(BSM_transitions$rangeAncestor[counter] == BSM_transitions$rangeStart[counter]) {
									BSM_transitions$transitionCladog[counter] = paste(BSM_transitions$rangeAncestor[counter],"->",BSM_transitions$rangeAncestor[counter])
								} else {
									# which range from which range - cladogenetic event
										BSM_transitions$transitionCladog[counter] = paste(BSM_transitions$rangeAncestor[counter], "->", BSM_transitions$rangeStart[counter])
								}
								#cladogenetic transition type
									BSM_transitions$transitionCladogType[counter] = state_by_state_Cladogenetic_events[BSM_transitions$rangeAncestor[counter],BSM_transitions$rangeStart[counter]]
								#cladogenetic transition type according to BGB	
									BSM_transitions$BGBeventTypeClado[counter] = RES_clado_events_tables[[b]]$clado_event_type[e]
								if(BSM_transitions$rangeStart[counter] == BSM_transitions$rangeEnd[counter]){
									BSM_transitions$transitionAnag[counter] = "–"
									BSM_transitions$transitionAnagType[counter] = "–"
									BSM_transitions$transitionTime[counter] = "–"
								} else {
									# which range from which range - anagenetic event
										BSM_transitions$transitionAnag[counter] = paste(BSM_transitions$rangeStart[counter], "->", BSM_transitions$rangeEnd[counter])
									#anagenetic transition type
										BSM_transitions$transitionAnagType[counter] = state_by_state_Anagenetic_events[BSM_transitions$rangeStart[counter],BSM_transitions$rangeEnd[counter]]
									#anagenetic transition type according to BGB
										if(length(which(RES_ana_events_tables[[b]]$node==unique_node))>0) {BSM_transitions$BGBeventTypeAnag[counter] <- RES_ana_events_tables[[b]]$event_type[which(RES_ana_events_tables[[b]]$node==unique_node)][1]}
									#time of anagenetic transition
										if(length(which(RES_ana_events_tables[[b]]$node==unique_node))>0) {BSM_transitions$transitionTime[counter] <- round(RES_ana_events_tables[[b]]$abs_event_time[which(RES_ana_events_tables[[b]]$node==unique_node)][1],n_digits_age)} else {BSM_transitions$transitionTime[counter] ="–"}
								} 
								if(length(nodes_of_interest)>0){								
									for (x in 1:length(node_interest_n)){
										# check if node of interest is monophyletic
										# gets all descendants of node of interest in current tree; if they include additional terminals, it is not monophyletic and the node is not marked in the current replication
											get_terminals <- getDescendants(tr, node_interest_n[x],curr=NULL)
											check_interest_monophyletic <- sort(tr$tip.label[get_terminals[get_terminals<length(tr$tip.label)]])
											if(prod(check_interest_monophyletic %in% nodes_of_interest[[x]])==1){
												#marks corresponding rows in table with name of nome
												if (node_interest_n[x] == RES_clado_events_tables[[b]]$ancestor[e]) {BSM_transitions$nodeOfInterest[counter] = paste("descendant-",names(node_interest_n[x]), sep="")}
												if (node_interest_n[x] == RES_clado_events_tables[[b]]$node[e]) {BSM_transitions$nodeOfInterest[counter] = paste("ancestor-",names(node_interest_n[x]), sep="")}
											}	
									}
								}
							}
						counter = counter+ 1	
					}
				  }
				}	
			
		########################################
		# count lineages in areas through time #
		########################################
		
			if (plot_LTT) {
				#make a time vector to store the time series 
					tMax = as.integer(max(node.depth.edgelength(tr)))
					if (tMax_allreps < tMax) tMax_allreps <- tMax # store the maximum tree height for all replications; this is used to plot results later
					TimePoint <- seq(from=0, to= tMax, by= time_period_size) 
					
				# make a matrix to place the count of lineages by time and states
					CountStatesByTime = matrix(NA, ncol = (length(statelabels)+4+num_areas), nrow = length(TimePoint))
	
				# count lineages through time from stochastic maps
						count_time_transitions <- BSM_transitions
						names(count_time_transitions)[which(names(count_time_transitions)=="transitionTime")] <- "transitionTimeRandom"
				
				# loop over time periods and counts lineages in that time period
				for (s in 1:n_maps_count_LTT) {
					
					# only look in rows of the transitions table that correspond to the current replication (or replication+BSM)
						start_transitions_this_rep <- min(which(BSM_transitions$replicate==r & BSM_transitions$BSMreplicate==s))
						end_transitions_this_rep <- min(which(BSM_transitions$replicate==r & BSM_transitions$BSMreplicate==s))+(n_edges)-1
					
					#for each time point...
					for(t in 1:length(TimePoint)){
						CountStates = c(rep(0,(length(statelabels)))) # clean the count vector
						#...loops through edges
							for(e in start_transitions_this_rep:end_transitions_this_rep){
								# if there is no anagenetic transition, sums +1 to current range for all time period between start and end of the branch
									if (count_time_transitions$rangeStart[e]==count_time_transitions$rangeEnd[e]){
									if (TimePoint[t] < count_time_transitions$ageStart[e] & TimePoint[t] >= count_time_transitions$ageEnd[e]){
											CountStates[which(statelabels==count_time_transitions$rangeEnd[e])] <- CountStates[which(statelabels==count_time_transitions$rangeEnd[e])]+1
										}
									} else {
									# if there has been an anagenetic transition...	
										if(count_time_transitions$transitionAnagType[e]!="–" & count_time_transitions$transitionTimeRandom[e]=="–") count_time_transitions$transitionTimeRandom[e] <- round(runif(1,count_time_transitions$ageEnd[e], count_time_transitions$ageStart[e]),n_digits_age) #some BSM results are lacking transition times -- get it randomly from scratch in these cases
										# ...sums +1 to start range to time slices before transition...
											if (TimePoint[t] < count_time_transitions$ageStart[e] & TimePoint[t] >= as.numeric(count_time_transitions$transitionTimeRandom[e])){
												CountStates[which(statelabels==count_time_transitions$rangeStart[e])] <- CountStates[which(statelabels==count_time_transitions$rangeStart[e])]+1
											}
										# ...and sums +1 to end range to time slices after transition
											if (TimePoint[t] < as.numeric(count_time_transitions$transitionTimeRandom[e]) & TimePoint[t] >= count_time_transitions$ageEnd[e]){
												CountStates[which(statelabels==count_time_transitions$rangeEnd[e])] <- CountStates[which(statelabels==count_time_transitions$rangeEnd[e])]+1
											}
									}
							}
							CountStatesByTime[t,]
							# sums the total number of species in each single area per time period -- including widespread species
								sum_single_areas <- vector(length= num_areas)
								for (y in 1:num_areas) {
									single_area <- statelabels[[y+1]]
									sum_single_area <- 0
									for (x in 1:length(statelabels)){
										if(str_contains(statelabels[[x]], single_area)){
											sum_single_area <- sum_single_area+as.numeric(CountStates[x])
										}
									}
									sum_single_areas[[y]]<-sum_single_area
								}
								
						BSM_n = s
						CountStatesByTime[t,] = c(TimePoint[t],r,BSM_n,CountStates, sum(CountStates), sum_single_areas)
					}
					
					# stores the counts from each rep in a table
						CountStatesByTime_all_reps <- rbind(CountStatesByTime_all_reps, CountStatesByTime)
				}
			}
			
	###########################################################################
	# calculates fit using AICc and stores parameters for current replication #
	###########################################################################			
			
			# calculates AICc
				n_parameters <- sum(BGB_res$inputs$BioGeoBEARS_model_object@params_table[,1]=="free")
				n_taxa <- length(species)
				LnL <- get_LnL_from_BioGeoBEARS_results_object(BGB_res)
				AICc <- aicc <- (-2*LnL)+((2*n_parameters)+(((2*n_parameters)*(n_parameters+1))/(n_taxa-n_parameters-1)))
			
			# stores parameters for each replication
				parameters_for_this_rep <- c(r,
				LnL,
				AICc,
				n_parameters,
				n_taxa,
				BGB_res$outputs@params_table["d","est"],
				BGB_res$outputs@params_table["e","est"],
				BGB_res$outputs@params_table["j","est"],
				BGB_res$outputs@params_table["w","est"],
				round(max(edgeAges$start),n_digits_age))
				for (k in 1:length(parameters_for_this_rep)){
					ReplicationParameters[r, k] <- parameters_for_this_rep[k]
				}
			print(paste("Completed replicate ",r,"/",n_trees,sep=""))
			toc() #prints total time spent in each replicate
		}
				
#################################################################################
# get the average number of transitions between ranges from the stochastic maps #
#################################################################################

	tic()
	print("Get the average number of transitions between ranges from the stochastic maps...")
	# get the transitions estimated by each BSM, divided by event type
			state_by_state_dispersal_BSM = state_by_state_transitions_template
			state_by_state_extinction_BSM = state_by_state_transitions_template
			state_by_state_vicar_BSM = state_by_state_transitions_template
			state_by_state_subset_sympatry_BSM = state_by_state_transitions_template
			state_by_state_insituspec_BSM = state_by_state_transitions_template
			state_by_state_founderevent_BSM = state_by_state_transitions_template
			state_by_state_dispersal_extinctions_BSM = state_by_state_transitions_template
		    state_by_state_transitions_average_BSM = state_by_state_transitions_template #all events
	
	# loops through BSM transitions table
	# stores transition events in all-ranges-by-all-ranges tables
			for (a in 1:(ColumnLenght-1)){
				if (BSM_transitions$nodeEnd[a] == first_node) next
				
				# gets ranges
					range_ancestor <- BSM_transitions$rangeAncestor[a]
					range_start <- BSM_transitions$rangeStart[a]
					range_end <- BSM_transitions$rangeEnd[a]
				
				# creates a table to check both descendants in case of vicariance/subset sympatry or wide sympatry/subset sympatry
				# this allows distinguishing between the two alternatives
					if(a == 1){
						check_subset_sympatry <- as.data.frame(matrix(0, ncol = 4, nrow = n_edges))
						colnames(check_subset_sympatry) <- c("nodeStart","CladogType", "rangeAncestor", "rangeStart")
						check_subset_sympatry[,1] <- BSM_transitions$nodeStart[a:(a+n_edges-1)]
						check_subset_sympatry[,2] <- BSM_transitions$transitionCladogType[a:(a+n_edges-1)]
						check_subset_sympatry[,3] <- BSM_transitions$rangeAncestor[a:(a+n_edges-1)]
						check_subset_sympatry[,4] <- BSM_transitions$rangeStart[a:(a+n_edges-1)]
					} else if ((BSM_transitions$BSMreplicate[a]!= BSM_transitions$BSMreplicate[a-1]) | (BSM_transitions$replicate[a]!= BSM_transitions$replicate[a-1])){
						check_subset_sympatry <- as.data.frame(matrix(0, ncol = 4, nrow = n_edges))
						colnames(check_subset_sympatry) <- c("nodeStart","CladogType", "rangeAncestor", "rangeStart")
						check_subset_sympatry[,1] <- BSM_transitions$nodeStart[a:(a+n_edges-1)]
						check_subset_sympatry[,2] <- BSM_transitions$transitionCladogType[a:(a+n_edges-1)]
						check_subset_sympatry[,3] <- BSM_transitions$rangeAncestor[a:(a+n_edges-1)]
						check_subset_sympatry[,4] <- BSM_transitions$rangeStart[a:(a+n_edges-1)]
					}
				
				# counts events by type of transition
				if(range_ancestor != "–" & range_start != "–" & range_end != "–") {
				#dispersal+extinction in the same branch
					if(BSM_transitions$transitionAnagType[a]==n_disp_ext) {state_by_state_dispersal_extinctions_BSM[range_start,range_end]=(state_by_state_dispersal_extinctions_BSM[range_start,range_end])+1}	
						
				#dispersal
					if(BSM_transitions$transitionAnagType[a]==n_dispersal) {state_by_state_dispersal_BSM[range_start,range_end]=(state_by_state_dispersal_BSM[range_start,range_end])+1}

				#extinction
					if (BSM_transitions$transitionAnagType[a] == n_extinction) {state_by_state_extinction_BSM[range_start,range_end]=(state_by_state_extinction_BSM[range_start,range_end])+1}							
						
				# distinguishes between vicariance and subset sympatry
					if (BSM_transitions$transitionCladogType[a] == n_vicsubsym) {
						# finds the two descendants of the same ancestor and their ranges
							both_descendants <- which(check_subset_sympatry$nodeStart==BSM_transitions$nodeStart[a])
							range1 <- check_subset_sympatry[both_descendants[1],4]
							range2 <- check_subset_sympatry[both_descendants[2],4]
						# checks if at least one of the ranges is the same as the ancestor and one consists of a single area
						# if so, subset sympatry
							if((nchar(range1) == 1 | nchar(range2) == 1) & (range1==range_ancestor | range2==range_ancestor)){
								state_by_state_subset_sympatry_BSM[range_ancestor,range_start]=(state_by_state_subset_sympatry_BSM[range_ancestor,range_start])+1
								BSM_transitions$transitionCladogType[a] <- n_subsym
						# otherwise, vicariance
							} else if (nchar(range1)<nchar(range_ancestor) & nchar(range2)<nchar(range_ancestor)) {
								state_by_state_vicar_BSM[range_ancestor,range_start]=(state_by_state_vicar_BSM[range_ancestor,range_start])+1
								BSM_transitions$transitionCladogType[a] <- n_vic
							}
					}
					
					# distinguishes between wide sympatry and subset sympatry
					if (BSM_transitions$transitionCladogType[a] == n_widesym_subsym) {
						# finds the two descendants of the same ancestor and their ranges
							both_descendants <- which(check_subset_sympatry$nodeStart==BSM_transitions$nodeStart[a])
							range1 <- check_subset_sympatry[both_descendants[1],4]
							range2 <- check_subset_sympatry[both_descendants[2],4]
						# checks if ranges of ancestor and both descendants are the same
						# if so, in situ speciation in wide sympatry
							if(range1 == range2 & range2 == range_ancestor){
								state_by_state_insituspec_BSM[range_ancestor,range_start]=(state_by_state_insituspec_BSM[range_ancestor,range_start])+1
								BSM_transitions$transitionCladogType[a] <- paste(n_insituesp, "_WIDE", sep="")
							} else {
								#get the smaller range
									if(nchar(range2)>nchar(range1)){
										smaller_range=range1
										larger_range=range2
									} else if(nchar(range1)>nchar(range2)){
										smaller_range=range2
										larger_range=range1
									}
								# if smaller range is contained in ancestral range, subset sympatry
									if (str_contains(range_ancestor, smaller_range)){
										state_by_state_subset_sympatry_BSM[range_ancestor,range_start]=(state_by_state_subset_sympatry_BSM[range_ancestor,range_start])+1
										BSM_transitions$transitionCladogType[a] <- n_subsym
								# if smaller range is a jump dispersal, and the larger range is the same as the ancestor range
								# transition is a wide in situ speciation
									} else if (nchar(smaller_range==1)& larger_range==range_ancestor) {
										state_by_state_insituspec_BSM[range_ancestor,range_start]=(state_by_state_insituspec_BSM[range_ancestor,range_start])+1
										BSM_transitions$transitionCladogType[a] <- paste(n_insituesp, "_WIDE", sep="")											
									}
							}
					}											
					
				#founder-event speciation
					if (BSM_transitions$transitionCladogType[a] == n_foundesp) {state_by_state_founderevent_BSM[range_ancestor,range_start]=(state_by_state_founderevent_BSM[range_ancestor,range_start])+1}
				
				#in situ speciation
					if (BSM_transitions$transitionCladogType[a] == n_insituesp) {state_by_state_insituspec_BSM[range_ancestor,range_start]=(state_by_state_insituspec_BSM[range_ancestor,range_start])+1}				
				}
			}
			
			#also counts when a branch has gone through dispersal AND extinction
				state_by_state_dispersal_BSM = state_by_state_dispersal_BSM + state_by_state_dispersal_extinctions_BSM
				state_by_state_extinction_BSM = state_by_state_extinction_BSM + state_by_state_dispersal_extinctions_BSM
			
			#makes a table with all transition types
				state_by_state_transitions_average_BSM  = state_by_state_dispersal_BSM + state_by_state_extinction_BSM + state_by_state_vicar_BSM + state_by_state_subset_sympatry_BSM + state_by_state_insituspec_BSM + state_by_state_founderevent_BSM - state_by_state_dispersal_extinctions_BSM
			
			#gets average matrices
				state_by_state_dispersal_average_BSM = round(state_by_state_dispersal_BSM,nchar(n_trees*n_stochastic_maps)-2)
				state_by_state_extinction_average_BSM = round(state_by_state_extinction_BSM,nchar(n_trees*n_stochastic_maps)-2)
				state_by_state_vicar_average_BSM = round(state_by_state_vicar_BSM,nchar(n_trees*n_stochastic_maps)-2)
				state_by_state_subset_sympatry_average_BSM = round(state_by_state_subset_sympatry_BSM,nchar(n_trees*n_stochastic_maps)-2)
				state_by_state_insituspec_average_BSM = round(state_by_state_insituspec_BSM,nchar(n_trees*n_stochastic_maps)-2)
				state_by_state_founderevent_average_BSM = round(state_by_state_founderevent_BSM,nchar(n_trees*n_stochastic_maps)-2)
				state_by_state_transitions_average_BSM = round(state_by_state_transitions_average_BSM/(n_trees*n_stochastic_maps),nchar(n_trees*n_stochastic_maps)-1)
			
			toc()

############################################
# average number of events per replication #
############################################

	tic()
	print("Get the average number of events for each replicate...")

	# create table for each replicate and event type
		total_events = matrix(0, ncol = 2+length(transition_types), nrow = n_trees*n_stochastic_maps)
		colnames(total_events) <- c("tree_replicate", "BSM_replicate", transition_types)
		BSM_transitions <- as.data.frame(BSM_transitions)
		counter_average <- 1
	
	# get total number of events in each replicate
		for(a in 1:n_trees){
			for(b in 1:n_stochastic_maps){
				total_events[counter_average, 1] <- a #number of tree replicate
				total_events[counter_average, 2] <- b #number of BSM replicata
				#filters the table to include just current tree/BSM replicate
					temp <- BSM_transitions %>% filter(BSM_transitions$replicate==a & BSM_transitions$BSMreplicate==b) 
					for(c in 1:length(transition_types)){ #for each event type...
						#sums instances of that event type for that tree/BSM replicate and populates the table
							total_events[counter_average, 2+c] <- sum(temp$transitionCladogType == transition_types[c])+ sum(temp$transitionAnagType == transition_types[c]) 
					}
				counter_average <- counter_average + 1
			}
		}
	
	# get mean and SD values for all BSM for each replicates, and populates ReplicateParameters table
		total_events <- as.data.frame(total_events)
		#for each tree...
			for(a in 1:n_trees){
				#filter table to include only data for that tree
					temp2 <- total_events %>% filter (total_events$tree_replicate == a)
				for(b in 1:length(transition_types)){
					col_number <- which(colnames(ReplicationParameters)==paste(transition_types[b], "_mean", sep=""))
					ReplicationParameters[a, col_number] <- round(mean(as.numeric(unlist(temp2[b+2]))),3) #calculates mean
					col_number <- which(colnames(ReplicationParameters)==paste(transition_types[b], "_SD", sep=""))
					ReplicationParameters[a, col_number] <- round(sd(as.numeric(unlist(temp2[b+2]))),3) #calculates SD
				}
		}
		
	write.csv(total_events, file = paste(dir_name,"/Number_of_events_in_each_BSM.csv",sep=""))

	toc()
	
###########################################################
# get average number of lineages per area per time period #
###########################################################

	if (plot_LTT) {
		CountStatesByTime_all_reps <- as.data.frame(CountStatesByTime_all_reps)
		# get the time periods up to the oldest tree age
			TimePoint <- seq(from= 0, to= tMax_allreps , by= time_period_size)
			names_single_areas_ <- c(names_single_areas, "LTT")
		
		# create an empty matrix to store results
			averageLTTbyarea = as.data.frame(matrix(NA, ncol = length(names_single_areas_), nrow = length(TimePoint)))
			colnames(averageLTTbyarea) <- names_single_areas_
			rownames(averageLTTbyarea) <- TimePoint
			averageLTTbyarea_95_lower <- averageLTTbyarea
			averageLTTbyarea_95_upper <- averageLTTbyarea
			counter_LTT_average <- 1
		
		for (t in seq(0, tMax_allreps, time_period_size)) {
			for (n in 1:length(names_single_areas_)){
				col_position <- which(colnames(CountStatesByTime_all_reps)==names_single_areas_[[n]])
				row_time_period <- which(CountStatesByTime_all_reps[,"TimePoint"]==t)
				mean_by_time_period <- vector(length=0)
				for (w in 1:length(row_time_period)){
					mean_by_time_period <- c(mean_by_time_period, as.numeric(CountStatesByTime_all_reps[row_time_period[[w]],col_position]))
				}
				mean_LTT <- mean(mean_by_time_period)
				averageLTTbyarea[counter_LTT_average,names_single_areas_[[n]]] <- mean_LTT
				if (length(row_time_period)>1){
					error_LTT <- sd(mean_by_time_period)/sqrt(length(row_time_period))
					low95 <- mean_LTT -(2*error_LTT)
					if(low95 < 0) {averageLTTbyarea_95_lower[counter_LTT_average,names_single_areas_[[n]]] <- 0} else {averageLTTbyarea_95_lower[counter_LTT_average,names_single_areas_[[n]]] <- low95}
					averageLTTbyarea_95_upper[counter_LTT_average,names_single_areas_[[n]]] <- mean_LTT +(2*error_LTT)
				} else {
					averageLTTbyarea_95_upper[counter_LTT_average,names_single_areas_[[n]]] <- mean_LTT
					averageLTTbyarea_95_lower[counter_LTT_average,names_single_areas_[[n]]] <- mean_LTT
				}					
			}
			counter_LTT_average <- 1 + counter_LTT_average
		}
	}
	
############################################################
# plot average number of lineages per area per time period #
############################################################
	
	if (plot_LTT) {
			averageLTTbyarea_print <- averageLTTbyarea
			averageLTTbyarea[averageLTTbyarea==0] <- NA
			averageLTTbyarea_95_upper[averageLTTbyarea_95_upper==0] <- NA
			averageLTTbyarea_95_lower[averageLTTbyarea_95_lower==0] <- NA
			for (a in 1:length(averageLTTbyarea[1,])){
				if(length(which(is.na(averageLTTbyarea[, a])))==0) next
				first_na <- min(which(is.na(averageLTTbyarea[,a])))
				if(first_na ==1) next
				averageLTTbyarea[first_na, a] <- 0
				averageLTTbyarea_95_upper[first_na, a] 
				averageLTTbyarea_95_lower[first_na, a] 
			}
		
		# get same state colors as in the BioGeoBEARS plot...
			labels_legend <- lapply(names_single_areas, sub, pattern="SUM_", replacement="")
			colors <- get_colors_for_numareas(length(names_single_areas))
			color_legends=vector(length=length(colors[1,]))
			
		# plots graph
			while (!is.null(dev.list())) {dev.off()} #cleans plotting area
			win.metafile(file=paste(dir_name, "/Average_LTT_by_area.emf", sep =""))				
			plot(x=TimePoint, y=rep(-1,length(TimePoint)),ylim= range(c(0,max(averageLTTbyarea[1,1:length(names_single_areas)], na.rm=T))), xlim=rev(range(TimePoint)), xlab="Age", ylab="Lineages", type="l", col=rgb(0,0,0,0),frame.plot=FALSE, cex.axis=1.5, cex.lab=2, yaxs="i")
				#adds vertical bars for time slices
					if(time_stratified) {
						time_slices <- head(scan(time_stratified_fn),-1)
						abline(v=time_slices, col ="gray", lwd=2.5, lty=2)
					}
				#adds 95% interval for LTT for each area -- if you have chosen to do so
					if(plot_95) {
						for (a in 1:length(names_single_areas)){
							if(is.na(averageLTTbyarea[1,a])) next
							lines(x=TimePoint, y=averageLTTbyarea_95_upper[,a],lwd=1, lty=2, col=rgb(colors[[1,a]],colors[[2,a]],colors[[3,a]], max=255))
						}
						for (a in 1:length(names_single_areas)){
							if(is.na(averageLTTbyarea[1,a])) next
							lines(x=TimePoint, y=averageLTTbyarea_95_lower[,a],lwd=1, lty=2, col=rgb(colors[[1,a]],colors[[2,a]],colors[[3,a]], max=255))
						}
					}
				#adds lines with average LTT for each area
					for (a in 1:length(names_single_areas)){
						color_legends[a]<-rgb(colors[[1,a]],colors[[2,a]],colors[[3,a]], max=255)
						if(is.na(averageLTTbyarea[1,a])) next
						lines(x=TimePoint, y=averageLTTbyarea[,a],lwd=2, col=rgb(colors[[1,a]],colors[[2,a]],colors[[3,a]], max=255))
					}
				#adds legend
					legend(x=max(TimePoint), y=max(averageLTTbyarea[1,1:length(names_single_areas)], na.rm=T), legend=labels_legend, fill=color_legends, cex=1)
			while (!is.null(dev.list())) {dev.off()}
		}
		
########################################
# get and sort most common transitions #
########################################

	#create tables to store most common events
		common_transitions <- matrix(NA, ncol = 5, nrow= length(statelabels)^2)
		colnames(common_transitions) <- c("starting_range", "ending_range", "average_n_events", "SD", "transition_type")
	
	#parse tables with all transitions
		i=1
		for (a in 1:length(statelabels)) {
			for (b in 1:length(statelabels)) {
					range1 <- statelabels[a]
					range2 <- statelabels[b]
					common_transitions[i,1] <- range1
					common_transitions[i,2] <- range2
					common_transitions[i,3] <- state_by_state_transitions_average_BSM[a, b]
					common_transitions[i,5] <- paste(state_by_state_Cladogenetic_events[range1, range2], " or ", state_by_state_Anagenetic_events[range1, range2], sep="")
				i <- i+1
			}
		}
	
	#sort table by count of events
	#drop events that have not been recorded
		common_transitions <- as.data.frame(common_transitions)
		common_transitions <- common_transitions[order(as.numeric(common_transitions$average_n_events), decreasing=TRUE),]
		common_transitions <-common_transitions[common_transitions$average_n_events>0.05,]
		strings_to_replace <- c ("X or ", "- or ", " or –", " or X", "_") #drops meaningless / impossible transitions
		replace_by <- c ("", "", "", "", " ")
		for (a in 1:length(strings_to_replace)){
			common_transitions$transition_type <- gsub(strings_to_replace[a], replace_by[a], common_transitions$transition_type)
		}

#################################
# get SD for most common events #
#################################

	if (get_SD_common_events) {
		tic()
		print("Calculate SD for most common events...")
		print("This is slow... sorry! :(")
	
		# create table for each replicate and transition type, but including only those among most common transitions
			number_transitions_SD <- length(common_transitions[,1])
			SD_common_events = matrix(0, ncol = 2+(number_transitions_SD), nrow = n_trees*n_stochastic_maps)
			SD_common_events <- as.data.frame(SD_common_events)
			SD_col_names <- c("tree_replicate", "BSM_replicate")
			for (p in 1:number_transitions_SD){ #small detour to get all the names...
				SD_col_names <- c(SD_col_names, paste(common_transitions[p,1],">",common_transitions[p,2],sep=""))
			}
			colnames(SD_common_events) <- SD_col_names
			counter_SD <- 1
		
		# get total number of each transition in each replicate
			for(a in 1:n_trees){
				for(b in 1:n_stochastic_maps){
					SD_common_events [counter_SD, 1] <- a
					SD_common_events [counter_SD, 2] <- b
					#filter table to include only currently tree/replication
					temp4 <- BSM_transitions %>% filter(BSM_transitions$replicate==a & BSM_transitions$BSMreplicate==b)	
					for(c in 1:number_transitions_SD){
						#define the transition to look for
						look_transition <- paste(common_transitions[c,1], "->", common_transitions[c,2])
						#gets the sum of this particular transition for this tree/BSM replicate
						SD_common_events[counter_SD, 2+c] <- sum(temp4$transitionCladog==look_transition) + sum(temp4$transitionAnag==look_transition)
					}
					counter_SD <- counter_SD + 1
				}
				print(paste(a, "/", n_trees, " completed...", sep=""))
			}
			
			#calculates SD for each transition
				for(i in 1:number_transitions_SD){
					common_transitions[i,4] <- round(sd(as.numeric(unlist(SD_common_events[i+2]))),3)
				}
			toc()
	}
	
#################
# write outputs #
#################
	
	#AIC, d, e, LnL, etc. for each of the trees
		write.csv(ReplicationParameters, file = paste(dir_name,"/Replication_parameters.csv",sep=""))
	# giant table including all biogeographic events, their ages, and the phylogenetic node / branch where they occcur
		write.csv(BSM_transitions, file = paste(dir_name,"/All_transitions_by_node_BSM.csv",sep=""))
	# table with all states vs. all states, and the observed average number of transitions (transitions are row state to column state)
		write.csv(state_by_state_transitions_average_BSM, file = paste(dir_name,"/Average_state_by_state_all_transitions_BSM.csv",sep=""))
	# same as above, but broken down by type of transition
		write.csv(state_by_state_dispersal_average_BSM, file = paste(dir_name,"/state_by_state_by_event_tables/Average_state_by_state_dispersal_BSM.csv",sep=""))
		write.csv(state_by_state_extinction_average_BSM, file = paste(dir_name,"/state_by_state_by_event_tables/Average_state_by_state_extinction_BSM.csv",sep=""))
		write.csv(state_by_state_vicar_average_BSM, file = paste(dir_name,"/state_by_state_by_event_tables/Average_state_by_state_vicar_BSM.csv",sep=""))
		write.csv(state_by_state_subset_sympatry_average_BSM, file = paste(dir_name,"/state_by_state_by_event_tables/Average_state_by_state_subset_sympatry_BSM.csv",sep=""))
		write.csv(state_by_state_insituspec_average_BSM, file = paste(dir_name,"/state_by_state_by_event_tables/Average_in_situ_speciation_BSM.csv",sep=""))
		write.csv(state_by_state_founderevent_average_BSM, file = paste(dir_name,"/state_by_state_by_event_tables/Average_state_by_state_jump_dispersal_BSM.csv",sep=""))
	# number of lineages in each area in each time_period, for each BSM replicate
		if (plot_LTT) write.csv(CountStatesByTime_all_reps, file = paste(dir_name,"/CountStatesByTime_all_reps.csv",sep=""))
	# average number of lineages in each area in each time_period, average over all BSM replicates
		if (plot_LTT) write.csv(averageLTTbyarea_print, file = paste(dir_name,"/Average_LTT_by_area.csv",sep=""))
	# most likely states for the root node
		write.csv(RootState, file = paste(dir_name,"/RootState.csv",sep=""))
	# a user-friendly summary of the most common transitions, summarized directly from state_by_state_transitions_average_BSM
		write.csv(common_transitions, file = paste(dir_name,"/Common_transitions.csv",sep=""))
	# total number of events in each BSM
		write.csv(total_events, file = paste(dir_name,"/tables/Number_of_events_in_each_BSM.csv",sep=""))
	}}

#hope it didn't take TOO long...
	print("Total time elapsed:")
	Sys.time() - start_time 

sink(file = NULL)
#acabou!