#! /usr/bin/env python

#################################################################################
#	amp.py v0.98					 											#
#									 											#
#	Copyright Jeremy M. Brown, 2010-2017										#
#	jeremymbrown@gmail.com														#
#																				#
#  This program is free software; you can redistribute it and/or modify			#
#  it under the terms of the GNU General Public License as published by			#			
#  the Free Software Foundation; either version 3 of the License, or			#
#  (at your option) any later version.											#
#																				#
#  This program is distributed in the hope that it will be useful,				#
#  but WITHOUT ANY WARRANTY; without even the implied warranty of				#
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the				#
#  GNU General Public License for more details.									#
#																				#
#  You should have received a copy of the GNU General Public License along		#
#  with this program. If not, see <http://www.gnu.org/licenses/>.				#
#																				#
#################################################################################

# Version changes since 0.91:
#	- Implements branch-specific likelihood ratio test statistic
#
# Version changes since 0.92:
#	- Corrects the implementation of the branch-specific likelihood ratio test statistic to
#		calculate the likelihood ratio between a positively constrained search and a negatively
#		constrained search for each branch, rather than between an unconstrained search and a
#		negatively constrained search
# 
# Version changes since 0.93:
#	- Corrects bug in 0.91-0.92 that failed to order the vector of RF distances when calculating
#       quantile-based test statistics.
#
# Version changes since 0.94:
#   - Corrects bug in 0.91-0.93 that failed to order the vector of RF distances when calculating
#		interquartile-based test statistics
#
# Version changes since 0.95:
#   - Implements Bayes factor test statistic for reading output of model-switch thermodynamic
#     integration program written by me
#
# Version changes since 0.96:
#	- Implements Bayes factor test statistic for reading output of steppingstone sampling from 
#     MrBayes 3.2.1. Expects two kinds of runs in a folder: some that positively constrain each 
#     branch in a set of branches and some that negatively constrain each branch in that set.
#
# Version changes since 0.97:
#	- Now pre-calculates the ordered vector of RF distances between trees whenever quantile
#	  or IQR tests will be performed.  Can also provide multiple quantile values on a single
#	  command line to avoid having to rerun AMP for different positions.
#	- Now implementing a variance in tree length statistic.

"""
Program for the calculation of phylogenetic model adequacy test statistics based on analyses of
posterior predictive datasets.


NOTE: Requires the previous installation of Dendropy v3
"""

import math
import traceback
import getopt
import sys
import dendropy
from dendropy import treecalc
import signal
import time
import threading


## Use this bool to turn verbose debugging output on and off
debug = False

## Use this bool to turn threading on and off
threadBool = False

## Use this bool to turn on and off the display of elapsed time for reading v. calculating
displayTiming = True

def main(argv):
	
	print """
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                                 #
#   AMP: Assessing Model adequacy with Predictive distributions   #
#                                                                 #
#                             v0.98                               #
#                                                                 #
#                        Jeremy M. Brown                          #
#                     jeremymbrown@gmail.com                      #
#                                                                 #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
	"""
	
	# Sends Ctrl-C interrupt signal to sigint_handler method
	signal.signal(signal.SIGINT,sigint_handler)
	
	# Parsing command-line arguments and options
	try:
		opts, args = getopt.getopt(argv,"q:eip:TVm:cluto:vL:b:B:",[])
	except:	# If command-line input not proper, print usage statement and exit
		print "Problem with command-line input!"
		usage()
		if debug:
			traceback.print_exc()
		sys.exit(1)
		
	# Evaluates command-line input and stores selections in "values" tuple
	values = parse_in(opts,args)
	
	# Data structures to store test output for pp datasets
	pp_quantiles = [ ]
	pp_entropies = [ ]
	pp_iqrs = [ ]
	pp_bpps = [ ]
	pp_tls = [ ]
	pp_tlVars = [ ]
	triBipart = [ ]
	likeRatios = [ ]  # This will be a list of lists, where each component list corresponds to all
					  #  likelihood ratio values for a particular dataset.
	allBFs = [  ] # As with likelihood ratios, this will be a list of lists, where each component
				  # list corresponds to all Bayes factors for a particular dataset.
	allmbSSbfs = [  ]  # As above for model-switch TI BFs
	
	"""
	# QUARTET  --  NEVER FULLY IMPLEMENTED
	# Calculates taxon_sets for three bipartitions (based on quartets in empirical greedy consensus 
	#	tree) for use in calculating the partition-specific entropy statistic
	if values[5]:
	
		# Reads in empirical tree list
		emp_tree_list = read_emp_trees(values[17],  # basename
								   int(values[19]), # number of replicate analyses
								   values[8],	# boolean for manual burn-in 
								   int(values[9]))  # manual burn-in (default = 0)
		
		# Calculates empirical greedy consensus tree
		empGreedyCon = emp_tree_list.consensus(min_freq=0)
		
		# Calculates appropriate taxon sets
		triBipart = getBiparts(empGreedyCon,	# Empirical greedy consensus tree
								values[6]))		# List of taxon names
	"""
	
	# Variables to hold amount of time spent reading trees versus calculating test statistics
	read_time = 0
	calc_time = 0
	
	# Read in, calculate, and store test values for each pp dataset

	print """
	Calculating test statistics for pp datasets...
	
		Replicates completed: """,
	sys.stdout.flush()
	
	for i in range(int(values[18])+1)[1:]:		# Iterates over 1,...,data_no
		
		############# Reading in data (either posteriors or likelihoods) #############
		
		read_start = time.clock()
		
		pp_tree_list = [ ]
		likeScores = [ ]
		bfs = [ ]
		mbSSbfs = [ ]
		
		# Checks to make sure a marginal test statistic has been selected
		if (values[0] | values[3] | values[4] | values[5] | values[7] | values[26]):
			# Read in posterior predictive tree distributions from file (excluding burn-in)
			pp_tree_list = read_trees(values[17], # basename
					   				   i, # relevant posterior predictive dataset number
			   						   int(values[19]), # number of replicate analyses per dataset
			   						   values[8], # boolean for manual burn-in (always opposite of boolean for MrC burnin)
			   						   int(values[9])) # manual burn-in value (default = 0)
		
		# Checks to see if bipartition-specific LR statistic has been selected
		if (values[20]):
			# Read in likelihood scores for unconstrained and any constrained analyses
			likeScores = read_likelihoods(values[17],	# Basename
										  values[21],	# Number of likelihood constraints
										  i)		# Relevant posterior predictive dataset number
			
		# Checks to see if bipartition-specific model-switch TI BF statistic has been selected
		if values[22]:
			# Read in Bayes factors for each replicate
			bfs = readBFs(values[17],	# Basename
						  values[23],	# Number of bf constraints
						  i)			# Relevant posterior predictive dataset number
		
		# Checks to see if bipartition-specific MrBayes steppingstone BF statistic has been selected
		if values[24]:
			# Read in marginal likelihoods from pos and neg log files and calculates the corresponding BF
			mbSSbfs = readmbSSbfs(values[17], # Basename
								  values[25], # Number of MrBayes steppingstone constraints
								  i)		  # Relevant posterior predictive dataset number
			
		read_end = time.clock()
		
		read_time += (read_end-read_start)
		
		########################### Finished Reading Data ###########################
		
		"""
		Next section calculates all requested test statistics for pp datasets.
	
		To add new test statistics, just:
			(1) update input parser and add new entry to be returned to values (above)
			(2) add boolean test below and corresponding method (both for pp and emp datasets)
			(3) add additional call to p_value()
		"""
		
		calc_start = time.clock()
		
		# Calculate ordered RF vector if either quantile-position or IQR stats selected
		if (values[0] | values[4]):
			pp_rf_dists = calcAndSortTreeDists(dendropy.DataSet(pp_tree_list))
		
		# Calculate quantile-based test statistic
		if values[0]:
			if threadBool:	# Threaded version
				ppQuantCalc = ppQuantileThread(dendropy.DataSet(pp_tree_list), # Copies and passes DataSet object with one TreeList
												  values[1],	# int -- k-th quantile
												  values[2],	# int -- p total quantiles	
												  pp_quantiles)
				ppQuantCalc.start()
			else:			# Simple serial version
				pp_quantiles.append(quantile_test(pp_rf_dists, # Passes vector of ordered RF distances
												  values[1],	# list of ints -- list of k-th quantiles
												  values[2]))	# list of ints -- list of p total quantiles			
	
		# Calculate entropy-based test statistic
		if values[3]:
			pp_entropies.append(entropy_test(dendropy.DataSet(pp_tree_list)))
				# pp_tree_list is a DataSet object with a single TreeList 
				#   for the relevant pp data set -- creates copy of this and passes to entropy_test()
	
		# Calculate interquartile range test statistic
		if values[4]:
			pp_iqrs.append(iqr_test(pp_rf_dists))	# Passes vector of ordered RF distances
	
		# Calculate partition-specific entropy test statistic
		if values[5]:
			pp_bpps.append(partition_test(dendropy.DataSet(pp_tree_list),
										values[6])) # Taxa used to define partition
		
		# QUARTET
		#if values[5]:
		#	pp_bpps.append(partition_test(dendropy.DataSet(pp_tree_list),triBipart))	
	
		# Calculate mean treelength test statistic
		if values[7]:
			pp_tls.append(treelength_test(dendropy.DataSet(pp_tree_list)))	

		# Calculate treelength variance test statistic
		if values[26]:
			pp_tlVars.append(tlVar_test(dendropy.DataSet(pp_tree_list)))

		# Calculate partition-specific likelihood ratio test statistic
		if values[20]:
			tempLikeRatios = [ ]
			for j in range(len(likeScores[1:])+1): # Iterates through branches and stores LRs in tempLikeRatios
				if j % 2 == 0:	# Skips odd numbers (corresponding to neg constraints)
					ratio = likeScores[j]-likeScores[j+1] # Calcs ln(LR) as ln(posConL)-ln(negConL)   
					if ratio < 0.001 and ratio > -0.001: # Rounds down to zero when two scores are very similar
						ratio = 0
					tempLikeRatios.append( ratio )
			likeRatios.append( tempLikeRatios )
	
		# Calculate partition-specific Bayes factor test statistic
		# 	This method is trivial, but I haven't bundled it with the reading of the values from file
		# 	to maintain consistency with other test statistic calculations
		if values[22]:
			allBFs.append( bfs )
			
		# Calculate partition-specific BF test statistic - uses marg. likelihoods from MrBayes steppingstone
		#   sampling of constrained runs. As above, this method is trivial but maintained here for consistency.		
		if values[24]:
			allmbSSbfs.append( mbSSbfs )
			
		## Waits here for any spawned threads to finish
		if values[0] and threadBool:	# Checking for quantile test stat and threading bools
			while ppQuantCalc.isAlive():
				time.sleep(1)
		
		calc_end = time.clock()
		
		calc_time += (calc_end-calc_start)

		print "%d" % (i),
		sys.stdout.flush()

	print ""

	if displayTiming:
		print """
		Time spent reading trees: %f s""" % read_time
		print "		Time spent calculating test statistics: %f s" % calc_time

	if debug:
		print """pp_entropies:
		"""
		print pp_entropies

	print """
	Calculating test statistics for empirical dataset...
	"""
	
	emp_tree_list = [ ]
	empLikeRatios = [ ]
	
	# Checks to make sure a marginal test statistic has been selected
	if (values[0] | values[3] | values[4] | values[5] | values[7] | values[26]):
		emp_tree_list = read_emp_trees(values[17],  # basename
									   int(values[19]), # number of replicate analyses
									   values[8],	# boolean for manual burn-in 
									   int(values[9]))  # manual burn-in (default = 0)

	# Checks to see if bipartition-specific LR statistic has been selected
	if (values[20]):
	
		# Reads empirical likelihood
		empLikeScores = read_emp_likelihoods(values[17],	# Basename
										  	 values[21])	# Number of constraints
		
		# Calculates LRs for empirical data
		empLikeRatios = [ ]
		for i in range(len(empLikeScores[1:])+1):
			if i % 2 == 0:
				empRatio = empLikeScores[i]-empLikeScores[i+1]
				if empRatio < 0.001 and empRatio > -0.001:
					empRatio = 0
				empLikeRatios.append( empRatio )
	
	# Checks to see if bipartition-specific model-switch TI Bayes factor statistic has been selected
	if (values[22]):
		# Reads empirical Bayes factors
		empBFs = read_emp_bfs(values[17],	# Basename
							 values[23])	# Number of bf constraints
	
	# Checks to see if MrBayes bipartition-specific Bayes factor statistic has been selected
	if (values[24]):
		# Reads empirical marginal likelihoods and calculates corresponding Bayes factors
		empmbSSbfs = read_emp_mbSSbfs(values[17],	# Basename
									  values[25])	# Number of mbSSbfConstraints

	# Calculate ordered RF vector if either quantile-position or IQR stats selected
	if (values[0] | values[4]):
		emp_rf_dists = calcAndSortTreeDists(dendropy.DataSet(emp_tree_list))

	# Calculates quantile-based test statistic
	if values[0]:
		emp_quantile = quantile_test(emp_rf_dists,   # Passes ordered vector of RF distances
									 values[1],		 # int -- k-th quantile
									 values[2])		 # int -- p total quantiles

	# Calculates empirical entropy test statistic
	if values[3]:	
		emp_entropy = entropy_test(dendropy.DataSet(emp_tree_list))

	# Calculates interquartile range test statistic
	if values[4]:
		emp_iqr = iqr_test(emp_rf_dists)
		
	# Calculates partition-specific entropy test statistic
	if values[5]:
		emp_bpp = partition_test(dendropy.DataSet(emp_tree_list),
								values[6])		# Taxa used to define partition
		
	# Calculates mean treelength test statistic
	if values[7]:
		emp_tl = treelength_test(dendropy.DataSet(emp_tree_list))

	# Calculates treelength variance test statistic
	if values[26]:
		emp_tlVar = tlVar_test(dendropy.DataSet(emp_tree_list))

	print """
	Calculating P-values...
	"""
	if values[14]:
		logfile = open(values[15],'w')
	
	if values[0]:
		for i in range(len(values[1])):
			p_value(emp_quantile[i],
					[row[i] for row in pp_quantiles], # Extracts column i from pp_quantiles
					("%s-th %s-quantile Test Statistic" % (str(values[1][i]),str(values[2][i]))),
					values[11],
					values[12],
					values[13],
					values[14],
					values[15],
					values[16],
					logfile)
	
	if values[3]:  # Calculates (and outputs) appropriate p-values for statistical entropy test statistic
		p_value(emp_entropy,
				pp_entropies,
				"Entropy-based Test Statistic",
				values[11],		# Bool for lower one-tailed p-value
				values[12],		# Bool for upper one-tailed p-value
				values[13],		# Bool for two-tailed p-value
				values[14],		# Bool for output to file
				values[15],		# Output file name (string)
				values[16],		# Bool for verbosity
				logfile)

 	if values[4]:  # Calculates (and outputs) specified p-values for interquartile range test statistic
		p_value(emp_iqr,
				pp_iqrs,
				"Interquartile Range Test Statistic",
				values[11],
				values[12],
				values[13],
				values[14],
				values[15],
				values[16],
				logfile)
	
	if values[5]:
		p_value(emp_bpp,
				pp_bpps,
				"Bipartition Frequency Test Statistic",
				values[11],
				values[12],
				values[13],
				values[14],
				values[15],
				values[16],
				logfile)
	
	if values[7]:
		p_value(emp_tl,
				pp_tls,
				"Mean Treelength Test Statistic",
				values[11],
				values[12],
				values[13],
				values[14],
				values[15],
				values[16],
				logfile)
	
	if values[26]:
		p_value(emp_tlVar,
				pp_tlVars,
				"Treelength Variance Test Statistic",
				values[11],
				values[12],
				values[13],
				values[14],
				values[15],
				values[16],
				logfile)				
	
	if values[20]:
		for i in range(len(empLikeRatios)):
			p_value(empLikeRatios[i],
					 [j[i] for j in likeRatios],
					 "Bipartition-Specific Likelihood Ratio Constraint %d Test Statistic" % (i+1,),
					 values[11],		# Bool for lower one-tailed p-value
					 values[12],		# Bool for upper one-tailed p-value
					 values[13],		# Bool for two-tailed p-value
					 values[14],		# Bool for output to file
					 values[15],		# Output file name (string)
					 values[16],		# Bool for verbosity
					 logfile)
	
	if values[22]:
		for i in range(len(empBFs)):
			p_value(empBFs[i],
					[j[i] for j in allBFs],
					"Bipartition-Specific Bayes Factor Constraint %d Test Statistic" % (i+1,),
					values[11],		# Bool for lower one-tailed p-value
					values[12],		# Bool for upper one-tailed p-value
					values[13],		# Bool for two-tailed p-value
					values[14],		# Bool for output to file
					values[15],		# Output file name (string)
					values[16],		# Bool for verbosity
					logfile)
	
	if values[24]:
		for i in range(len(empmbSSbfs)):
			p_value(empmbSSbfs[i],
					[j[i] for j in allmbSSbfs],
					"Bipartition-Specific MrBayes Steppingstone Bayes Factor Constraint %d Test Statistic" % (i+1,),
					values[11],		# Bool for lower one-tailed p-value
					values[12],		# Bool for upper one-tailed p-value
					values[13],		# Bool for two-tailed p-value
					values[14],		# Bool for output to file
					values[15],		# Output file name (string)
					values[16],		# Bool for verbosity
					logfile)
	
	if values[14]:
		logfile.close()

	print """
	Program execution complete.  Exiting...
	"""


def p_value(emp_stat,pp_stats,test_stat_name,lower,upper,two,tofile,outfile,verbose,logfile):
	"""
	Calculates specified p-values
	"""
	
	less_than_count = 0
	for i in pp_stats:
		if i <= emp_stat:
			less_than_count += 1
	lower_p = float(less_than_count)/float(len(pp_stats))
	if lower and debug:
		print "Lower One-tailed P-value: %f" % lower_p
			
	greater_than_count = 0
	for i in pp_stats:
		if i >= emp_stat:
			greater_than_count += 1
	upper_p = float(greater_than_count)/float(len(pp_stats))
	if upper and debug:
		print "Upper One-tailed P-value: %f" % upper_p
		
	two_p = min(2*min(lower_p,upper_p),1)
	if two and debug:
		print "Two-tailed P-value: %f" % two_p
		
	if tofile:
	
		logfile.write("****** %s ******\n" % test_stat_name)
		logfile.write('\n')
	
		if verbose:
			## Outputs all posterior predictive test statistic values
			logfile.write("Posterior Predictive %ss:\n" % test_stat_name)
			logfile.write('\n')
			for i in pp_stats:
				logfile.write(str(i)+'\n')
				
			## Outputs empirical test statistic value
			logfile.write('\n')
			logfile.write("Empirical %s:\n" % test_stat_name)
			logfile.write('\n')
			logfile.write(str(emp_stat)+'\n')
			
		## Outputs appropriate p-values
		logfile.write('\n')
		logfile.write("P-values:\n")
		logfile.write('\n')
		if lower:
			logfile.write("Lower One-tailed P-value: %f\n" % lower_p)
		if upper:
			logfile.write("Upper One-tailed P-value: %f\n" % upper_p)
		if two:
			logfile.write("Two-tailed P-value: %f\n" % two_p)
		logfile.write('\n')

def read_emp_trees(basename,rep_no,manual,burnin):
	"""
	Reads in trees resulting from analysis of the original empirical data set
	"""

	# Instantiates taxon set object for tree list
	taxa = dendropy.TaxonSet()
	
	# Instantiates data set object to hold tree list
	emp_trees = dendropy.DataSet()
	emp_trees.add_taxon_set(taxa)
	
	try:
		if not manual:
			burnin = empGetMrCburn(basename)
		temp_tree_list = dendropy.TreeList()
		for j in range(rep_no+1)[1:]:
			filename = "%s_emp_r%d.t" % (basename,j)
			temp_tree_list.read_from_path(filename,"nexus",tree_offset=burnin,as_unrooted=True,taxon_set=taxa)
		emp_trees.add(temp_tree_list)
	except:
		print "Problem reading in empirical trees. Exiting..."
		if debug:
			traceback.print_exc()
		sys.exit(1)

	# Returns a data set object (for consistency with read_trees()) with a single tree list
	return(emp_trees)

def empGetMrCburn(basename):
	"""
	empGetMrCburn() retrieves appropriate burnin for empirical MrConverge log file
	"""
	
	try:
		mrcin = open("%s_emp.log" % (basename,))	# Opens log file for reading
		temp = mrcin.readline()							# Reads first line
		while temp.find("!CONVERGENCE REACHED!") == -1:	# Iterates through lines until convergence
			if temp.find("BURNIN") != -1:				# Stores burnin values as they are found
				burnin = int(temp.split()[3])
			temp = mrcin.readline()
		mrcin.close()
	except:
		print """
		Problem getting burn-in value from empirical MrConverge log file. Exiting...
		"""
		if debug:
			traceback.print_exc()
		sys.exit(1)
		
	if debug:	
		print "empirical MrC burn: %d" % burnin
		
	return(burnin)
	

def read_trees(basename,data_no,rep_no,manual,burnin):
	"""
	Reads in trees resulting from analyses of posterior predictive datasets.
	"""

	# Instantiates taxon set object to be shared by tree lists
	taxa = dendropy.TaxonSet()
	
	# Instantiates data set object to hold all lists of trees
	pp_trees = dendropy.DataSet()
	pp_trees.add_taxon_set(taxa)

	try:
		if not manual:
			burnin = getMrCburn(basename,data_no)
		temp_tree_list = dendropy.TreeList()  # Instantiates temp tree list to hold current data set trees
		for j in range(rep_no+1)[1:]:	# Iterates over 1,...,rep_no
			filename = "%s_%d_r%d.t" % (basename,data_no,j)
			temp_tree_list.read_from_path(filename,"nexus",tree_offset=burnin,as_unrooted=True,taxon_set=taxa)
		pp_trees.add(temp_tree_list)	# Adds new tree list to the "pp_trees" data set object
	except:
		print """
	Problem reading in trees.  Exiting...
		"""
		if debug:
			traceback.print_exc()
		sys.exit(1)
		
	## Currently some error with the DataSet.unify_taxa() method [1.7.10]
	# pp_trees.unify_taxa()	# Just to make sure all trees share the same TaxonSet
	
	# Returns a dataset object containing separate treelists for each pp dataset
	return(pp_trees)
	
def	read_emp_likelihoods(basename,numCon):
	"""
	Reads in likelihoods for positively and negatively constrained searches for empirical dataset.  
	Returns a list containing all likelihood scores for the empirical dataset.
	
	Example for N constraints:
	[posCon1,negCon1,posCon2,negCon2,...,posConN,negConN]		
	"""
	
	likelihoods = []
	
	try:
		for i in range(int(numCon)+1)[1:]:
			
			# Reading in positive constraint likelihoods
			likeIn = open("%s_emp_bp%d.pos.constraint.best.tre" % (basename,i))
			tempLine = likeIn.readline()
			while (tempLine.find("!GarliScore") == -1):
				tempLine = likeIn.readline()
			scoreList = tempLine.split("][")
			score = float(scoreList[1].split(" ")[1])
			likelihoods.append(score)
			
			# Reading in negative constraint likelihoods
			likeIn = open("%s_emp_bp%d.neg.constraint.best.tre" % (basename,i))
			tempLine = likeIn.readline()
			while (tempLine.find("!GarliScore") == -1):
				tempLine = likeIn.readline()
			scoreList = tempLine.split("][")
			score = float(scoreList[1].split(" ")[1])
			likelihoods.append(score)				
	except:
		print "Problem reading in likelihood scores.  Exiting...."
		if debug:
			traceback.print_exc()
		sys.exit(1)
	
	return likelihoods
	
def	read_likelihoods(basename,numCon,data_no):
	"""
	Reads in likelihoods for positively and negatively constrained searches for each dataset.  
	Returns a list containing all likelihood scores for each replicate dataset.
	
	Example for N constraints:
	[posCon1,negCon1,posCon2,negCon2,...,posConN,negConN]		
	"""
	
	likelihoods = []
	
	try:
		for i in range(int(numCon)+1)[1:]:	# Iterates from 1 to numCon
			
			# Storing positive constraint likelihoods
			likeIn = open("%s_%d_bp%d.pos.constraint.best.tre" % (basename,data_no,i))
			tempLine = likeIn.readline()
			while (tempLine.find("!GarliScore") == -1):
				tempLine = likeIn.readline()
			scoreList = tempLine.split("][")
			score = float(scoreList[1].split(" ")[1])
			likelihoods.append(score)	
			
			# Storing negative constraint likelihoods
			likeIn = open("%s_%d_bp%d.neg.constraint.best.tre" % (basename,data_no,i))
			tempLine = likeIn.readline()
			while (tempLine.find("!GarliScore") == -1):
				tempLine = likeIn.readline()
			scoreList = tempLine.split("][")
			score = float(scoreList[1].split(" ")[1])
			likelihoods.append(score)				
	except:
		print "Problem reading in likelihood scores.  Exiting...."
		if debug:
			traceback.print_exc()
		sys.exit(1)
	
	return likelihoods
	
def read_emp_bfs(basename,numCon):
	"""
	Reads in ln(BF) estimated by model-switching thermodynamic integration for the empirical dataset.
	Returns a list containing all Bayes factors for the empirical dataset.
	
	Example for N constraints:
	[ln(BF)1,ln(BF)2,...,ln(BF)N]
	"""
	bfs = [ ]
	
	try:
		for i in range(int(numCon)+1)[1:]:	# Iterates from 1 to numCon
			
			bfIn = open("%s_emp_bp%d.bf.out" % (basename,i))
			tempLine = bfIn.readline()
			while (tempLine.find("Estimated ln(Bayes Factor)") == -1):
				tempLine = bfIn.readline()
			bfVal = float(tempLine.split(" ")[4])
			bfs.append(bfVal)
			
	except:
		print "Problem reading in Bayes factors. Exiting..."
		if debug:
			traceback.print_exc()
		sys.exit(1)
	
	return bfs

	
def readBFs(basename,numCon,data_no):
	"""
	Reads in Bayes factors (actually ln(BFs)) estimated from model-switch thermodynamic 
	integration between a model that positively constrains some branch and one that negatively
	constrains it.  List simply contains ln(BF) values for each constraint.
	
	Example for N constraints:
	[ln(BF)1,ln(BF)2,...,ln(BF)N]
	"""
	bfs = [ ]
	
	try:
		for i in range(int(numCon)+1)[1:]:	# Iterates from 1 to numCon
			
			bfIn = open("%s_%d_bp%d.bf.out" % (basename,data_no,i))
			tempLine = bfIn.readline()
			while (tempLine.find("Estimated ln(Bayes Factor)") == -1):
				tempLine = bfIn.readline()
			bfVal = float(tempLine.split(" ")[4])
			bfs.append(bfVal)
			
	except:
		print "Problem reading in Bayes factors. Exiting..."
		if debug:
			traceback.print_exc()
		sys.exit(1)
	
	return bfs
	
def read_emp_mbSSbfs(basename,numCon):
	
	mbSSbfs = [ ]
	
	try:
		for i in range(int(numCon)+1)[1:]:	# Iterates from 1 to numCon
			posIn = open("%s_emp.con%d.pos.log" % (basename,i))
			negIn = open("%s_emp.con%d.neg.log" % (basename,i))
			posTempLine = posIn.readline()
			negTempLine = negIn.readline()
			while (posTempLine.find("Mean:") == -1):
				posTempLine = posIn.readline()
			while (negTempLine.find("Mean:") == -1):
				negTempLine = negIn.readline()
			posLike = float(posTempLine.strip().split()[1])
			negLike = float(negTempLine.strip().split()[1])
			bfVal = float(posLike - negLike)
			mbSSbfs.append(bfVal)
			posIn.close()
			negIn.close()
	except:
		print "Problem reading in MrBayes marginal likelihood values for empirical dataset. Exiting..."
		if debug:
			traceback.print_exc()
		sys.exit(1)
	
	return mbSSbfs
	
def readmbSSbfs(basename,numCon,data_no):

	mbSSbfs = [ ]

	try:
		for i in range(int(numCon)+1)[1:]:	# Iterates from 1 to numCon
			posIn = open("%s_%d.con%d.pos.log" % (basename,data_no,i))
			negIn = open("%s_%d.con%d.neg.log" % (basename,data_no,i))
			posTempLine = posIn.readline()
			negTempLine = negIn.readline()
			while (posTempLine.find("Mean:") == -1):
				posTempLine = posIn.readline()
			while (negTempLine.find("Mean:") == -1):
				negTempLine = negIn.readline()
			posLike = float(posTempLine.strip().split()[1])
			negLike = float(negTempLine.strip().split()[1])
			bfVal = float(posLike - negLike)
			mbSSbfs.append(bfVal)
			posIn.close()
			negIn.close()
	except:
		print "Problem reading in MrBayes marginal likelihood values for empirical dataset. Exiting..."
		if debug:
			traceback.print_exc()
		sys.exit(1)

	return mbSSbfs
	
def getMrCburn(basename,data_no):
	"""
	getMrCburn() retrieves appropriate burnins from MrConverge log files.
	
	Files should be named as basename_data#.log.
	"""
	try:
		mrcin = open("%s_%d.log" % (basename,data_no))	# Opens log file for reading
		temp = mrcin.readline()							# Reads first line
		while temp.find("!CONVERGENCE REACHED!") == -1:	# Iterates through lines until convergence
			if temp.find("BURNIN") != -1:				# Stores burnin values as they are found
				burnin = int(temp.split()[3])
			temp = mrcin.readline()
		mrcin.close()
	except:
		print "Problem getting burn-in value from MrConverge log file. Exiting..."
		if debug:
			traceback.print_exc()
		sys.exit(1)
		
	if debug:	
		print "MrC burn for file %s_%d.log: %d" % (basename,data_no,burnin)
		
	return(burnin)

"""
	# Quartet
def getBiparts(greedyConTree,taxa):
	"""
	#Function to find three bipartitions based on quartets surrounding a given bipartition in the 
	#empirical greedy consensus tree
"""
	try:
		# Defines the target bitmask based on the user input lists of taxa
		searchMask = bit_mask(taxa,greedyConTree.leaf_nodes())
		
		# Defines bitmasks for each internal node in the tree
		for i in greedyConTree.internal_nodes():
			i.label = bit_mask(i.leaf_nodes(),greedyConTree.leaf_nodes())
		
		# Locates the node corresponding to the specified bipartition using bitmasks
		focalNode = greedyConTree.find_node_with_label(searchMask)
		
		# Reroots the tree at the focal node so that the leaf_nodes() function can be used to 
		#	define the taxon sets in each part of the quartet
		greedyConTree.reroot_at(focalNode)
		
		# Defines taxon sets based on the quartet surrounding the bipartition of interest
		taxSet1 = greedyConTree.seed_node.child_nodes()[0].leaf_nodes()
		taxSet2 = greedyConTree.seed_node.child_nodes()[1].leaf_nodes()
		taxSet3 = greedyConTree.seed_node.child_nodes()[2].child_nodes()[0].leaf_nodes()
		taxSet4 = greedyConTree.seed_node.child_nodes()[2].child_nodes()[1].leaf_nodes()
		
		# NEED TO DEFINE THREE DIFFERENT PARTITIONS BASED ON THESE TAXON SETS AND RETURN THEM
		
	except:
		print "Problem getting three bipartitions from empirical greedy consensus tree. Exiting..."
		if debug:
			traceback.print_exc()
		sys.exit(1)	
"""


###############  Begin Threaded Classes ###############

# Annoyingly, python multi-threading doesn't take advantage of multiple CPU's.  Will have to recode
#   using something like parallel python (http://www.parallelpython.com/) to use multi-CPU's

class ppQuantileThread( threading.Thread ):

	def __init__ ( self,tree_data,k,p,pp_quantiles ):
		self.tree_data = tree_data
		self.k = k
		self.p = p
		self.pp_quantiles = pp_quantiles
		threading.Thread.__init__ ( self )
		
	def run ( self ):
		self.pp_quantiles.append(quantile_test( self.tree_data,self.k,self.p,debug ))

###############  Begin Test Statistic Functions ###############

def calcAndSortTreeDists(tree_data):
	"""
	A function to calculate RF distances from a set of trees and sort them.  This sorted
	vector will then be used when calculating quantile-based values or interquantile
	distances.
	"""
	try:
		# Create initial variables, including vector to hold pairwise RF values
		rf_dists = [ ]
		trees = tree_data.tree_lists[0]
		total_trees = len(trees)
		
		# Calculate and store all pairwise RF values
		# SLOW! Could approximate (and speed up) by sampling a large number of trees (but not all)?
		while len(trees) > 1:
			for i in range(len(trees))[:len(trees)-1]:
				# Calculates RFs b/w first tree and all others
				rf_dists.append(treecalc.symmetric_difference(trees[0],trees[i+1]))
			trees.pop(0);	# Removes first tree
		
		# Orders RF values (implemented in v0.94)
		rf_dists.sort()

		return rf_dists

	except:
		print "Problem calculating and ordering RF distances..."
		if debug:
			traceback.print_exc()
		sys.exit(1)

def quantile_test(rf_dists,k,p):
	"""
	Calculates positions of a set of k-th p-quantiles in ordered vector of RF distances 
	from one posterior distribution.
	"""
	
	try:
		quantile = []
		for i in range(len(k)):
			# Iterate and find k-th p-quantiles
			# 	See definition of k-th p-quantile position in my manuscript
			g = (len(rf_dists)*int(k[i])) % int(p[i])
			j = (len(rf_dists)*int(k[i])) / int(p[i])
			if g == 0:
				quantile.append((float(rf_dists[j])+float(rf_dists[j+1]))/2)
			else:
				quantile.append(float(rf_dists[j+1]))
		return quantile

	except:
		print "Problem calculating quantile-based test statistic..."
		if debug:
			traceback.print_exc()
		sys.exit(1)

def entropy_test(trees):
	"""
	Calculates entropy-based test statistic for one posterior distribution
	"""
	
	try:
		temp_tree_list = trees.tree_lists[0]
		total_trees = len(temp_tree_list)
		unique_topos = 0
		topo_freqs = [ ]
		while len(temp_tree_list) > 0: # Iterates over trees within a single list
			unique_topos += 1
			rf_dists = [ ]
				
			# Look for and record # of topologies identical to j
			for k in range(len(temp_tree_list)): # Iterates over all remaining trees
				rf_dists.append(treecalc.symmetric_difference(temp_tree_list[0],temp_tree_list[k]))
				
			if debug:
				print rf_dists
				
			# Appends frequency of new unique topology
			topo_freqs.append(float((rf_dists.count(0)))/float(total_trees)) 
				
			if debug:   
				print topo_freqs
				
			# Delete duplicate topologies
			while rf_dists.count(0) > 0:
				temp_tree_list.pop(rf_dists.index(0))
				rf_dists.remove(0)

		######  Calculates entropy test statistic for the posterior predictive data sets  #####

		#	Calculates change in entropy from the prior to the posterior
		#	-- Assuming a uniform prior on topologies, change in entropy can be calculated as:
		#		T(X) = (sum across topologies i=1 to N: post_i * ln(post_i)) - ln(prior_i)
			
		entropy = 0
		for j in range(len(topo_freqs)):
			entropy = entropy + (topo_freqs[j]*math.log(topo_freqs[j]))
		
		no_taxa = len(trees.taxon_sets[0])		
		if debug:
			print "# Taxa: %d" % no_taxa
				
		no_topologies = int(factorial(2*no_taxa-4))/int((math.pow(2,(no_taxa-2)) * factorial(no_taxa-2)))
		if debug:
			print "# Topologies: %d" % no_topologies
				
		# Could remove subtraction by final constant if running into overflow problems in 
		#	calculating the number of topologies.
		entropy = entropy - math.log(float(1)/float(no_topologies))
		
		if debug:
			print """
			Entropy = %f
			""" % entropy
	
		return(entropy)
	
	except:
		print "Problem calculating entropy-based test statistic..."
		if debug:
			traceback.print_exc()
		sys.exit(1)

	

def iqr_test(rf_dists):
	"""
	Calculates interquartile distance (1st to 3rd quartile) for one posterior distribution
	"""
	try:
	
		# Find 1st quartile
		# 	See definition of k-th p-quantile position in my manuscript
		g = (len(rf_dists)*1) % 4
		j = (len(rf_dists)*1) / 4
		if g == 0:
			first_quant = (float(rf_dists[j])+float(rf_dists[j+1]))/2
		else:
			first_quant = float(rf_dists[j+1])

		# Find 3rd quartile
		# 	See definition of k-th p-quantile position in my manuscript
		g = (len(rf_dists)*3) % 4
		j = (len(rf_dists)*3) / 4
		if g == 0:
			third_quant = (float(rf_dists[j])+float(rf_dists[j+1]))/2
		else:
			third_quant = float(rf_dists[j+1])	
	
		iqr = third_quant - first_quant
		
		return iqr
		
	except:
		print """
	Problem calculating interquartile range test statistic...
		"""
		if debug:
			traceback.print_exc()
		sys.exit(1)

def partition_test(tree_data,taxa):
	"""
	Calculates entropy of posterior distribution with two categories: (i) trees with a split or
	(ii) trees without a split.
	"""
	
	try:
		trees = tree_data.tree_lists[0]
		p = trees.frequency_of_split(labels=taxa)
		if p != 1 and p != 0:
			ent_stat = -(p*math.log(p)+((1-p)*math.log((1-p))))
		else:
			ent_stat = 0
		return ent_stat
	
	except:
		print """
	Problem calculating bipartition frequency test statistic...
		"""
		if debug:
			traceback.print_exc()
		sys.exit(1)	

"""
	# QUARTET
def bit_mask(focal,ref):
	"""
	#Used by partition_test
	
	#Generic function to create bitmask corresponding to a particular bipartition. "focal" is a list 
	#of tips on one side of a bipartition, while "ref" is a list of all tips in a tree
"""
	bitmask = []
	for i in ref:
		if i == ref[0]:
			bitmask.append(1)
		elif (i in focal) == (ref[0] in focal):
			bitmask.append(1)
		else:
			bitmask.append(0)
	return bitmask
"""

def treelength_test(tree_data):
	"""
	Calculates mean treelength for one posterior distribution
	"""
	try:
		sum = 0.0
		trees = tree_data.tree_lists[0]
		for tr in trees:
			sum += tr.length()
		mean_tl = sum/len(trees)
		
		return mean_tl
	
	except:
		print """
	Problem calculating mean treelength test statistic...
		"""
		if debug:
			traceback.print_exc()
		sys.exit(1)

def mean(vals):
	"""
	Calculates the mean of a list of numbers.
	"""
	sum = 0.0
	for i in vals:
		sum += i
	return float(sum)/float(len(vals))

def variance(vals):
	"""
	Calculates the variance of a list of numbers.
	"""
	sumSquares = 0.0
	meanVal = mean(vals)
	for i in vals:
		sumSquares += math.pow((i-meanVal),2)
	return sumSquares/float(len(vals))

def tlVar_test(tree_data):
	"""
	Calculates the variance in treelength across one posterior distribution
	"""
	try:
		tls = [ ]
		trees = tree_data.tree_lists[0]
		for tr in trees:
			tls.append(tr.length())
		return variance(tls)
	
	except:
		print """
	Problem calculating treelength variance test statistic...
		"""
		if debug:
			traceback.print_exc()
		sys.exit(1)

	
###############  End Test Statistic Functions ###############
	
	
	
###############  Begin Utility Functions ###############
	
def factorial(x):
	"""
	Home-cooked function to calculate a factorial
	"""
	factorial = 1
	for i in range(x):
		factorial = factorial * (i+1)
	return factorial
	
def quickSort(list):
	"""
	quickSort algorithm for ordering values in a list
	Code taken from: http://en.literateprograms.org/Quicksort_(Python)
	
	MUCH slower than python's built-in list sorting
	"""
	if list == [ ]:
		return [ ]
	else:
		pivot = list[0]
		lesser = quickSort([x for x in list[1:] if x < pivot])
		greater = quickSort([x for x in list[1:] if x >= pivot])
		return lesser + [pivot] + greater


###############  End Utility Functions ###############


def parse_in(opts,args):
	"""
	Responsible for interpreting command-line input
	"""
	
	# Status update
	print """
	Parsing input...
	"""
	
	# Debugging -- make sure that opts and args are correct
	# print opts
	# print args
	
	# Default values (all statistics turned off by default)
	quantile = False
	quant_focal_breaks = []
	quant_bin_nos = []
	ent = False
	inter = False
	part = False
	part_list = []
	length = False
	manual = False
	burn = 0
	mrc = True
	l = False
	u = False
	t = False
	out = False
	outfile = "amp_out.txt"
	v = False
	like = False
	likeConstraints = 1
	bf = False
	bfConstraints = 1
	mbSSbf = False
	mbSSbfConstraints = 1
	tlVar = False
	
	# Sets up a dictionary where options are keys and arguments are values
	opt_dict = {}
	for o, a in opts:
		opt_dict[o] = a
		# print (o,a)	# For debugging flags
		
	# Mild dummy-proofing for command-line input
	if len(opt_dict) == 0: # If no input, print out usage
		if len(args) == 0:
			print "	No command-line input detected!"
			usage()
	if len(opt_dict) == 0: # Error message, b/c no test statistics to calculate
		print """
		Nothing to do if you don't specify test statistics! Exiting...
		"""
		usage()
	if len(args) == 3:
		repNum = args[2]
	else:
		repNum = 0
	if len(args) < 2:
		print"""
		Some required command-line argument is missing! Exiting...
		"""
		usage()
		
	# Parsing optional input
	if opt_dict.has_key("-q"):
		quantile = True
		try:
			quant_string = opt_dict["-q"]
			if quant_string.find(",") == -1:
				raise Exception;
			if (len(quant_string.split(",")) % 2 != 0): # Nums not in pairs
				raise Exception
			quant_vals = quant_string.split(",")
			for i in range(len(quant_vals)/2):
				quant_bin_nos.append(quant_vals.pop())
				quant_focal_breaks.append(quant_vals.pop())
		except:
			print """
	Problem with values given for quantile-based test statistic. Exiting...
			"""
			if debug:
				traceback.print_exc()
			sys.exit(1)
		for i in range(len(quant_focal_breaks)):
			if int(quant_focal_breaks[i]) >= int(quant_bin_nos[i]):
				print """
	All focal quantiles must be smaller than the corresponding overall number of bins! Exiting...
				"""
				sys.exit(1)
	ent = opt_dict.has_key('-e')
	inter = opt_dict.has_key('-i')
	if opt_dict.has_key('-p'):
		part = True
		part_list = opt_dict['-p']
		part_list = part_list.split(",")
	length = opt_dict.has_key('-T')
	if opt_dict.has_key('-m'):
		manual = True
		burn = opt_dict['-m']
		mrc = False
	l = opt_dict.has_key('-l')
	u = opt_dict.has_key('-u')
	t = opt_dict.has_key('-t')
	if opt_dict.has_key('-o'):
		out = True
		if opt_dict['-o'] != "":
			outfile=opt_dict['-o']
	v = opt_dict.has_key('-v')
	if opt_dict.has_key('-L'):
		like = True
		likeConstraints = opt_dict['-L']
	if opt_dict.has_key('-b'):
		bf = True
		bfConstraints = opt_dict['-b']
	if opt_dict.has_key('-B'):
		mbSSbf = True
		mbSSbfConstraints = opt_dict['-B']
	tlVar = opt_dict.has_key('-V')
	
	"""
	Need to return (bool[q], list of ints[focal breaks], list of ints[total #'s of bins],
					bool[e],
					bool[i],
					bool[p], list[p -- taxon names],
					bool[T],
					bool[m], int[m - burnin],
					bool[c],
					bool[l],
					bool[u],
					bool[t],
					bool[o],
					string[outfile],
					bool[v],
					string[basename],
					int[#.datasets],
					int[#.reps],
					bool[L], int[#.like.constraints],
					bool[bf], int[#.bf.constraints],
					bool[mbSSbf], int[#.mbSSbf.constraints],
					bool[V])
	"""
	
									# index
	return (quantile,				# 0
			quant_focal_breaks,		# 1
			quant_bin_nos,			# 2
			ent,					# 3
			inter,					# 4
			part,					# 5
			part_list,				# 6
			length,					# 7
			manual,					# 8
			burn,					# 9
			mrc,					# 10
			l,						# 11
			u,						# 12
			t,						# 13
			out,					# 14
			outfile,				# 15
			v,						# 16
			args[0],				# 17
			args[1],				# 18
			repNum, 				# 19
			like,					# 20
			likeConstraints,		# 21 
			bf,						# 22
			bfConstraints,			# 23
			mbSSbf,					# 24
			mbSSbfConstraints,		# 25
			tlVar)					# 26
	
	
def usage():
	"""
	Prints a nice summary of the program usage.
	"""
	print """
	
	Proper program usage:
	
	Usage: python amp0.98.py [-q##eip#TVL#b#B#] [-m#c] [-lut] [-o#v] basename #.datasets #.reps
	
	[ ] -- square brackets denote optional arguments
	
	# -- some argument is needed if this option is invoked (details below)
	
	basename: The portion of the filename common to all analysis files.  Tree files should
			  have the following name structure:
			  
			  <basename>_<treefile#>_<rep#>.t, 
			  
			  where the parts within < and > should be substituted as necessary. This name can 
			  include a path if files are not in the same directory as this script. The treefiles 
			  resulting from analysis of the original (empirical) data should be named: 
			  
			  <basename>_emp_<rep#>.t  
			  
			  The rep# will not be used if only the bipartition-specific likelihood ratio test 
			  statistic is selected.  When using the bipartition-specific likelihood ratio test 
			  statistic, output ML tree files from constrained searches (only Garli supported right 
			  now) should have this structure: 
			  
			  <basename>_<treefile#>_bp<constraint#>.<pos/neg>.constraint.best.tre, 
			  
			  substituting pos or neg as appropriate in the name to specify a positive or negative 
			  constraint.  When using the bipartition-specific Bayes factor test statistic estimated 
			  using model-switch thermodynamic integration (e.g., in-house software written by JMB), 
			  output files for Bayes factors should be named:
			  
			  <basename>_<treefile#>_bp<constraint#>.bf.out
			  
			  and should contain the following line:
			  
			  Estimated ln(Bayes Factor) = <ln(BF) value>
			  
			  When using the bipartition-specific Bayes factor test based on marginal likelihoods
			  estimated with steppingstone sampling in MrBayes v3.2.1, two types of output files will
			  be used, based on runs positively or negatively constraining a particular branch. These
			  output files should have the form:
			  
			  <basename>_<treefile#>.con<constraint#>.<pos/neg>.log
			  
			  The average marginal likelihood (taken from the line containing the word "Mean") in each
			  log file will be used to calculate the Bayes factor.
						 
	#.datasets: The total number of posterior predictive datasets that were analyzed.
	
	#.reps: The number of replicate analyses run for each posterior predictive dataset.  This value
			will be ignored if only the bipartition-specific likelihood ratio test statistic has 
			been chosen.
	
	Test Statistic Options [-q##eip#TVL#b#B#]:
		q: quantile-based test statistic
		    [Requires at least two numbers: k and p. These correspond to the k-th p-quantile. If only 
		    two, these numbers should be given as "-q k,p". If you wish to calculate multiple such statistics,
		    continue adding pairs of numbers separated by commas.  For example, "-q k1,p1,k2,p2,k3,p3,...".]
		
		e: entropy-based test statistic
		
		i: interquartile range test statistic
		
		p: partition-specific entropy test statistic [Requires taxon names on one side of partition.]
		
		T: mean treelength test statistic
		
		V: treelength variance test statistic
		
		L: bipartition-specific likelihood ratio [Requires number of constraint trees used.]
		
		b: bipartition-specific Bayes factor estimated via model-switch TI
		   [Requires number of constraint trees used.]
		
		B: bipartition-specific Bayes factor estimated via steppingstone sampling of marginal 
		   likelihoods for positively and negatively constrained searches in MrBayes v3.2.1
		   [Requires number of constraint trees used.]
		
	Burnin Options [-m#c]
		m: manual burn-in [Needs to be followed by an integer value]
		
		c: find burn-in values in MrConverge log files
		   [Requires log files named as: basename_dataset#.log]
		
	P-value Type [-lut]
		l: lower one-tailed p-value
		
		u: upper one-tailed p-value
		
		t: two-tailed p-value
	
	Output Options [-o#v]
		o: direct output to a file
		   [Requires an output file name]
		
		v: output test statistic values for all posterior predictive datasets
	
	Defaults:  No test statistics are turned on by default.
			   Default burn-in is 0.
			   Default quantile is the 9th 10-quantile.
			   Default output filename is pma_out.txt.
	"""
	sys.exit(0)



def sigint_handler(x,y):
	print """
	
	Program execution interrupted by the user.  Exiting..."""
	sys.exit(1)



if __name__ == "__main__":
	main(sys.argv[1:])
