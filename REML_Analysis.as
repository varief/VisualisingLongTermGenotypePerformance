title: Two-stage analysis
 year  25  !I				# Year
 trial 1387 !I				# Trial
 env 1362 !I				# Combination of years and locations
 gen  686 !I				# Genotype
 BLUE    				# BLUE from the first-stage analysis
 wt    					# Weight from the first-stage analysis (Cullis et.l. 1996)


Kgiv.giv !skip 1			# Inverse of the kinship matrix
Data.csv !skip 1  !WORKSPACE 8000 !NODISPLAY   !CONTINUE !MAXIT 100	# Input file contains factors, BLUE and weight



# Models for analysis
# ID model - no kinship and FA
#BLUE !WT wt ~ mu  !r   yea.gen env.gen !f trial mv		

# K model - kinship only
#BLUE !WT wt ~ mu  !r   yea.giv(gen,1) env.gen !f trial mv		

# KFA model - fit kinship and FA2 model
#BLUE !WT wt ~ mu  !r   diag(yea).giv(gen,1) env.gen !f trial mv		# DIAG model - provide initial value to run FA model				
#BLUE !WT wt ~ mu  !r   xfa(yea,1).giv(gen,1) env.gen !f trial mv	# FA1 model - kinship + FA1 - provide initial value to run FA2 model 			
#BLUE !WT wt ~ mu  !r   xfa(yea,2).giv(gen,1) env.gen !f trial mv	# Final KFA model	


# G and R structure
1 1 0
0 0 0  !s2== 0.44		# Fit a fixed weight


# Obtain BLUP for each genotype in each year
predict yea.gen 1:100 !IGNORE mu trial
predict yea.gen 101:200 !IGNORE mu trial
predict yea.gen 201:300 !IGNORE mu trial
predict yea.gen 301:400 !IGNORE mu trial
predict yea.gen 401:500 !IGNORE mu trial
predict yea.gen 501:600 !IGNORE mu trial
predict yea.gen 601:686 !IGNORE mu trial

