#!/usr/bin/env python

#Author: Jose Lugo-Martinez
#Seeded from Jun Ding's alignment code for gene expression profiles
#File: getAlignmentsIBD.py
#Date: October 31, 2017
#Advisor Profs. Ziv Bar-Joseph and Giri Narasimhan
#Description: This function aligns temporal metagenomic samples using a linear time warping algorithm over relative abundance across multiple taxa.
#The program reports a set of taxa that best agrees with the global alignment.

#Last Modified: October 08, 2018

#Example calls: python getAlignmentsIBD.py human_ibd_microbiota.tsv True human_idb_microbiota_taxon_pairwise_all_subject_ranking.tsv

import sys, copy, math, random
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
from scipy.stats import pearsonr
from scipy.stats import spearmanr
from scipy.stats import kendalltau

#Parameters
METABOLITES_OFFSET = 3 #Index offset in dataset where taxon data starts 
MINIMUN_NUMBER_MEASURED_TIMEPOINTS = 5 #Minimum number of required measured timepoints per subject/sample
TOP_K = 25 #Maximum number features (e.g., taxa) for temporal alignment 
OVERLAP_THRESHOLD = 0.50 #Minimum overlap allowed between measured points of reference sample
UPPER_BOUND = 1.0 #Maximum value in relation to relative abundance for bounding continuous representation range
LOWER_BOUND = 0.0 #Minimum value in relation to relative abundance for bounding continuous representation range
SAMPLING_RATE = 1.0 #1, 2, 3, 5, 7, 14
SCALING_FACTOR = 10.0

SUBJECT_MAP = {'C3009': 9, 'C3008': 8, 'C3003': 3, 'C3002': 2, 'C3001': 1, 'C3007': 7, 'C3006': 6, 'E5023': 44, 'C3004': 4, 'P6009': 116,
	       'M2028': 86, 'C3015': 14, 'M2021': 81, 'M2024': 82, 'M2025': 83, 'M2026': 84, 'M2027': 85, 'H4012': 53, 'H4013': 54, 'H4010': 51,
	       'H4011': 52, 'H4016': 57, 'H4017': 58, 'H4014': 55, 'H4015': 56, 'M2103': 114, 'H4018': 59, 'H4019': 60, 'M2048': 92, 'C3005': 5,
	       'M2047': 91, 'E5002': 35, 'E5022': 43, 'H4035': 70, 'M2059': 94, 'M2058': 93, 'H4028': 66, 'H4027': 65, 'H4024': 64, 'H4023': 63,
	       'H4022': 62, 'H4020': 61, 'E5009': 40, 'E5008': 39, 'E5003': 36, 'H4044': 76, 'E5001': 34, 'H4040': 73, 'H4043': 75, 'E5004': 37,
	       'H4038': 71, 'H4039': 72, 'C3016': 15, 'C3017': 16, 'C3010': 10, 'C3011': 11, 'C3012': 12, 'C3013': 13, 'H4030': 67, 'H4031': 68,
	       'H4032': 69, 'M2042': 90, 'C3019': 17, 'M2041': 89, 'P6012': 118, 'M2079': 105, 'M2072': 102, 'M2071': 101, 'M2077': 104, 'M2075': 103,
	       'P6028': 126, 'M2091': 112, 'P6025': 125, 'P6024': 124, 'M2097': 113, 'M2064': 97, 'P6010': 117, 'M2067': 98, 'M2060': 95, 'M2061': 96,
	       'M2068': 99, 'M2069': 100, 'H4045': 77, 'P6013': 119, 'E5013': 41, 'P6014': 120, 'E5019': 42, 'P6038': 130, 'P6016': 121, 'M2082': 107,
	       'M2083': 108, 'M2081': 106, 'M2086': 111, 'P6037': 129, 'M2084': 109, 'E5006': 38, 'P6018': 123, 'P6017': 122, 'H4042': 74, 'M2010': 79,
	       'M2014': 80, 'C3029': 25, 'C3028': 24, 'C3024': 22, 'C3027': 23, 'C3021': 19, 'C3020': 18, 'C3023': 21, 'C3022': 20, 'M2008': 78,
	       'P6033': 127, 'C3036': 32, 'C3037': 33, 'C3034': 30, 'C3035': 31, 'C3032': 28, 'C3033': 29, 'C3030': 26, 'C3031': 27, 'M2085': 110,
	       'P6035': 128, 'M2039': 88, 'M2034': 87, 'H4001': 45, 'P6005': 115, 'H4004': 46, 'H4007': 48, 'H4006': 47, 'H4009': 50, 'H4008': 49}

class timepoint:
	#constructor ---
	def __init__(self, offsetID, ID, taxaNames, abundanceValuesPerTaxa, splineParametersPerTaxa):
		self.offsetID = offsetID #Use to offset timepoints from different timescales (not always necessary). 
		self.ID = ID
		relativeAbundance = {}
		splineParameters = {}
		for taxaIndex in xrange(len(taxaNames)):
			relativeAbundance[taxaNames[taxaIndex]] = abundanceValuesPerTaxa[taxaIndex]
			splineParameters[taxaNames[taxaIndex]] = splineParametersPerTaxa[taxaIndex]
		self.relativeAbundance = relativeAbundance
		self.splineParameters = splineParameters

class taxa:
	#constructor----
	def __init__(self, ID, timePoints, relativeAbundance, splineParameters):
		self.ID = ID
		self.timePoints = timePoints
		self.relativeAbundance = relativeAbundance
		self.splineParameters = splineParameters
		
	def getMean(self):
		if len(self.relativeAbundance) < 1:
			return 0.0
		
		return sum(self.relativeAbundance) / float(len(self.relativeAbundance))
	
	def getVariance(self):
		variance = 0.0
		meanAbundanceValue = sum(self.relativeAbundance) / float(len(self.relativeAbundance))
		for currentAbundanceValue in self.relativeAbundance:
			variance += (currentAbundanceValue - meanAbundanceValue)**2
		variance = variance / float(len(self.relativeAbundance))
		self.variance = variance

		return variance

	def getMeanSpline(self):
		abundance = interpolate.splev(self.timePoints, self.splineParameters)

		return sum(abundance) / float(len(abundance))

	def getVarianceSpline(self):
		variance = 0.0
		abundance = interpolate.splev(self.timePoints, self.splineParameters)
		meanAbundanceValue = sum(abundance) / float(len(abundance))
		for currentAbundanceValue in abundance:
			variance += (currentAbundanceValue - meanAbundanceValue)**2
		variance = variance / float(len(abundance))
		self.variance = variance

		return variance

def buildTaxon(taxonSample, taxonSplines):
	timepointHeaders = taxonSample[0][1:]
	taxon = []
	for currTaxa in taxonSample[1:]:
		currentTaxa = taxa(currTaxa[0], timepointHeaders, [float(relativeAbundance) for relativeAbundance in currTaxa[1:]], taxonSplines[currTaxa[0]])
		taxon.append(currentTaxa)

	return taxon

def filterTaxon(taxonReferenceSample, taxonCurrentSample, useSplines):
	outTaxonReferenceSample = []
	outTaxonCurrentSample = []
	taxonCurrentSampleIDs = [taxaCurrentSample.ID for taxaCurrentSample in taxonCurrentSample]
	for currentTaxaReferenceSample in taxonReferenceSample:
		if currentTaxaReferenceSample.ID in taxonCurrentSampleIDs:
			currentTaxaIndexCurrentSample = taxonCurrentSampleIDs.index(currentTaxaReferenceSample.ID)
			currentTaxaCurrentSample = taxonCurrentSample[currentTaxaIndexCurrentSample]
			if useSplines == True:
				meanTaxaReferenceSample = currentTaxaReferenceSample.getMeanSpline()
				varianceTaxaReferenceSample = currentTaxaReferenceSample.getVarianceSpline()
				meanTaxaCurrentSample = currentTaxaCurrentSample.getMeanSpline()
				varianceTaxaCurrentSample = currentTaxaCurrentSample.getVarianceSpline()
			else:
				meanTaxaReferenceSample = currentTaxaReferenceSample.getMean()
				varianceTaxaReferenceSample = currentTaxaReferenceSample.getVariance()
				meanTaxaCurrentSample = currentTaxaCurrentSample.getMean()
				varianceTaxaCurrentSample = currentTaxaCurrentSample.getVariance()
			#NOTE: This removes taxa whose relative abundance profiles are either (1) too low (mean < 1%), or (2) unchanged (variance = 0) in at least one of the two time-series samples.
			if meanTaxaReferenceSample >= 0.01 and varianceTaxaReferenceSample > 0.0 and meanTaxaCurrentSample >= 0.01 and varianceTaxaCurrentSample > 0.0:
				outTaxonReferenceSample.append([varianceTaxaReferenceSample, currentTaxaReferenceSample])
				outTaxonCurrentSample.append([varianceTaxaCurrentSample, currentTaxaCurrentSample])
	outTaxonReferenceSample.sort(reverse=True)
	outTaxonCurrentSample.sort(reverse=True)
	
	outTaxonReferenceSample = [taxaReferenceSample[1] for taxaReferenceSample in outTaxonReferenceSample]
	outTaxonCurrentSample = [taxaCurrentSample[1] for taxaCurrentSample in outTaxonCurrentSample]
	taxonCurrentSampleIDs = [taxaCurrentSample.ID for taxaCurrentSample in outTaxonCurrentSample]
	filteredTaxonCurrentSample = []
	for currentTaxaReferenceSample in outTaxonReferenceSample:
		currentTaxaIndexCurrentSample = taxonCurrentSampleIDs.index(currentTaxaReferenceSample.ID)
		filteredTaxonCurrentSample.append(outTaxonCurrentSample[currentTaxaIndexCurrentSample])
	filteredTaxonReferenceSample = outTaxonReferenceSample[0:TOP_K]
	filteredTaxonCurrentSample = filteredTaxonCurrentSample[0:TOP_K]

	return [filteredTaxonReferenceSample, filteredTaxonCurrentSample]

def buildTimepointsProfile(taxonSample, sampleInfo, useSplines):
	sampleTimepoints = taxonSample[0].timePoints
	taxonNames = [taxaSample.ID for taxaSample in taxonSample]
	taxonAbundances = [taxaSample.relativeAbundance for taxaSample in taxonSample]
	taxonSplineParameters = [taxaSample.splineParameters for taxaSample in taxonSample]
	taxonSampleTimepoints = []
	firstSample = float(sampleInfo[0])
	lastSample = float(sampleInfo[1])
	#Add week of first sample obtained as a timepoint via continuous representation (if enabled)
	if useSplines and not (firstSample in sampleTimepoints):
		abundances = []
		for taxaIndex in xrange(len(taxonNames)):
			abundances.append(interpolate.splev(firstSample, taxonSplineParameters[taxaIndex]))
		currTimepoint = timepoint(firstSample, firstSample, taxonNames, abundances, taxonSplineParameters)
		currTimepoint.offsetID = float(firstSample) #- float(firstSample) #Alternatively, one can just hard-code 0.0
		currTimepoint.ID = float(firstSample)
		taxonSampleTimepoints.append(currTimepoint)
	else:
		firstSample = copy.copy(sampleTimepoints[0])
	#Process originally measured timepoints	
	for timepointIndex in xrange(len(sampleTimepoints)):
		currTimepoint = timepoint(sampleTimepoints[timepointIndex], sampleTimepoints[timepointIndex], taxonNames, [taxaAbundances[timepointIndex] for taxaAbundances in taxonAbundances], taxonSplineParameters)
		currTimepoint.offsetID = float(sampleTimepoints[timepointIndex]) #- float(firstSample)
		currTimepoint.ID = float(sampleTimepoints[timepointIndex])
		taxonSampleTimepoints.append(currTimepoint)
	#Add week of last sample obtained as a timepoint via continuous representation (if enabled)
	if useSplines and not (lastSample in sampleTimepoints):
		abundances = []
		for taxaIndex in xrange(len(taxonNames)):
			abundances.append(interpolate.splev(lastSample, taxonSplineParameters[taxaIndex]))
		currTimepoint = timepoint(lastSample, lastSample, taxonNames, abundances, taxonSplineParameters)
		currTimepoint.offsetID = float(lastSample) #- float(lastSample) #Alternatively, one can just hard-code 0.0
		currTimepoint.ID = float(lastSample)
		taxonSampleTimepoints.append(currTimepoint)
	else:
		lastSample = copy.copy(sampleTimepoints[-1])
	
	return taxonSampleTimepoints

def compareTimepoint(timepointReferenceSample, timepointCurrentSample, a, b, useSplines, taxonWorkingSet, method='ssd'):
	abundanceValuesReferenceSample = []
	abundanceValuesCurrentSample = []
	for currTaxaReferenceSample in timepointReferenceSample.relativeAbundance:
		if currTaxaReferenceSample in taxonWorkingSet:
			currTaxaCurrentSample = taxonWorkingSet[currTaxaReferenceSample]
			if currTaxaCurrentSample in timepointCurrentSample.relativeAbundance:
				if useSplines == True:
					currTaxaReferenceSampleAbundanceValue = interpolate.splev(timepointReferenceSample.offsetID, timepointReferenceSample.splineParameters[currTaxaReferenceSample])
					currTaxaCurrentSampleAbundanceValue = interpolate.splev(warpFunctionInverse(a, b, timepointCurrentSample.ID), timepointCurrentSample.splineParameters[currTaxaCurrentSample])
				else:
					currTaxaReferenceSampleAbundanceValue = timepointReferenceSample.relativeAbundance[currTaxaReferenceSample] 
					currTaxaCurrentSampleAbundanceValue = timepointCurrentSample.relativeAbundance[currTaxaCurrentSample]
				abundanceValuesReferenceSample.append(currTaxaReferenceSampleAbundanceValue)
				abundanceValuesCurrentSample.append(currTaxaCurrentSampleAbundanceValue)

	abundanceValuesReferenceSample = truncateAbundanceValues(abundanceValuesReferenceSample)
	abundanceValuesCurrentSample = truncateAbundanceValues(abundanceValuesCurrentSample)	
	if method == 'pearson':
		value = pearsonr(abundanceValuesReferenceSample, abundanceValuesCurrentSample)[0]
	elif method == 'spearman':
		value = spearmanr(abundanceValuesReferenceSample, abundanceValuesCurrentSample)[0]
	elif method == 'kendalltau':
		value = kendalltau(abundanceValuesReferenceSample, abundanceValuesCurrentSample)[0]
	else:
		#Get sum of squared differences
		value = getSSD(abundanceValuesReferenceSample, abundanceValuesCurrentSample)
		
	return value

def getAgreementPerTimepoint(timepointsListReferenceSample, timepointsListCurrentSample, a, b, useSplines, taxonWorkingSet, method='ssd'):
	P = []
	for currTimepointCurrentSample in timepointsListCurrentSample:
		PI = []
		for currTimepointReferenceSample in timepointsListReferenceSample:
			pij = compareTimepoint(currTimepointReferenceSample, currTimepointCurrentSample, a, b, useSplines, taxonWorkingSet, method)
			PI.append(pij)
		P.append(PI)

	return P

def warpFunction(a, b, s, warpType='linear'):
	if warpType == 'exponential':
		return np.exp((s - b) / a)
	else:
		return (s - b) / a

def warpFunctionInverse(a, b, t, warpType='linear'):
	if warpType == 'exponential':
		return a * np.log(t) + b
	else:
		return (a * t) + b

def getAlignmnetError(a, b, alpha, beta, timepointsListReferenceSample, timepointsListCurrentSample, taxonWeights, useSplines):
	timepointsReferenceSample = [timepointReferenceSample.offsetID for timepointReferenceSample in timepointsListReferenceSample]
	timepointsReferenceSampleSplineParameters = [timepointReferenceSample.splineParameters for timepointReferenceSample in timepointsListReferenceSample]
	timepointsCurrentSample = [timepointCurrentSample.offsetID for timepointCurrentSample in timepointsListCurrentSample]
	timepointsCurrentSampleSplineParameters = [timepointCurrentSample.splineParameters for timepointCurrentSample in timepointsListCurrentSample]
	filteredTaxon = timepointsListReferenceSample[0].relativeAbundance.keys()
	alignmentErrorPerTaxa = {}
	for currTaxa in filteredTaxon:
		alignmentErrorPerTaxa[currTaxa] = 0.0
	if useSplines == True:
		timepointsReferenceSample = np.arange(alpha, (beta + 1.0), SAMPLING_RATE)
	referenceSampleSplineParameters = timepointsReferenceSampleSplineParameters[0]
	currentSampleSplineParameters = timepointsCurrentSampleSplineParameters[0]
	for currentTimepoint in xrange(len(timepointsReferenceSample)):
		timepointReferenceSample = timepointsReferenceSample[currentTimepoint] #reference timpepoint s according to Bar-Joseph et al. (2003)
		if timepointReferenceSample < alpha or timepointReferenceSample > beta:
			continue
		timepointCurrentSampleTransformed = warpFunction(a, b, timepointReferenceSample) #T(s) according to Bar-Joseph et al. (2003)
		for currTaxa in filteredTaxon:
			if useSplines == True:
				relativeAbundanceTimepointReferenceSample = interpolate.splev(timepointReferenceSample, referenceSampleSplineParameters[currTaxa])
				relativeAbundanceTimepointCurrentSample = interpolate.splev(timepointCurrentSampleTransformed, currentSampleSplineParameters[currTaxa])
##			else:
##				relativeAbundanceTimepointReferenceSample = timepointReferenceSample.relativeAbundance[currTaxa] 
##				relativeAbundanceTimepointCurrentSample = timepointCurrentSample.relativeAbundance[currTaxa] #Need to map timepointCurrentSampleTransformed to closest point
			if relativeAbundanceTimepointReferenceSample > UPPER_BOUND:
				relativeAbundanceTimepointReferenceSample = UPPER_BOUND
			elif relativeAbundanceTimepointReferenceSample < LOWER_BOUND:
				relativeAbundanceTimepointReferenceSample = LOWER_BOUND
			if relativeAbundanceTimepointCurrentSample > UPPER_BOUND:
				relativeAbundanceTimepointCurrentSample = UPPER_BOUND
			elif relativeAbundanceTimepointCurrentSample < LOWER_BOUND:
				relativeAbundanceTimepointCurrentSample = LOWER_BOUND

			alignmentErrorPerTaxa[currTaxa] += ((relativeAbundanceTimepointReferenceSample - relativeAbundanceTimepointCurrentSample)**2)

	alignmentErrorTaxon = 0.0
	for taxaIndex in xrange(len(filteredTaxon)):
		currTaxa = filteredTaxon[taxaIndex]
		alignmentErrorTaxon += (alignmentErrorPerTaxa[currTaxa] / (beta - alpha)) * taxonWeights[taxaIndex]

	return [alignmentErrorTaxon, a, b]

def getOptimalMapping(timepointsListReferenceSample, timepointsListCurrentSample, taxonWeights, useSplines):
	ReferenceSampleT = [timepointReferenceSample.offsetID for timepointReferenceSample in timepointsListReferenceSample]
	CurrentSampleT = [timepointCurrentSample.offsetID for timepointCurrentSample in timepointsListCurrentSample]
	if useSplines == True:
		ReferenceSampleT = np.arange(ReferenceSampleT[0], (ReferenceSampleT[-1] + 1.0), SAMPLING_RATE)
		CurrentSampleT = np.arange(CurrentSampleT[0], (CurrentSampleT[-1] + 1.0), SAMPLING_RATE)
##		print ReferenceSampleT
##		print CurrentSampleT
	timepointReferenceSampleMin = min(ReferenceSampleT)
	timepointReferenceSampleMax = max(ReferenceSampleT)

	optimalAlignmentParameters = []
	for a in np.arange(0.01, 2.01, 0.01): #This paramater needs to be adjusted according to data set properties as well as warp function type
		for b in np.arange(-5.0, 5.5, 0.5): #This paramater needs to be adjusted according to data set properties as well as warp function type
			T = [warpFunction(a, b, timepointReferenceSample.offsetID) for timepointReferenceSample in timepointsListReferenceSample]
##			T_inverse = [warpFunctionInverse(a, b, timepointCurrentSample.offsetID) for timepointCurrentSample in timepointsListCurrentSample]
			timepointCurrentSampleMin = warpFunctionInverse(a, b, min(CurrentSampleT))
			timepointCurrentSampleMax = warpFunctionInverse(a, b, max(CurrentSampleT))
			alpha = max(timepointReferenceSampleMin, timepointCurrentSampleMin)
			beta = min(timepointReferenceSampleMax, timepointCurrentSampleMax)
			overlap =  (beta - alpha) / (timepointReferenceSampleMax - timepointReferenceSampleMin)
			if overlap > OVERLAP_THRESHOLD and alpha < beta:
				[alignmentError, a, b] = getAlignmnetError(a, b, alpha, beta, timepointsListReferenceSample, timepointsListCurrentSample, taxonWeights, useSplines)
				if len(optimalAlignmentParameters) == 0 or optimalAlignmentParameters[0] > alignmentError:
					optimalAlignmentParameters = [alignmentError, a, b, alpha, beta, overlap]

	return optimalAlignmentParameters

def getAlignmentAgreementScorePerTaxa(timepointsListReferenceSample, timepointsListCurrentSample, a, b, alpha, beta, taxonWeights, method='ssd'):
	filteredTaxon = timepointsListReferenceSample[0].relativeAbundance.keys()
	timepointsReferenceSample = [timepointReferenceSample.offsetID for timepointReferenceSample in timepointsListReferenceSample]
	timepointsReferenceSampleSplineParameters = [timepointReferenceSample.splineParameters for timepointReferenceSample in timepointsListReferenceSample]
	timepointsCurrentSampleSplineParameters = [timepointCurrentSample.splineParameters for timepointCurrentSample in timepointsListCurrentSample]

	timepointsReferenceSample = np.arange(alpha, (beta + 1.0), SAMPLING_RATE)
	timepointsCurrentSampleAligned = [warpFunction(a, b, timepointReferenceSample)  for timepointReferenceSample in timepointsReferenceSample]
	referenceSampleSplineParameters = timepointsReferenceSampleSplineParameters[0]
	currentSampleSplineParameters = timepointsCurrentSampleSplineParameters[0]
	alignmentAgreementScoresPerTaxa = {}
	for currTaxa in filteredTaxon:
		relativeAbundancesReferenceSample = interpolate.splev(timepointsReferenceSample, referenceSampleSplineParameters[currTaxa])
		relativeAbundancesReferenceSample = truncateAbundanceValues(relativeAbundancesReferenceSample)
		relativeAbundancesCurrentSampleAligned = interpolate.splev(timepointsCurrentSampleAligned, currentSampleSplineParameters[currTaxa])
		relativeAbundancesCurrentSampleAligned = truncateAbundanceValues(relativeAbundancesCurrentSampleAligned)

##		if method == 'pearson':
##			alignmentScore = pearsonr(relativeAbundancesCurrentSampleAligned, relativeAbundancesReferenceSample)[0]
##		elif method == 'kendall':
##			alignmentScore = kendalltau(relativeAbundancesCurrentSampleAligned, relativeAbundancesReferenceSample)[0]
##		elif method == 'spearman':
##			alignmentScore = spearmanr(relativeAbundancesCurrentSampleAligned, relativeAbundancesReferenceSample)[0]
##		elif method == 'rsquared':
##			maxError = ((UPPER_BOUND - LOWER_BOUND)**2) * float(len(timepointsReferenceSample))
##			alignmentError = getSSD(relativeAbundancesReferenceSample, relativeAbundancesCurrentSampleAligned)
##			alignmentScore = 1.0 - (alignmentError / maxError)
##		else:
##			alignmentScore = getSSD(relativeAbundancesReferenceSample, relativeAbundancesCurrentSampleAligned)
		
		alignmentScorePearson = pearsonr(relativeAbundancesReferenceSample, relativeAbundancesCurrentSampleAligned)[0]
		alignmentScoreSpearman = spearmanr(relativeAbundancesReferenceSample, relativeAbundancesCurrentSampleAligned)[0]
		alignmentScoreSSD = getSSD(relativeAbundancesReferenceSample, relativeAbundancesCurrentSampleAligned)
		maxError = ((UPPER_BOUND - LOWER_BOUND)**2) * float(len(timepointsReferenceSample))
		alignmentScoreRsquared = 1.0 - (alignmentScoreSSD / maxError)
##		alignmentAgreementScoresPerTaxa[currTaxa] = [alignmentScorePearson, alignmentScoreSpearman, alignmentScoreSSD, alignmentScoreRsquared]
		alignmentAgreementScoresPerTaxa[currTaxa] = [alignmentScoreSSD]

	return alignmentAgreementScoresPerTaxa

def plotAlignment(timepointsListReferenceSample, timepointsListCurrentSample, a, b, alpha, beta, taxaName, referenceSampleID, currentSampleID, taxaAlignmentScore):
	timepointsReferenceSample = [timepointReferenceSample.offsetID for timepointReferenceSample in timepointsListReferenceSample]
	timepointsReferenceSampleSplineParameters = [timepointReferenceSample.splineParameters for timepointReferenceSample in timepointsListReferenceSample]
	timepointsCurrentSample = [timepointCurrentSample.ID for timepointCurrentSample in timepointsListCurrentSample]
	timepointsCurrentSampleSplineParameters = [timepointCurrentSample.splineParameters for timepointCurrentSample in timepointsListCurrentSample]
	referenceSampleSplineParameters = timepointsReferenceSampleSplineParameters[0]
	currentSampleSplineParameters = timepointsCurrentSampleSplineParameters[0]

	timepointsReferenceSampleOriginal = np.arange(timepointsReferenceSample[0], (timepointsReferenceSample[-1] + 1.0), SAMPLING_RATE)
	timepointsCurrentSampleOriginal = np.arange(timepointsCurrentSample[0], (timepointsCurrentSample[-1] + 1.0), SAMPLING_RATE)
	relativeAbundancesReferenceSampleOriginal = interpolate.splev(timepointsReferenceSampleOriginal, referenceSampleSplineParameters[taxaName])
	relativeAbundancesReferenceSampleOriginal = truncateAbundanceValues(relativeAbundancesReferenceSampleOriginal)
	relativeAbundancesCurrentSampleOriginal = interpolate.splev(timepointsCurrentSampleOriginal, currentSampleSplineParameters[taxaName])
	relativeAbundancesCurrentSampleOriginal = truncateAbundanceValues(relativeAbundancesCurrentSampleOriginal)
	
	timepointsReferenceSample = np.arange(alpha, (beta + 1.0), SAMPLING_RATE)
	timepointsCurrentSample = np.arange(alpha, (beta + 1.0), SAMPLING_RATE)
	timepointsCurrentSampleAligned = [warpFunction(a, b, timepointReferenceSample) for timepointReferenceSample in timepointsReferenceSample]
	timepointsCurrentSampleInverse = [warpFunctionInverse(a, b, timepointCurrentSample) for timepointCurrentSample in timepointsCurrentSample]
	relativeAbundancesReferenceSample = interpolate.splev(timepointsReferenceSample, referenceSampleSplineParameters[taxaName])
	relativeAbundancesReferenceSample = truncateAbundanceValues(relativeAbundancesReferenceSample)
	relativeAbundancesCurrentSample = interpolate.splev(timepointsCurrentSample, currentSampleSplineParameters[taxaName])
	relativeAbundancesCurrentSample = truncateAbundanceValues(relativeAbundancesCurrentSample)
	relativeAbundancesCurrentSampleAligned = interpolate.splev(timepointsCurrentSampleAligned, currentSampleSplineParameters[taxaName])
	relativeAbundancesCurrentSampleAligned = truncateAbundanceValues(relativeAbundancesCurrentSampleAligned)

	fig = plt.figure() #plt.figure(figsize=(3, 6))
	plt.plot(timepointsReferenceSampleOriginal, relativeAbundancesReferenceSampleOriginal, '--b', timepointsReferenceSample, relativeAbundancesReferenceSample, '-b',
		 timepointsCurrentSampleOriginal, relativeAbundancesCurrentSampleOriginal, '--g', timepointsReferenceSample, relativeAbundancesCurrentSampleAligned, '-g')
#	title = 'Alignment of ' + str(taxaName) + ' for ' + referenceSampleID + ' to ' + currentSampleID + ' [a = ' + str(a) + ', b = ' + str(b) + ' | Alignment score [Pearson, Spearman, SSD, Rsquared]:' + str(taxaAlignmentScore) + ']'
	title = 'Alignment of ' + str(taxaName) + ' for ' + referenceSampleID + ' to ' + currentSampleID + ' (a = ' + str(a) + ', b = ' + str(b) + ' | Alignment interval: [' + str(alpha) + ', ' + str(beta) + '])'
	plt.title(title)
	plt.legend(['reference sample unaligned', 'reference sample aligned', 'candidate sample unaligned', 'candidate sample aligned'])
	plt.show()
#	fig.savefig(taxaName + '_' + referenceSampleID + '_vs_' + currentSampleID + '.png', dpi=fig.dpi)
	
	return 

def getSamples(dataFilename):
	try:
		#Open input file
		infile = open(dataFilename, "r")
	except(IOError), e:
		print "<<ERROR>> Unable to open the file", dataFilename, "\nThis program will be quiting now.", e
		sys.exit()

	headers = infile.readline().strip().split(',')
	metabolomeNames = copy.copy(headers[METABOLITES_OFFSET:])
	print "# of Metabolites", len(metabolomeNames)

	subjectIDs = []
	samplesPerSubject = {}
	samplesPerSubjectInfo = {}
	currSubjectSample = {}
	previousSubjectID = ''
	#Iterate over file
	for line in infile:
		line = line.strip()
		tokens = line.split(',')
		patientID = tokens[1]
		subjectID = str(SUBJECT_MAP[patientID])
		if not (subjectID in subjectIDs):
			subjectIDs.append(subjectID)
		if previousSubjectID != '' and previousSubjectID != subjectID:
			samplesPerSubject[previousSubjectID] = copy.copy(currSubjectSample)
			sampleSizePerSubject = len(currSubjectSample[metabolomeNames[0]])
			samplesPerSubjectInfo[previousSubjectID] = (weekOfLifeFirstSample, weekOfLifeLastSample, sampleSizePerSubject)
			currSubjectSample = {}
		weekSampleObtained = float(tokens[2])
		currentAbundancePerTaxa = copy.copy(tokens[METABOLITES_OFFSET:])
		for taxaIndex in xrange(len (metabolomeNames)):
			taxaName = metabolomeNames[taxaIndex]
			abundance = float(currentAbundancePerTaxa[taxaIndex].strip())
			if not (taxaName in currSubjectSample):
				currSubjectSample[taxaName] = [(weekSampleObtained, abundance)]
				weekOfLifeFirstSample = weekSampleObtained
			else:
				currSubjectSample[taxaName].append((weekSampleObtained, abundance))
				weekOfLifeLastSample = weekSampleObtained	
		previousSubjectID = subjectID
	#Close file
	infile.close()

	samplesPerSubject[previousSubjectID] = copy.copy(currSubjectSample)
	sampleSizePerSubject = len(currSubjectSample[metabolomeNames[0]])
	samplesPerSubjectInfo[previousSubjectID] = (weekOfLifeFirstSample, weekOfLifeLastSample, sampleSizePerSubject)

##	for subjectID, sampleInfo in samplesPerSubjectInfo.iteritems():
##		print subjectID, '\t', sampleInfo[0], '\t', sampleInfo[1], '\t', sampleInfo[2]

	taxonSplinesPerSubject = {}
	for subjectID in subjectIDs:
##		print "Processing subjectID -> ", subjectID
		abundanceByTaxa = samplesPerSubject[subjectID]
		splinesPerTaxa = {}
		weekFirstSample, weekLastSample, numSamples = samplesPerSubjectInfo[subjectID]
		if numSamples < MINIMUN_NUMBER_MEASURED_TIMEPOINTS:
##			print "\t", subjectID
			del samplesPerSubject[subjectID]
			continue
		#Get splines for each taxa across timepoints
		for taxaName, abundanceLevelPerTimepoint in abundanceByTaxa.iteritems():
			timepoints = []
			relativeAbundances = []
			for timepoint, abundance in abundanceLevelPerTimepoint:
				if not (timepoint in timepoints):
					timepoints.append(timepoint)
					relativeAbundances.append(abundance)
			mean = getMean(relativeAbundances)
			variance = getVariance(relativeAbundances)
			#Use B-spline to extrapolate values. NOTE: Parameters s must be adjusted appropriately to avoid over-fitting.
			tck = interpolate.splrep(timepoints, relativeAbundances, k=3, s=0.001, xb=weekFirstSample, xe=weekLastSample)
			splinesPerTaxa[taxaName] = copy.copy(tck)
##			#This code plots the original relative abudance data for a specific taxa along with a linear and two cubic B-splines approximations as long as the variance is above an arbitrary threshold. 
##			if variance > 0.0:
####				weights = [1/(math.sqrt(variance)) for i in xrange(len(timepoints))]
####				weights = [weight/sum(weights) for weight in weights]
##				t, c, k = interpolate.splrep(timepoints, relativeAbundances, k=3, s=0.001, xb=weekFirstSample, xe=weekLastSample)
##				sampleLength = weekLastSample - weekFirstSample + 1.0
####                                timepointsNew = np.arange(weekFirstSample, (weekLastSample + 1.0), 1.0)
##				timepointsNew = np.linspace(weekFirstSample, weekLastSample, num = sampleLength, endpoint = True)
##				relativeAbundancesSplev = interpolate.splev(timepointsNew, tck)
##				relativeAbundancesSplev = truncateAbundanceValues(relativeAbundancesSplev)
##				spline = interpolate.BSpline(t, c, k, extrapolate = False)
##				relativeAbundancesBspline = spline(timepointsNew)
##				relativeAbundancesBspline = truncateAbundanceValues(relativeAbundancesBspline)
##				fig = plt.figure() #plt.figure(figsize=(3, 6))
##				plt.plot(timepoints, relativeAbundances, 'x', timepointsNew, relativeAbundancesSplev, '-b', timepointsNew, relativeAbundancesBspline, '-g', timepoints, relativeAbundances, '--')
##				title = 'Relative abundance of ' + taxaName + ' for subject ' + subjectID
##				plt.legend(['Data', 'Splev', 'BSpline', 'Linear'])
##				plt.title(title)
##				plt.show()
####				fig.savefig(taxaName + '_' + subjectID + '.png', dpi=fig.dpi)
####				plt.close()
		taxonSplinesPerSubject[subjectID] = copy.copy(splinesPerTaxa)

##	print len(samplesPerSubject)

	taxonSamplesPerSubject = {}
	for subjectID in samplesPerSubject.keys():
		sample = samplesPerSubject[subjectID]
		abundanceLevelPerTimepoint = sample[metabolomeNames[0]]
		headers = ['TaxaName']
		for timepoint, abundance in abundanceLevelPerTimepoint:
			headers.append(timepoint)
		sampleTaxonAbundances = [headers]
		for taxaName, abundanceLevelPerTimepoint in sample.iteritems():
			currTaxaAbundances = [taxaName]
			for timepoint, abundance in abundanceLevelPerTimepoint:
				currTaxaAbundances.append(abundance)
			sampleTaxonAbundances.append(currTaxaAbundances)
		taxonSamplesPerSubject[subjectID] = sampleTaxonAbundances

##	print taxonSamplesPerSubject['M4']

	return taxonSamplesPerSubject, taxonSplinesPerSubject, samplesPerSubjectInfo

def getAlignmentsGivenReference(referenceSampleID, taxonSamples, splinesPerSubject, samplesPerSubjectInfo, useSplines, outfilename):
	outfile = open(outfilename, 'a')
##	outline = 'Reference SampleID' + '\t' + 'Aligned SampleID' + '\t' + 'TaxonError' + '\t' + 'a' + '\t' + 'b' + '\t' + 'alpha' + '\t' + 'beta' + '\t' + 'overlap' + '\t' + 'Taxa Names' + '\t' + 'Alignment Scores (Pearson, Spearman, SSD, R-squared)' + '\n'
	outline = 'Reference SampleID' + '\t' + 'Aligned SampleID' + '\t' + 'TaxonError' + '\t' + 'a' + '\t' + 'b' + '\t' + 'alpha' + '\t' + 'beta' + '\t' + 'overlap' + '\t' + 'Taxa Names' + '\t' + 'Alignment Scores (SSD)' + '\n'
	outfile.writelines(outline)
	subjectIDs = splinesPerSubject.keys()
	referenceSampleweekFirstSample, referenceSampleweekLastSample, referenceSampleSize = samplesPerSubjectInfo[referenceSampleID]
	for i in xrange(0, len(subjectIDs)):
		currentSampleID = copy.copy(subjectIDs[i])
		currentSampleweekFirstSample, currentSampleweekLastSample, currentSampleSize = samplesPerSubjectInfo[currentSampleID]
		if currentSampleID == referenceSampleID:
			continue
		#Get taxon info for reference sample
		taxonAbundancesReferenceSample = taxonSamples[referenceSampleID]
		taxonSplinesReferenceSample = splinesPerSubject[referenceSampleID]
		#Get taxon info for candidate sample
		taxonAbundancesCurrentSample = taxonSamples[currentSampleID]
		taxonSplinesCurrentSample = splinesPerSubject[currentSampleID]
		
		outline = referenceSampleID + '\t' + currentSampleID
		print 'Processing current alignment between samples', referenceSampleID, 'and', currentSampleID 

		taxonAgreementScorePostAlignment, taxonError, a, b, alpha, beta, overlap = getPairwiseAlignment(referenceSampleID, taxonAbundancesReferenceSample, taxonSplinesReferenceSample, currentSampleID, taxonAbundancesCurrentSample, taxonSplinesCurrentSample, samplesPerSubjectInfo, useSplines)
		if len(taxonAgreementScorePostAlignment) < 1:
			outline += '\n'
			outfile.writelines(outline)
			continue
		outline += '\t' + str(taxonError) + '\t' + str(a) + '\t' + str(b) + '\t' + str(alpha) + '\t' + str(beta) + '\t' + str(overlap)
		for taxaName, taxaScores in taxonAgreementScorePostAlignment.iteritems():
			lineScores = []
			for taxaScore in taxaScores:
				lineScores.append(str(taxaScore))                                     
			outline += '\t' + taxaName + '\t' + ','.join(lineScores)
		outline += '\n'
		outfile.writelines(outline)
	#Close output file
	outfile.close()

	return

def getAllPairwiseAlignments(taxonSamples, splinesPerSubject, samplesPerSubjectInfo, useSplines, outfilename):
	outfile = open(outfilename, 'a')
##	outline = 'Reference SampleID' + '\t' + 'Aligned SampleID' + '\t' + 'TaxonError' + '\t' + 'a' + '\t' + 'b' + '\t' + 'alpha' + '\t' + 'beta' + '\t' + 'overlap' + '\t' + 'Taxa Names' + '\t' + 'Alignment Scores(Pearson, Spearman, SSD, R-squared)' + '\n'
	outline = 'Reference SampleID' + '\t' + 'Aligned SampleID' + '\t' + 'TaxonError' + '\t' + 'a' + '\t' + 'b' + '\t' + 'alpha' + '\t' + 'beta' + '\t' + 'overlap' + '\t' + 'Taxa Names' + '\t' + 'Alignment Scores (SSD)' + '\n'
	outfile.writelines(outline)
	subjectIDs = splinesPerSubject.keys()
#	referenceSampleID = '48'
	for i in xrange(0, len(subjectIDs) - 1):
		sample1SubjectID = copy.copy(subjectIDs[i])
		sample1weekFirstSample, sample1weekLastSample, sample1Size = samplesPerSubjectInfo[sample1SubjectID]
		for j in xrange(i + 1, len(subjectIDs)):
			sample2SubjectID = copy.copy(subjectIDs[j])
			sample2weekFirstSample, sample2weekLastSample, sample2Size = samplesPerSubjectInfo[sample2SubjectID]
			if sample1Size >= sample2Size:
				referenceSampleID = copy.copy(sample1SubjectID)
				currentSampleID = copy.copy(sample2SubjectID)
			else:
				referenceSampleID = copy.copy(sample2SubjectID)
				currentSampleID = copy.copy(sample1SubjectID)
			if currentSampleID == referenceSampleID:
				continue
			#Get taxon info for reference sample
			taxonAbundancesReferenceSample = taxonSamples[referenceSampleID]
			taxonSplinesReferenceSample = splinesPerSubject[referenceSampleID]
			#Get taxon info for candidate sample
			taxonAbundancesCurrentSample = taxonSamples[currentSampleID]
			taxonSplinesCurrentSample = splinesPerSubject[currentSampleID]
			
			outline = referenceSampleID + '\t' + currentSampleID
			print 'Processing current alignment between samples', referenceSampleID, 'and', currentSampleID

			taxonAgreementScorePostAlignment, taxonError, a, b, alpha, beta, overlap = getPairwiseAlignment(referenceSampleID, taxonAbundancesReferenceSample, taxonSplinesReferenceSample, currentSampleID, taxonAbundancesCurrentSample, taxonSplinesCurrentSample, samplesPerSubjectInfo, useSplines)
			if len(taxonAgreementScorePostAlignment) < 1:
				outline += '\n'
				outfile.writelines(outline)
				continue
			outline += '\t' + str(taxonError) + '\t' + str(a) + '\t' + str(b) + '\t' + str(alpha) + '\t' + str(beta) + '\t' + str(overlap)
			for taxaName, taxaScores in taxonAgreementScorePostAlignment.iteritems():
				lineScores = []
				for taxaScore in taxaScores:
					lineScores.append(str(taxaScore))                                     
				outline += '\t' + taxaName + '\t' + ','.join(lineScores)
			outline += '\n'
			outfile.writelines(outline)
	#Close output file
	outfile.close()

	return

def getPairwiseAlignment(referenceSampleID, taxonAbundancesReferenceSample, taxonSplinesReferenceSample, currentSampleID, taxonAbundancesCurrentSample, taxonSplinesCurrentSample, sampleInfo, useSplines):
	taxonWorkingSet = {}
	for currRow in taxonAbundancesReferenceSample:
		taxonWorkingSet[currRow[0]] = copy.copy(currRow[0])

	taxonReferenceSample = buildTaxon(taxonAbundancesReferenceSample, taxonSplinesReferenceSample)
	taxonCurrentSample = buildTaxon(taxonAbundancesCurrentSample, taxonSplinesCurrentSample)

	[filteredTaxonReferenceSample, filteredTaxonCurrentSample] = filterTaxon(taxonReferenceSample, taxonCurrentSample, useSplines)
	if len(filteredTaxonReferenceSample) < 2 or len(filteredTaxonCurrentSample) < 2:
##		print "\tCurrent alignment skipped due to lack of shared taxa after filtering ... "
		return [], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0

	timepointsListReferenceSample = buildTimepointsProfile(filteredTaxonReferenceSample, sampleInfo[referenceSampleID], useSplines)
	timepointsListCurrentSample = buildTimepointsProfile(filteredTaxonCurrentSample, sampleInfo[currentSampleID], useSplines)

##	return [], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0

	taxonWeights = [1.0 for taxa in xrange(len(timepointsListReferenceSample[0].relativeAbundance))]
	taxonWeights = [taxaWeight / sum(taxonWeights) for taxaWeight in taxonWeights]

	#Get optimal alignment between reference sample and current sample
	optimalAlignmentInfo = getOptimalMapping(timepointsListReferenceSample, timepointsListCurrentSample, taxonWeights, useSplines)
	if len(optimalAlignmentInfo) < 1:
##		print "\tLinear warp method failed to find an alignment with at least", (OVERLAP_THRESHOLD * 100), "% overlap between samples. Consider increasing the range of a or b."
		return [], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
##	print optimalAlignmentInfo[0], optimalAlignmentInfo[1], optimalAlignmentInfo[2], optimalAlignmentInfo[3], optimalAlignmentInfo[4], optimalAlignmentInfo[5]
	#For each taxa, compute the agreement with optimal alignment
	taxonAgreementScorePostAlignment = getAlignmentAgreementScorePerTaxa(timepointsListReferenceSample, timepointsListCurrentSample, optimalAlignmentInfo[1], optimalAlignmentInfo[2], optimalAlignmentInfo[3], optimalAlignmentInfo[4], taxonWeights, method='rsquared')

	filteredTaxon = timepointsListReferenceSample[0].relativeAbundance.keys()
	alignmentWeigthedAgreementPerTaxa = [[taxonAgreementScorePostAlignment[filteredTaxon[taxaIndex]], filteredTaxon[taxaIndex]] for taxaIndex in xrange(len(taxonAgreementScorePostAlignment))] #Ignoring weights for now
	alignmentWeigthedAgreementPerTaxa.sort(reverse=True)

##	#Plot samples pre- and post-alignment for each taxa in decreasing order of agreement score 
##	for alignmentWeigthedAgreementScore, taxaName in alignmentWeigthedAgreementPerTaxa:
##		plotAlignment(timepointsListReferenceSample, timepointsListCurrentSample, optimalAlignmentInfo[1], optimalAlignmentInfo[2], optimalAlignmentInfo[3], optimalAlignmentInfo[4], taxaName, referenceSampleID, currentSampleID, taxonAgreementScorePostAlignment[taxaName])
			
	return taxonAgreementScorePostAlignment, optimalAlignmentInfo[0], optimalAlignmentInfo[1], optimalAlignmentInfo[2], optimalAlignmentInfo[3], optimalAlignmentInfo[4], optimalAlignmentInfo[5]

def getMean(sampleValues):
	meanValue = sum(sampleValues) / float(len(sampleValues))

	return meanValue

def getVariance(sampleValues):
	variance = 0.0
	meanValue = sum(sampleValues) / float(len(sampleValues))
	for currentValue in sampleValues:
		variance += (currentValue - meanValue)**2
	variance = variance / float(len(sampleValues))

	return variance

def getSSD(referenceSampleValues, currentSampleValues):
	ssd = 0.0
	if len(referenceSampleValues) != len(currentSampleValues):
		print "<<ERROR>> Temporal data is not of the same length. Consider re-sampling one time series to the length of the other and try again."
		sys.exit()
	for i in xrange(len(referenceSampleValues)):
		ssd += (referenceSampleValues[i] - currentSampleValues[i])**2

	return ssd

def truncateAbundanceValues(sampleAbundanceValues):
	sampleAbundanceValues = np.asarray(sampleAbundanceValues)
	lowValues = sampleAbundanceValues < LOWER_BOUND
	highValues = sampleAbundanceValues > UPPER_BOUND
	sampleAbundanceValues[lowValues] = LOWER_BOUND
	sampleAbundanceValues[highValues] = UPPER_BOUND

	return sampleAbundanceValues

def main(argv):
	if (len(argv) == 4):
		dataFilename = argv[1]
		useSplines = bool(argv[2])
		outfilename = argv[3]
	else:
		print "<<ERROR>> Invalid number of parameters!"
		return

	#Read dataset and prepare corresponding data structures
	taxonSamplesPerSubject, taxonSplinesPerSubject, samplesPerSubjectInfo = getSamples(dataFilename)

	referenceSampleID = '102'
	print "Processing pairwise alignments across all subjects against reference sample", referenceSampleID, "..."
#	#Get all pairwise alignments agaist a given reference sample.
	getAlignmentsGivenReference(referenceSampleID, taxonSamplesPerSubject, taxonSplinesPerSubject, samplesPerSubjectInfo, useSplines, outfilename)

#	print "Processing all pairwise alignments between subjects ..."
#	#Get all pairwise alignments between subjects.
#	getAllPairwiseAlignments(taxonSamplesPerSubject, taxonSplinesPerSubject, samplesPerSubjectInfo, useSplines, outfilename)

if __name__ == '__main__':
##	main(sys.argv)
	#Alternate call, if you are using IDLE. 
	main(['getAlignmentsIBD.py', 'HMDBMetabolitesNormalizedRows.csv', 'True', 'human_ibd_microbiota_metabolites_pairwise_subject_ranking.tsv'])
#	main(['getAlignmentsIBD.py', 'HMDBMetabolitesNormalizedRows.csv', 'True', 'human_ibd_microbiota_metabolites_pairwise_all_subject_ranking.tsv'])

