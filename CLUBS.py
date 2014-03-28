#!/usr/bin/env python
#
# Author: Elena Graverini
# Contact: elena.graverini@cern.ch
#
# Implementation of the M-CLUBS algorithm
# by E. Masciari, G.M. Mazzeo, C. Zaniolo
#
# Original documentation in:
# "Analysing microarray expression data
# through effective clustering",
# Information Sciences, Volume 262 (2014),
# pages 32-45
#
# The program takes a dictionary as input data set:
#     key = tuple of coordinates (binned)
#     value = bin content
# It divides the data set into clusters and 
# returns their centers of gravity.

import numpy as np
import sys
from random import randint

dataSet = None
ndim = None
power = 0.8
pq = None
done = False
avgDeltaSSQ = None
listOfMergeablePairs = None
keylist = None

class newCluster():
	""" Class holding cluster data (edges and statistics) """
	def __init__(self, limitsLow, limitsHigh, status = None):
		# Either merge two clusters (bottom-up)
		# or create a new micro-cluster (top-down)
		if status is "merging":
			# The first two input parameters are two newCluster() objects
			first = limitsLow
			second = limitsHigh
			self.limitsLow = [None]*ndim
			self.limitsHigh = [None]*ndim
			self.weight = first.weight + second.weight
			self.Sum = first.Sum + second.Sum
			for i in xrange(ndim):
				self.limitsLow[i] = min(first.limitsLow[i], second.limitsLow[i])
				self.limitsHigh[i] = max(first.limitsHigh[i], second.limitsHigh[i])
			self.sqSum = [0]*ndim
			self.CoG = [None]*ndim
			self.SSQ = 0
		else:
			# The first two parameters are the edges of the cluster
			self.limitsLow = limitsLow
			self.limitsHigh = limitsHigh
			self.weight = computeWeightIn(self.limitsLow, self.limitsHigh)
			self.Sum = computeSum(self.limitsLow, self.limitsHigh)
			self.sqSum = [0]*ndim
			self.CoG = [None]*ndim
			self.SSQ = 0
	def computeSSQ(self):
		# Compute SSQ of the cluster
		self.computeSqSum()
		for i in xrange(ndim):
			self.SSQ += ( self.sqSum[i] - pow(self.Sum[i],2)/self.weight )
		return self.SSQ
	def computeCoG(self):
		# Find center of gravity of the cluster 
		self.CoG = tuple(x/self.weight for x in self.Sum)
		return self.CoG
	def computeSqSum(self):
		# For each dimension, compute sum of square coordinates
		keys = keylist
		self.sqSum = list(self.sqSum)
		for i in xrange(ndim):
			keys = [key for key in keys if key[i] >= self.limitsLow[i] and key[i] <= self.limitsHigh[i]]
		for i in xrange(ndim):
			for key in keys:
				# A bin with weight W counts as W times that bin
				self.sqSum[i] += ( pow(key[i],2) * dataSet[key] )
		# Make the vector sum of squares an immutable object
		self.sqSum = tuple(self.sqSum)
		return self.sqSum

class newTreeNode():
	""" Auxiliary binary tree for top-down splitting """
	def __init__(self, limitsLow, limitsHigh):
		self.left = None
		self.right = None
		self.clusterData = newCluster(limitsLow, limitsHigh)
	def setChildren(self, left, right):
		# Sets two clusters as children of the current node
		self.left = newTreeNode(left.limitsLow, left.limitsHigh)
		self.right = newTreeNode(left.limitsLow, left.limitsHigh)

class newPriorityQueue():
	""" Selects the cluster with the highest SSQ """
	def __init__(self):
		self.queue = []
		self.listOfClusters = []
	def add(self, node):
		self.queue.append((node, node.clusterData.computeSSQ()))
	def get(self):
		# Returns the cluster with the highest SSQ
		self.queue = sorted(self.queue, key = lambda x: x[1])
		return self.queue[-1]
	def delete(self, cl):
		self.queue.remove(cl)
	def makeListOfClusters(self):
		# For the bottom-up part of the algorithm
		self.listOfClusters = [item[0].clusterData for item in self.queue]
		print self.listOfClusters
	def addCluster(self, cl):
		self.listOfClusters.append(cl)
	def deleteCluster(self, cl):
		self.listOfClusters.remove(cl)

def computeWeightIn(limitsLow, limitsHigh):
	""" Compute content of cluster """
	keys = keylist
	weight = 0
	for i in xrange(ndim):
		keys = [key for key in keys if key[i] >= limitsLow[i] and key[i] <= limitsHigh[i]]
	for key in keys:
		#print key, dataSet[key]
		weight += dataSet[key]
	return weight

def computeSum(limitsLow, limitsHigh):
	""" Compute vector sum of the cluster """
	keys = keylist
	vecSum = [0]*ndim
	for i in xrange(ndim):
		keys = [key for key in keys if key[i] >= limitsLow[i] and key[i] <= limitsHigh[i]]
	for i in xrange(ndim):
		for key in keys:
			# A bin with weight W counts as W times that bin
			vecSum[i] += (key[i] * dataSet[key])
	# Make the vector sum an immutable object
	return tuple(vecSum)

def computeDeltaSSQ(c1, c2):
	""" Variation in SSQ due to the splitting of a cluster or to the merging of two sub-clusters"""
	dssq = 0
	#print c1.weight, c2.weight
	factor = (c1.weight * c2.weight) / (c1.weight + c2.weight)
	for i in xrange(ndim):
		dssq += pow( c1.Sum[i]/c1.weight - c2.Sum[i]/c2.weight, 2 )
	return factor * dssq

def areAdjacent(c1, c2):
	""" Find out if two clusters are adjacent or overlapping """
	# Case 1: one edge in common
	for i in xrange(ndim):
		if (c1.limitsHigh[i] == c2.limitsLow[i]) or (c2.limitsHigh[i] == c1.limitsLow[i]):
			return True
	# Case 2: partially overlapping (in at least ndim-1 dimensions)
	dimCount = 0
	for i in xrange(ndim):
		if testOverlapping(c1, c2) or testOverlapping(c2, c1):
			dimCount += 1
	if dimCount >= (ndim-1):
		return True
	# In any other case:
	return False

def testOverlapping(c1, c2):
	""" Check if two clusters overlap in one dimension """
	c1.computeCoG()
	c2.computeCoG()
	spatialOverlapping = bool(c2.limitsLow[i] < c1.limitsHigh[i])
	proximity = bool(dist(c1.CoG, c2.CoG) <= (np.abs(c1.CoG[i]-c1.limitsHigh[i]) + np.abs(c2.CoG[i]-c2.limitsLow[i])))
	if (spatialOverlapping and proximity):
		return True
	else:
		return False

def dist(p1, p2):
	""" Compute the distance between two points """
	dist = 0
	for i in xrange(ndim):
		dist += pow(p1[i]-p2[i], 2)
	return np.sqrt(dist)

def findClustersToMerge():
	""" Find the pair of adjacent clusters that yiealds the least SSQ increase """
	global listOfMergeablePairs
	listOfMergeablePairs = sorted(listOfMergeablePairs, key = lambda x: x[2])
	if len(listOfMergeablePairs) is 0:
		return False
	bestPair = listOfMergeablePairs[0]
	return bestPair[0], bestPair[1]

def splitCluster(cl, avgDeltaSSQ):
	""" Returns the margins of the two children or False """
	keys = keylist
	refDssq = 0
	leftLimitsLow = [None]*ndim
	leftLimitsHigh = [None]*ndim
	rightLimitsHigh = [None]*ndim
	rightLimitsLow = [None]*ndim
	subClusters = [None]*2
	# Select only the range of the current cluster
	for i in xrange(ndim):
		keys = [key for key in keys if key[i] > cl.limitsLow[i] and key[i] < cl.limitsHigh[i]]
	for i in xrange(ndim):
		for key in keys:
			leftLimitsLow = cl.limitsLow
			leftLimitsHigh = cl.limitsHigh
			leftLimitsHigh[i] = key[i]
			rightLimitsLow = cl.limitsLow
			rightLimitsHigh = cl.limitsHigh
			rightLimitsLow[i] = key[i]
			# Create the two subclusters from the mother
			left = newCluster(leftLimitsLow, leftLimitsHigh)
			right = newCluster(rightLimitsLow, rightLimitsHigh)
			# Compute weighted Delta SSQ
			#print key
			if left.weight>0 and right.weight>0:
				newDssq = pow(computeDeltaSSQ(left, right), power)
				if newDssq > refDssq:
					# Look for the maximum weighted Delta SSQ
					refDssq = newDssq
					subClusters = [left, right]
	if refDssq > avgDeltaSSQ:
		return subClusters
	else:
		return False

def mergeClusters(c1, c2):
	""" Merge two adjacent clusters """
	return newCluster(c1, c2, "merging")

def initStructures():
	""" Initialize priority queue and binary tree """
	global avgDeltaSSQ, pq, keylist
	keylist = dataSet.keys()
	# Edges of the data set
	marginsLow = [None]*ndim
	marginsHigh = [None]*ndim
	for i in xrange(ndim):
		sortedkeys = sorted(keylist, key = lambda x: x[i])
		marginsLow[i] = sortedkeys[0][i]
		marginsHigh[i] = sortedkeys[-1][i]
	# Initialize the root of the binary tree as the full data set
	root = newTreeNode(marginsLow, marginsHigh)
	pq.add(root)
	#avgDeltaSSQ = root.clusterData.computeSSQ() / root.clusterData.weight
	avgDeltaSSQ = root.clusterData.computeSSQ() / np.sqrt(len(keylist))
	return True

def topDownSplitting():
	""" Split the data set into micro-clusters """
	global pq, done
	while not done:
		currentNode = pq.get()
		currentCluster = currentNode[0].clusterData
		children = splitCluster(currentCluster, avgDeltaSSQ)
		if not children:
			done = True
		else:
			currentNode[0].setChildren(children)
			pq.delete(currentNode)
			pq.add(currentNode[0].left)
			pq.add(currentNode[0].right)
	print "Data set split into %s micro-clusters."%len(pq.queue)
	return True

def bottomUpClustering():
	""" Merge micro-clusters to provide the final clustering """
	global listOfMergeablePairs, pq
	pq.makeListOfClusters()
	for c1 in pq.listOfClusters:
		for c2 in pq.listOfClusters:
			if (c2 is not c1) and areAdjacent(c1, c2):
				dssq = computeDeltaSSQ(c1, c2)
				listOfMergeablePairs.append( (c1, c2, dssq) )
	clustersToMerge = findClustersToMerge() # a pair of two clusters
	if not clustersToMerge:
		return True
	minSSQincrease = computeDeltaSSQ(clustersToMerge)
	while minSSQincrease < avgDeltaSSQ:
		largerCluster = mergeClusters(clustersToMerge)
		listOfMergeablePairs = [item for item in listOfMergeablePairs if item[0] is not clustersToMerge[0] and item[0] is not clustersToMerge[1] and item[1] is not clustersToMerge[0] and item[1] is not clustersToMerge[1]]
		pq.deleteCluster(clustersToMerge[0])
		pq.deleteCluster(clustersToMerge[1])
		pq.addCluster(largerCluster)
		for cl in pq.listOfClusters:
			if (cl is not largerCluster) and areAdjacent(largerCluster, cl):
				dssq = computeDeltaSSQ(largerCluster, cl)
				listOfMergeablePairs.append( (largerCluster, cl, dssq) )
		clustersToMerge = findClustersToMerge()
		minSSQincrease = computeDeltaSSQ(clustersToMerge)
	print "Micro-clusters regrouped into %s clusters."%len(pq.listOfClusters)
	return True

def findClusters():
	""" Main program: returns a list of clusters in pq """
	if not initStructures():
		print "Error: initStructures"
		sys.exit(1)
	#sys.exit(0)
	if not topDownSplitting():
		print "Error: topDownSplitting"
		sys.exit(1)
	if not bottomUpClustering():
		print "Error: bottomUpClustering"
		sys.exit(1)
	# Return the centers of gravity of found clusters
	listOfCoGs = []
	for cl in pq.listOfClusters:
		listOfCoGs.append(cl.computeCoG())
	return listOfCoGs

def CLUBSclustering(data, nd):
	""" Execute findClusters() and handle I/O """
	global dataSet, ndim, power, pq, done, avgDeltaSSQ, listOfMergeablePairs
	dataSet = data
	ndim = nd
	power = 0.8
	pq = newPriorityQueue()
	done = False
	avgDeltaSSQ = 0
	listOfMergeablePairs = []
	return findClusters()

	
if __name__ == '__main__':
	"""sample usage"""
	dim = 50
	DS = dict([((x,y),0) for x in range(dim) for y in range(dim)])	
	for i in range(1800):
	    DS[(randint(30,dim-1),randint(0,dim-1))] = randint(0,20)
	for i in range(700):
	    DS[(randint(0,20),randint(0,dim-1))] = randint(0,20)
	clusters = CLUBSclustering(DS,2)
	print clusters


