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
import time
from random import randint
import matplotlib.pyplot as plt

dataSet = None
ndim = None
power = 0.8
#power = 1
pq = None
done = False
#avgDeltaSSQ = None
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
			for i in xrange(ndim):
				self.limitsLow[i] = min(first.limitsLow[i], second.limitsLow[i])
				self.limitsHigh[i] = max(first.limitsHigh[i], second.limitsHigh[i])
			self.keys = self.findKeys()
			self.weight = first.weight + second.weight
			self.Sum = np.add(first.Sum, second.Sum)
			self.sqSum = np.add(first.sqSum, second.sqSum)#self.computeSqSum()
			self.SSQ = self.computeSSQ()
			self.CoG = [None]*ndim
		else:
			# The first two parameters are the edges of the cluster
			self.limitsLow = limitsLow
			self.limitsHigh = limitsHigh
			self.keys = self.findKeys()
			self.weight = computeWeightIn(self.limitsLow, self.limitsHigh)
			self.Sum = computeSum(self.limitsLow, self.limitsHigh)
			self.sqSum = self.computeSqSum()
			self.SSQ = self.computeSSQ()
			self.CoG = [None]*ndim
	def computeSSQ(self):
		# Compute SSQ of the cluster
		self.SSQ = 0
		self.computeSqSum()
		if self.weight is not 0:
			for i in xrange(ndim):
				self.SSQ += ( self.sqSum[i] - pow(self.Sum[i],2)/self.weight )
		return self.SSQ
	def computeAvgDeltaSSQ(self):
		#self.avgDeltaSSQ = self.SSQ / self.weight

		self.avgDeltaSSQ = self.SSQ / len(self.keys)

		#keys = keylist[:]
		#for i in xrange(ndim):
		#	ckeys = [key for key in keys if key[i] >= self.limitsLow[i] and key[i] <= self.limitsHigh[i]]
		#	keys = ckeys
		#n = [0]*ndim
		#for i in xrange(ndim):
		#	n[i] = len(set([key[i] for key in keys]))
		#	self.avgDeltaSSQ += self.SSQ / n[i]
		##self.avgDeltaSSQ = sum(avgDeltaSSQperDim)
		return self.avgDeltaSSQ
	def computeCoG(self):
		# Find center of gravity of the cluster 
		self.CoG = tuple(x/self.weight for x in self.Sum)
		return self.CoG
	def findKeys(self):
		keys = keylist[:]
		for i in xrange(ndim):
			ckeys = [key for key in keys if key[i] >= self.limitsLow[i] and key[i] <= self.limitsHigh[i]]
			keys = ckeys
		self.keys = keys[:]
		return self.keys
	def computeSqSum(self):
		# For each dimension, compute sum of square coordinates
		self.sqSum = [0]*ndim
		keys = keylist[:]
		self.sqSum = list(self.sqSum)
		for i in xrange(ndim):
			ckeys = [key for key in keys if key[i] >= self.limitsLow[i] and key[i] <= self.limitsHigh[i]]
			keys = ckeys
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
	def setChildren(self, children):
		left,right = children
		# Sets two clusters as children of the current node
		self.left = newTreeNode(left.limitsLow, left.limitsHigh)
		self.right = newTreeNode(right.limitsLow, right.limitsHigh)

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
		#print self.listOfClusters
	def addCluster(self, cl):
		self.listOfClusters.append(cl)
	def deleteCluster(self, cl):
		self.listOfClusters.remove(cl)

def computeWeightIn(limitsLow, limitsHigh):
	""" Compute content of cluster """
	keys = keylist[:]
	weight = 0
	for i in xrange(ndim):
		ckeys = [key for key in keys if key[i] >= limitsLow[i] and key[i] <= limitsHigh[i]]
		keys = ckeys
	for key in keys:
		#print key, dataSet[key]
		weight += dataSet[key]
	return weight

def computeSum(limitsLow, limitsHigh):
	""" Compute vector sum of the cluster """
	keys = keylist[:]
	vecSum = [0]*ndim
	for i in xrange(ndim):
		ckeys = [key for key in keys if key[i] >= limitsLow[i] and key[i] <= limitsHigh[i]]
		keys = ckeys
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
		if testOverlapping(c1, c2, i) or testOverlapping(c2, c1, i):
			dimCount += 1
	if dimCount >= (ndim-1):
		return True
	# In any other case:
	return False

def testOverlapping(c1, c2, i):
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
	bestPair = listOfMergeablePairs[0]
	return bestPair[0], bestPair[1]

def splitCluster(cl):
	""" Returns the margins of the two children or False """
	avgDeltaSSQ = cl.computeAvgDeltaSSQ()
	refDssq = 0
	leftLimitsLow = [None]*ndim
	leftLimitsHigh = [None]*ndim
	rightLimitsHigh = [None]*ndim
	rightLimitsLow = [None]*ndim
	subClusters = [None]*2
	# Select only the range of the current cluster
	keys = keylist[:]
	for i in xrange(ndim):
		ckeys = [key for key in keys if key[i] > cl.limitsLow[i] and key[i] < cl.limitsHigh[i]]
		keys = ckeys
	for i in xrange(ndim):
		for key in keys:
			leftLimitsLow = cl.limitsLow[:]
			leftLimitsHigh = cl.limitsHigh[:]
			leftLimitsHigh[i] = key[i]
			rightLimitsLow = cl.limitsLow[:]
			rightLimitsHigh = cl.limitsHigh[:]
			rightLimitsLow[i] = key[i]
			# Create the two subclusters from the mother
			left = newCluster(leftLimitsLow, leftLimitsHigh)
			right = newCluster(rightLimitsLow, rightLimitsHigh)
			# Compute weighted Delta SSQ
			#print key
			if left.weight>0 and right.weight>0:
				newDssq = pow(computeDeltaSSQ(left, right), power)
				#print computeDeltaSSQ(left, right), newDssq
				if (newDssq > refDssq) and (newDssq > avgDeltaSSQ):
					print computeDeltaSSQ(left,right), newDssq, avgDeltaSSQ
					# Look for the maximum weighted Delta SSQ
					refDssq = newDssq
					subClusters = [left, right]
					return subClusters
	return False

def areDifferentClusters(c1, c2):
	""" Check if two newCluster() objects refer to the same physical cluster """
	w1, w2 = c1.weight, c2.weight
	s1, s2 = c1.Sum, c2.Sum
	if (w2 is not w1) and (s2 is not s1):
		return True
	else:
		return False

def mergeClusters(c1, c2):
	""" Merge two adjacent clusters """
	return newCluster(c1, c2, "merging")

def generateListOfMergeablePairs():
	global listOfMergeablePairs
	listOfMergeablePairs = []
	for c1 in pq.listOfClusters:
		for c2 in pq.listOfClusters:
			if areDifferentClusters(c1, c2) and areAdjacent(c1, c2):
				dssq = computeDeltaSSQ(c1, c2)
				listOfMergeablePairs.append( (c1, c2, dssq) )

def initStructures():
	""" Initialize priority queue and binary tree """
	#global avgDeltaSSQ, pq, keylist
	global pq, keylist
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
	#avgDeltaSSQ = root.clusterData.computeSSQ() / np.sqrt(len(keylist))
	return True

def topDownSplitting():
	""" Split the data set into micro-clusters """
	global pq, done
	while not done:
		currentNode = pq.get()
		currentCluster = currentNode[0].clusterData
		children = splitCluster(currentCluster)
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
	global listOfMergeablePairs, pq, done
	done = False
	pq.makeListOfClusters()
	# Generate the list of pairs of clusters that can be merged
	generateListOfMergeablePairs()
	if len(listOfMergeablePairs) is 0:
		done = True
	if len(pq.listOfClusters) < 2:
		done = True
	while not done:
		c1, c2 = findClustersToMerge() # a pair of two clusters
		minSSQincrease = computeDeltaSSQ(c1, c2)
		largerCluster = mergeClusters(c1, c2)
		avgDeltaSSQ = np.average(largerCluster.computeAvgDeltaSSQ())
		# Confirm the merging only if minSSQincrease < avgDeltaSSQ
		if minSSQincrease < avgDeltaSSQ:
			print minSSQincrease, avgDeltaSSQ, c1.computeAvgDeltaSSQ(), c2.computeAvgDeltaSSQ()
			# Remove the two merged clusters from the list of current custers...
			pq.deleteCluster(c1)
			pq.deleteCluster(c2)
			# ...and add the newly created one
			pq.addCluster(largerCluster)
			# Regenerate the list of pairs of clusters that can be merged
			if len(pq.listOfClusters) < 2:
				done = True
			else:
				generateListOfMergeablePairs()
				if len(listOfMergeablePairs) is 0:
					done = True
				#else:
				#	c1, c2 = findClustersToMerge()
				#	minSSQincrease = computeDeltaSSQ(c1, c2)
				#	largerCluster = mergeClusters(c1, c2)
				#	avgDeltaSSQ = np.average(largerCluster.computeAvgDeltaSSQ())
		else:
			done = True
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
		listOfCoGs.append((cl.computeCoG(), cl.weight))
	return listOfCoGs

def CLUBSclustering(data, nd):
	""" Execute findClusters() and handle I/O """
	#global dataSet, ndim, power, pq, done, avgDeltaSSQ, listOfMergeablePairs
	global dataSet, ndim, power, pq, done, listOfMergeablePairs
	dataSet = data
	ndim = nd
	#power = 0.8
	pq = newPriorityQueue()
	done = False
	#avgDeltaSSQ = 0
	listOfMergeablePairs = []
	return findClusters()

	
if __name__ == '__main__':
	"""sample usage"""
	dim = 50
	DS = dict([((x,y),0) for x in range(dim) for y in range(dim)])	
	#for i in range(100):
	#    DS[(randint(30,dim-1),randint(0,25))] = 1
	#for i in range(60):
	#    DS[(randint(30,dim-1),randint(29,41))] = 1
	#for i in range(10):
	#    DS[(randint(30,dim-1),randint(43,dim-1))] = 1
	#for i in range(50):
	#    DS[(randint(0,10),randint(0,10))] = 1
	#for i in range(10):
	#    DS[(randint(0,10),randint(14,16))] = 1
	#for i in range(50):
	#    DS[(randint(0,10),randint(19,26))] = 1
	#for i in range(100):
	#    DS[(randint(0,10),randint(29,46))] = 1
	#for i in range(10):
	#    DS[(randint(0,10),randint(48,dim-1))] = 1

	#for i in range(10):
	#    DS[(randint(30,dim-1),randint(0,dim-1))] = 1
	#for i in range(100):
	#    DS[(randint(0,20),randint(0,dim-1))] = 1
	borders = 4, 45

	for i in range(10):
		center = [0]*2
		center[0] = randint(borders[0], borders[1])
		center[1] = randint(borders[0], borders[1])
		dim = randint(3,15)
		for j in xrange(dim):
			DS[(randint(center[0]-3,center[0]+3),randint(center[1]-3,center[1]+3))] = randint(3,7)


	x,y,w = [],[],[]
	for item in DS.keys():
		if DS[item] > 0:                  
			x.append(item[0])
			y.append(item[1])
			w.append(DS[item])
	#plt.ion()
	plt.scatter(x,y,c=w,s=200)
	plt.gray()
	plt.grid()
	plt.draw()
	#plt.pause(0.001)
	clusters = CLUBSclustering(DS,2)
	print clusters
	xc, yc, wc = [],[],[]
	for i in xrange(len(clusters)):
		xc.append(clusters[i][0][0])
		yc.append(clusters[i][0][1])
		wc.append(clusters[i][1])
		print xc[i], yc[i], wc[i]
	area = [200*np.sqrt(w) for w in wc]
	plt.scatter(xc,yc,c='red',s=area,alpha=0.5)
	plt.show()


