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
import matplotlib.patches
import pylab

dataSet = None
ndim = None
power = 0.8
pq = None
done = False
avgDeltaSSQ = 0
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
			self.findKeys()
			self.weight = first.weight + second.weight
			self.Sum = np.add(first.Sum, second.Sum)
			self.sqSum = np.add(first.sqSum, second.sqSum)#self.computeSqSum()
			self.computeSSQ()
			self.computeCoG()
		else:
			# The first two parameters are the edges of the cluster
			self.limitsLow = limitsLow
			self.limitsHigh = limitsHigh
			self.findKeys()
			self.computeWeight()
			self.computeSum()
			self.computeSqSum()
			self.computeSSQ()
			self.computeCoG()
	def computeSSQ(self):
		# Compute SSQ of the cluster
		self.SSQ = 0
		if self.weight is not 0:
			for i in xrange(ndim):
				self.SSQ += ( self.sqSum[i] - pow(self.Sum[i],2)/self.weight )
	def computeCoG(self):
		# Find center of gravity of the cluster 
		if self.weight > 0:
			self.CoG = tuple(x/self.weight for x in self.Sum)
		else:
			self.CoG = False
	def findKeys(self):
		keys = keylist[:]
		for i in xrange(ndim):
			ckeys = [key for key in keys if key[i] >= self.limitsLow[i] and key[i] <= self.limitsHigh[i]]
			keys = ckeys
		self.keys = keys[:]
	def computeSqSum(self):
		# For each dimension, compute sum of square coordinates
		self.sqSum = [0]*ndim
		self.sqSum = list(self.sqSum)
		for i in xrange(ndim):
			for key in self.keys:
				# A bin with weight W counts as W times that bin
				self.sqSum[i] += ( pow(key[i],2) * dataSet[key] )
		# Make the vector sum of squares an immutable object
		self.sqSum = tuple(self.sqSum)
	def computeWeight(self):
		# Compute content of cluster 
		self.weight = 0
		for key in self.keys:
			self.weight += dataSet[key]
	def computeSum(self):
		# Compute vector sum of the cluster
		self.Sum = [0]*ndim
		for i in xrange(ndim):
			for key in self.keys:
				# A bin with weight W counts as W times that bin
				self.Sum[i] += (key[i] * dataSet[key])
		# Make the vector sum an immutable object
		self.Sum = tuple(self.Sum)

		

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
		self.queue.append((node, node.clusterData.SSQ))
	def get(self):
		# Returns the cluster with the highest SSQ
		self.queue = sorted(self.queue, key = lambda x: x[1])
		return self.queue[-1]
	def delete(self, cl):
		self.queue.remove(cl)
	def makeListOfClusters(self):
		# For the bottom-up part of the algorithm
		self.listOfClusters = [item[0].clusterData for item in self.queue]
	def addCluster(self, cl):
		self.listOfClusters.append(cl)
	def deleteCluster(self, cl):
		self.listOfClusters.remove(cl)

def computeDeltaSSQ(c1, c2):
	""" Variation in SSQ due to the splitting of a cluster or to the merging of two sub-clusters"""
	dssq = 0
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
	refDssq = 0
	leftLimitsLow = [None]*ndim
	leftLimitsHigh = [None]*ndim
	rightLimitsHigh = [None]*ndim
	rightLimitsLow = [None]*ndim
	subClusters = [None]*2
	# Select only the range of the current cluster
	keys = cl.keys[:]
	keys = list(set(keys))
	for i in xrange(ndim):
		for j in xrange(1,len(keys)):
			if keys[j][i] != keys[j-1][i]:
				#print keys[j][i], keys[j-1][i], type(keys[j][i])
				leftLimitsLow = cl.limitsLow[:]
				leftLimitsHigh = cl.limitsHigh[:]
				leftLimitsHigh[i] = 0.5 * (keys[j-1][i] + keys[j][i])
				rightLimitsLow = cl.limitsLow[:]
				rightLimitsHigh = cl.limitsHigh[:]
				rightLimitsLow[i] = 0.5 * (keys[j-1][i] + keys[j][i])
				# Create the two subclusters from the mother
				left = newCluster(leftLimitsLow, leftLimitsHigh)
				right = newCluster(rightLimitsLow, rightLimitsHigh)
				# Compute weighted Delta SSQ
				if left.weight>0 and right.weight>0:
					newDssq = pow(computeDeltaSSQ(left, right), power)
					if (newDssq > refDssq) and (newDssq > avgDeltaSSQ):
						#print computeDeltaSSQ(left,right), newDssq, avgDeltaSSQ, i, j, leftLimitsLow, leftLimitsHigh, rightLimitsLow, rightLimitsHigh
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
	avgDeltaSSQ = root.clusterData.SSQ / root.clusterData.weight
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
		# Confirm the merging only if minSSQincrease < avgDeltaSSQ
		if minSSQincrease < avgDeltaSSQ:
			#print "minSSQincrease: %s avgDeltaSSQ: %s"%(minSSQincrease, avgDeltaSSQ)
			largerCluster = mergeClusters(c1, c2)
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
		else:
			done = True
	print "Micro-clusters regrouped into %s clusters."%len(pq.listOfClusters)
	return True

def findClusters():
	""" Main program: returns a list of clusters in pq """
	if not initStructures():
		print "Error: initStructures"
		sys.exit(1)
	if not topDownSplitting():
		print "Error: topDownSplitting"
		sys.exit(1)
	if not bottomUpClustering():
		print "Error: bottomUpClustering"
		sys.exit(1)
	# Return the centers of gravity of found clusters
	listOfCoGs = []
	for cl in pq.listOfClusters:
		listOfCoGs.append((cl.CoG, cl.weight))
	return listOfCoGs

def CLUBSclustering(data, nd):
	""" Execute findClusters() and handle I/O """
	global dataSet, ndim, power, pq, done, listOfMergeablePairs
	dataSet = data
	ndim = nd
	pq = newPriorityQueue()
	done = False
	listOfMergeablePairs = []
	return findClusters()

	
if __name__ == '__main__':
	""" Sample usage with toy clusters """
	dim = 50
	DS = dict([((x,y),0) for x in range(dim) for y in range(dim)])	
	# 10 clusters 
	borders = 4, 45
	for i in range(10):
		center = [0]*2
		center[0] = randint(borders[0], borders[1])
		center[1] = randint(borders[0], borders[1])
		cdim = randint(3,15)
		for j in xrange(cdim):
			DS[(randint(center[0]-3,center[0]+3),randint(center[1]-3,center[1]+3))] = randint(3,7)

	x,y,w = [],[],[]
	for item in DS.keys():
		if DS[item] > 0:                  
			x.append(item[0])
			y.append(item[1])
			w.append(DS[item])
	plt.ion()
	#plt.scatter(x,y,c=w,s=200)
	plt.scatter(x,y,c=w)
	plt.gray()
	plt.grid()
	plt.pause(0.001)
	clusters = CLUBSclustering(DS,2)
	print clusters

	xc, yc, wc = [],[],[]
	for i in xrange(len(clusters)):
		xc.append(clusters[i][0][0])
		yc.append(clusters[i][0][1])
		wc.append(clusters[i][1])
	#area = [200*np.sqrt(w) for w in wc]
	#area = [np.sqrt(w) for w in wc]
	area = wc

	normal = pylab.normalize(min(wc), max(wc))
	colors = pylab.cm.jet(normal(wc))
	rectangles = []
	for c,cl in zip(colors,pq.listOfClusters):
		bottomleft = cl.limitsLow
		extendright = cl.limitsHigh[0]-cl.limitsLow[0]
		extendtop = cl.limitsHigh[1]-cl.limitsLow[1]
		rectangles.append(matplotlib.patches.Rectangle(bottomleft,extendright,extendtop, color=c, alpha = 0.3))

	for r in rectangles:
		plt.gca().add_patch(r)
	plt.axis([0,dim,0,dim])
	plt.pause(0.001)

	plt.scatter(xc,yc,c='red',s=area,alpha=0.9)
	plt.pause(0.001)


