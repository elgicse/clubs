#!/usr/bin/env python
#
# Implementation of the M-CLUBS algorithm
# by E. Masciari, G.M.Mazzeo and C. Zaniolo
# ("A New, Fast and Accurate Algorithm for
# Hierarchical Clustering on Euclidean Distances",
# PAKDD 2, volume 7819 of Lecture Notes in Computer Science,
# page 111-122. Springer, 2013)
#
# Author: Elena Graverini
# Contact: elena.graverini@cern.ch

import numpy as np
import sys

dataSet = None
ndim = None
power = 0.8
pq = None
done = False
avgDeltaSSQ = None
listOfMergeablePairs = None

class newCluster():
	""" Class holding cluster data (edges and statistics) """
	def __init__(self, limitsLow, limitsHigh, status = None):
		""" Either merge two clusters (bottom-up) or create a new micro-cluster (top-down) """
		if status is "merging":
			first = limitsLow
			second = limitsHigh
			self.limitsLow = None
			self.limitsHigh = None
			self.weight = first.weight + second.weight
			self.Sum = first.Sum + second.Sum
			self.sqSum = 0
			self.CoG = 0
			self.SSQ = 0
		else:
			self.limitsLow = limitsLow
			self.limitsHigh = limitsHigh
			self.weight = computeWeightIn(limitsLow, limitsHigh)
			self.Sum = computeSum(limitsLow, limitsHigh)
			self.sqSum = 0
			self.CoG = 0
			self.SSQ = 0
	def computeSSQ(self):
		""" Compute SSQ of the cluster """
		self.computeSqSum()
		for i in xrange(ndim):
			self.SSQ += ( self.sqSum[i] - pow(self.Sum[i],2)/self.weight )
		return self.SSQ
	def computeCoG(self):
		""" Find center of gravity of the cluster """
		self.CoG = self.Sum / self.weight
		return self.CoG
	def computeSqSum(self):
		""" For each dimension, compute sum of square coordinates """
		pass

def computeWeightIn(limitsLow, limitsHigh):
	""" Compute content of cluster """
	pass

def computeSum(limitsLow, limitsHigh):
	""" Compute vector sum of the cluster """
	pass

class newTreeNode():
	""" Auxiliary binary tree for top-down splitting """
	def __init__(self, limitsLow, limitsHigh):
		self.left = None
		self.right = None
		self.clusterData = newCluster(limitsLow, limitsHigh)
	def setChildren(self, left, right):
		self.left = left
		self.right = right

class newPriorityQueue():
	""" Selects the cluster with the highest SSQ """
	def __init__(self):
		self.queue = []
	def add(self,cl):
		self.queue.append((cl,cl.SSQ))
	def get(self):
		self.queue = sorted(self.queue, key = lambda x: x[1])
		return self.queue[-1]
	def delete(self, cl):
		self.queue.remove(cl)

def computeDeltaSSQ(c1, c2):
	""" Variation in SSQ due to the splitting of a cluster or to the merging of two sub-clusters"""
	dssq = 0
	factor = (c1.weight * c2.weight) / (c1.weight + c2.weight)
	for i in xrange(ndim):
		dssq += pow( c1.Sum[i]/c1.weight - c2.Sum[i]/c2.weight, 2 )
	return factor * dssq

def areAdjacent(c1, c2):
	""" Find out if two clusters are adjacent or overlapping """
	pass

def findClustersToMerge():
	""" Find the pair of adjacent clusters that yiealds the least SSQ increase """
	global listOfMergeablePairs
	listOfMergeablePairs = sorted(listOfMergeablePairs, key = lambda x: x[2])
	bestPair = listOfMergeablePairs[0]
	return bestPair[0], bestPair[1]

def mergeClusters(c1, c2):
	""" Merge two adjacent clusters """
	return newCluster(c1, c2, "merging")

def initStructures():
	""" Initialize priority queue and binary tree """
	global avgDeltaSSQ, pq
	marginsLow = "MARGINI INFERIORI DEL DATA SET NELLE N DIMENSIONI"
	marginsHigh = "MARGINI SUPERIORI DEL DATA SET NELLE N DIMENSIONI"
	root = newTreeNode(marginsLow, marginsHigh)
	pq.add(root)
	avgDeltaSSQ = root.clusterData.computeSSQ() / root.clusterData.weight
	return True

def splitCluster(cl, avgDeltaSSQ):
	""" Returns the margins of the two children or False """
	pass

def topDownSplitting():
	""" Split the data set into micro-clusters """
	global pq, done
	while not done:
		currentCluster = pq.get()
		children = splitCluster(currentCluster, avgDeltaSSQ)
		if not children:
			done = True
		else:
			currentCluster.setChildren(children)
			pq.delete(currentCluster)
			pq.add(children[0])
			pq.add(children[1])
	return True

def bottomUpClustering():
	""" Merge micro-clusters to provide the final clustering """
	global listOfMergeablePairs
	for c1 in pq:
		for c2 in pq:
			if (c2 is not c1) and areAdjacent(c1, c2):
				dssq = computeDeltaSSQ(c1, c2)
				listOfMergeablePairs.append( (c1, c2, dssq) )
	clustersToMerge = findClustersToMerge() # a pair of two clusters
	minSSQincrease = computeDeltaSSQ(clustersToMerge)
	while minSSQincrease < avgDeltaSSQ:
		largerCluster = mergeClusters(clustersToMerge)
		listOfMergeablePairs = [item for item in listOfMergeablePairs if item[0] is clustersToMerge[0] or item[0] is clustersToMerge[1] or item[1] is clustersToMerge[0] or item[1] is clustersToMerge[1]]
		pq.delete(clustersToMerge[0])
		pq.delete(clustersToMerge[1])
		pq.add(largerCluster)
		for cl in pq:
			if (cl is not largerCluster) and areAdjacent(largerCluster, cl):
				dssq = computeDeltaSSQ(largerCluster, cl)
				listOfMergeablePairs.append( (largerCluster, cl, dssq) )
		clustersToMerge = findClustersToMerge()
		minSSQincrease = computeDeltaSSQ(clustersToMerge)
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
	return pq

def CLUBSclustering(dataSet, ndim):
	""" Execute findClusters() and handle I/O"""
	global dataSet, ndim, power, pq, done, avgDeltaSSQ, listOfMergeablePairs
	dataSet = dataSet
	ndim = ndim
	power = 0.8
	pq = newPriorityQueue()
	done = False
	avgDeltaSSQ = 0
	listOfMergeablePairs = []
	return findClusters()