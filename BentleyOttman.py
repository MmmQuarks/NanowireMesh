import sortedcontainers as sc
import random
import numpy as np
import NanowireMesh as nwm
import itertools
import pdb
from copy import copy, deepcopy

# annotaed program outline
# Input: a set S of line segments in the plane
# Output: the set of intersection points amont the segments in S and the
# segments that contain those points

# 1. Initialize empty event queue Q. Insert segment endpoints into Q.
#	store corresponding segment and whether endpoint is top or bottom endpoint
# 2.  initialize empty status structure T
# 3. while Q is not empty
	# 4. get next event point p in Q and delete it
	# 5. handle_event_point(p)

# handle_event_point(p) description
	# 1. Let U(p) be set of segments in Q whose upper endpoint is p
		# for horiz segments upper endpoint is left one
	# 2. find all segments in T that contain p.
	#	L(p) = segments in T with lower point = p
	#	C(p) = segments in T containing p on interior
	# 3. If union of L(p), U(p), C(p) contains more than one segment
	#	4. THEN report p as intersection along with the line segments in L(p), U(p), C(p)
	# 5. Delete segments in L(p) union C(p) from T
	# 6. insert segments in U(p) union C(p) into T
	#	order should be order in which they are intersected by sweep
	#	line just below p. If there is horiz segment it comes last.
	#	7. deleting and re-inserting in this way changes the order correctly 
	# 8. If U(p) union C(p) is empty:
		# 9. Then Let s_l and s_r be the left and right neighbors 
		#	of p in T
		# 10. find_new_event(s_l, s_r, p)
	# 11. else:  let s' be the leftmost segment of U(p) union C(p) in T
		# 12. let s_l be the left neighbor of s' in T
		# 13. find_new_event(s_l, s', p)
		# 14. let s'' be the rightmost segment of U(p) uion C(p) in T
		# 15. let s_r be the right neighbor of s'' in T
		# 16. find_new_events(s'', s_r, p)

# find_new_event(s_l, s_r, p)  description
	# 1. if (s_l and s_r intersect below the sweep line)  OR (they intersect on it AND to the right of p)
		# AND if the intersection point is not in Q as an event
			# 2. THEN insert intersection point into Q as event



def handle_event_point(p, Q, T, dictOfIntersections, dictOfEndpoints, probDict):
	# p is of type EventPoint
	# move sweep line to be at the height of p
	T.sweepLine = np.array([(T.sweepLine[0][0], p.y), (T.sweepLine[1][0], p.y)])
	segmentsInT = sc.SortedSet({segment.segment for segment in T})
	segs = segmentsInT
	ins = probDict['inserted']
	rem = probDict['removed']
	# 1. let U(p) be the set of segments in Q whose uper endpoint is p
	U = set(p.upperSegments)
	pNumpy = np.array((p.x, p.y))
	# 2. find all segments in T that contain p
	# L(p) is segments in T with lower point p
	L = set(p.lowerSegments)
	LWithSegmentData = {activeSegment for activeSegment in T if activeSegment.segment in L}
	# C(p) is segments in T that contain p on interior

#	
	CWithSegmentData = {seg for seg in T if seg.contains(p.x, p.y) and seg.segment not in L}

#	# making just the segments
	C = {activeSegment.segment for activeSegment in CWithSegmentData}

	# 3. If union LUC has more than one point, report p as intersection
	# along with all segments in LUC
	LUC = L | U | C
	print(LUC)
	if len(LUC) > 1:
		for pair in itertools.combinations(LUC, 2):
			#4. reporting intersections
			dictOfIntersections[pair] = p

	# 5. delete segments in L and C from T
	for segment in LWithSegmentData | CWithSegmentData:
		T.remove(segment)

	# 6. insert segments in U and C into T
	# should be sorted such that they are intersected by sweep line just below p
	# this is done at the (***)


	# then we convert the segments in U to be ActiveSegment objects
	UWithSegmentData = { ActiveSegment(segment, dictOfEndpoints[segment], sweepLine = T.sweepLine) for segment in U}
	# 7. then we add them in the status structure
	T.update(UWithSegmentData | CWithSegmentData)

	# (***) - sorting by point just below p
	epsilon = np.array((0, 10**(-8)))
	T.sweepLine = T.sweepLine - epsilon
	# this makes sure that the segments are sorted according to their ordering
	# at a point just below p
	# this also makes sure that the point p is never exactly in the status structure
	# so its left and right neighbors are unambiguously
	# left = T.bisect_left(activeSegmentP)
	# right = T.bisect_right(activeSegmentP)  = left + 1
	
	# 8. if the union of U and P is empty
	if U | C == {} or U | C == set():
		# 9. find left and right neighbors of p in T
		pNumpy = np.array((p.x, p.y))
		endpointsP = np.array([pNumpy + 3 * epsilon, pNumpy - 3 * epsilon])
		activeSegmentP = ActiveSegment(segment = -1, endpoints = endpointsP, sweepLine = T.sweepLine) 
		# making sure that T has elements to examine
		if len(T) > 1:
			# make sure that the left and right neighbors of activeSegmentP
			# exist in T (by making sure that p is between first and last elems of T)
			if T[0] < activeSegmentP < T[-1]:
				# get the index of the point p if it were to be inserted into T
				l_index = T.bisect_left(activeSegmentP) - 1
				r_index = l_index + 1
				s_l = T[l_index]
				try:
					s_r = T[r_index]
				except IndexError:
					pdb.set_trace()
				Q = find_new_event(s_l, s_r, p, T, Q)
	else:
		# if the union U, C isn't empty, do this
		if len(T) > 1:
			UCSorted = sc.SortedList(UWithSegmentData | CWithSegmentData)
			# 11. let s' be the leftmost segment of union(U,C) in T
			sp = UCSorted[0]

			# 12. let s_l be the left neighbor of s' in T
			s_l_index = T.bisect_left(sp) - 1 if sp in T else T.bisect_left(sp)
			# if sp in T, then bisect_left
			# will return the index of sp.
			# therefore the index of the left neighbor
			# of sp will be bisect_left(sp) - 1
			# if sp not in T, then sp will lie between
			# two values, a and b, with a < b.
			# Note that a or b or both may not exist.
			# In this case, bisect_left(sp) will return the index
			# of a, which is the correct left neighbor of sp.

			# there are multiple members of T
			# the index returned isn't the rightmost member of 

			# making sure that the left neighbor of sp exists
			if s_l_index >= 0:
				s_l = T[s_l_index]
				# 13. 
				Q = find_new_event(s_l, sp, p, T, Q)

			#14. let s'' be the rightmost segment of union(U,C) in T
			spp = UCSorted[-1]

			# 15. let s_r be the right neighbor of s'' in T
			s_r_index = T.bisect_right(spp) + 1 if spp in T else T.bisect_right(spp)
			# if spp in T, then bisect_right will return
			# the index of spp. So in this case the right neighbor of 
			# spp is bisect_right(spp) + 1
			# if spp not in T, it will reside between two values,
			# call them a and b (possible one or both won't exist)
			# let a < b.
			# then bisect_right(spp) will give the index of b which is 
			# correct right neighbor of spp. So in this case 
			# right neighbor of spp is bisec_rich(spp)
			# with no plus 1.
			
			# making sure that the right neighbor of sp exists
			if s_r_index < len(UCSorted):
				s_r = T[s_r_index]
				#16.
				Q = find_new_event(spp, s_r, p, T, Q)
	return Q, T, dictOfIntersections, probDict

def find_new_event(s_l, s_r, p, T, Q):
	# if s_l and s_r intersect below the sweep line or on it and to the right
	# of the current event point p
	seg0 = np.array((s_l.upperEndpoint, s_l.lowerEndpoint))
	seg1 = np.array((s_r.upperEndpoint, s_r.lowerEndpoint))
	intersection = _find_intersection(seg0, seg1)
	#making sure that the intersection exists
	# which means its type is not string
	if type(intersection) != str:
		intersectBelow = intersection[1] < T.sweepLine[0][1]
		intersectOnAndToRight = intersection[1] == T.sweepLine[0][1] and intersection[0] > p.x
		if intersectBelow or intersectOnAndToRight:
			# check to make sure event is not already in Q
			newEventPoint = EventPoint(x = intersection[0], y = intersection[1], kind = 'int')
			if newEventPoint.coordinates not in Q:
				print(newEventPoint)
				Q.insert(newEventPoint)
			else:
				print('discarded point because already in Q', intersection)
				print(Q)
	return Q


class EventPoint:
	# subclass of EventPoint
	class __EventPointCoordinates:
		# This is a class that overloads the < operator so that the binary
		# search tree class can propertly sort event points
		def __init__(self, x, y):
			self._x = x
			self._y = y
		
		def get_x(self):
			return self._x

		def get_y(self):
			return self._y

		# making these properties and therefore less editable
		# using https://inventwithpython.com/blog/2019/02/01/hashable-objects-must-be-immutable/
		# and https://www.programiz.com/python-programming/property
		x = property(fget = get_x)
		y = property(fget = get_y)

		#this ordering means that points are sorted from top to bottom
		# and then from left to right
		def __lt__(self, otherPoint):
			# self < otherPoint (earlier in the queue) if self is above otherPoint.
			# or if they are at the same height if self is left of otherPoint
			if self._y > otherPoint._y:
				return True
			elif self._y == otherPoint._y and self._x < otherPoint._x:
				return True
			else:
				return False

		def __repr__(self):
			return ''.join(['(x = ', str(self.x), ', y = ', str(self.y), ')'])

		def __eq__(self, other):
			if not isinstance(other, type(self)):
				return False
			else:
				return (self._x == other._x) and (self._y == other._y)
	
		# re-using the hash function from tuples
		def __hash__(self):
			return hash( (self._x, self._y) )


	def __init__(self, x, y, upperSegments = [], lowerSegments = [], kind = None):
		# adding coordinates as an __EventPoint object
		# so dict sorting can happen
		self.coordinates = EventPoint.__EventPointCoordinates(x = x, y = y)

		def _to_list(item):
			# converting item to list if it isn't already
			try:
				itemIterable = iter(item)
			except TypeError:
				return [item]
			else:
				return list(item)

		self.upperSegments = _to_list(upperSegments)
		self.lowerSegments = _to_list(lowerSegments)
		self.kind = kind




	def get_x(self):
		return self.coordinates.x

	def get_y(self):
		return self.coordinates.y

	# making the x and y coordinates not editable
	x = property(fget = get_x)
	y = property(fget = get_y)

	def __repr__(self):
		return str(self)

	def __str__(self):
		return ''.join(['EventPoint(x = ', str(self.x),', y = ', str(self.y),
				', upperPtSegs = ', str(self.upperSegments),', lowerPtSegs = ', str(self.lowerSegments),
				', kind = ', str(self.kind), ')'])
	


class EventPointQueue(sc.SortedDict):
	def __init__(self):
		super().__init__(self)

	def insert(self, newEventPoint):
		# handling the insertion of new event points
		if newEventPoint.coordinates in self:
			# if this point exists, append segments to upperSegments
			# note that if newEventPoint is not an upper segment
			# we will append [] which has no effect
			self[newEventPoint.coordinates] += newEventPoint.upperSegments
			self[newEventPoint.coordinates] += newEventPoint.lowerSegments
		else:
			# if this is a new event point, add it
			self[newEventPoint.coordinates] = newEventPoint

	def __repr__(self):
		points = [str(point) for point in self.values()]
		return 'EventPointQueue(' + '\n'.join(points) + ')'

class StatusStructure(sc.SortedSet):
	def __init__(self, sweepLine):
		super().__init__()
		self._sweepLine = sweepLine

	def get_sweep_line(self):
		return self._sweepLine
	def set_sweep_line(self, newSweepLine):
		self._sweepLine = newSweepLine
		# setting the sweep lines of all the active segments
		# this will automatically recalculate all the sweepLineCrossings
#		for activeSegment in self:
#			# remove the item from the list
##			self.discard(activeSegment)
#			# set the new sweep line value
#			activeSegment.sweepLine = newSweepLine
#			# re-insert the item into the list
##			self.add(activeSegment)
		segments = []
		while self:
			# get non self container of segments
			segments.append(self.pop())
		for segment in segments:
			# set the new sweepline
			segment.sweepLine = newSweepLine
			# add thesegment back into status structure
			self.add(segment)

	def insert(self, newSegment):
		if newSegment in self:
			pass
		else:
			self.add(newSegment)
	def __repr__(self):
		segments = [str(segment) for segment in self]
		return 'StatusStructure(' + '\n'.join(segments) + ')'

	sweepLine = property(fget = get_sweep_line, fset = set_sweep_line)

class ActiveSegment(object):
	def __init__(self, segment, endpoints, sweepLine):
		self.segment = segment
		x1, y1 = endpoints[0]
		x2, y2 = endpoints[1]
		# assigning the endpoints names as they are sorted
		# in the EventPoint definition
		if y2 > y1 or (y2 == y1 and x2 < x1):
			self.upperEndpoint = np.array((x2, y2))
			self.lowerEndpoint = np.array((x1, y1))
		else:
			self.upperEndpoint = np.array((x1, y1))
			self.lowerEndpoint = np.array((x2, y2))

		# calculating the sweep line point
		self.sweepLineCrossing = _find_intersection([self.upperEndpoint, self.lowerEndpoint], sweepLine)

		# setting sweepLine
		self._sweepLine = sweepLine
	
	def get_sweep_line(self):
		return self._sweepLine

	def set_sweep_line(self, newSweepLine):
		self._sweepLine = newSweepLine
		self.sweepLineCrossing = _find_intersection([self.upperEndpoint, self.lowerEndpoint], self.sweepLine)

	sweepLine  = property(fget = get_sweep_line, fset = set_sweep_line)

	# tests whether the segment contains the point (x,y)
	def contains(self,x,y):
		# test if the segments are aligned
		endpointsDisplacement = self.upperEndpoint - self.lowerEndpoint
		pointDisplacement = np.array((x,y)) - self.lowerEndpoint
		# calculate cross product
		cross = np.cross(endpointsDisplacement, pointDisplacement)
		if np.linalg.norm(cross) <= 10**(-10):
			# okay they're parallel
			# now test to make sure that the point lies on 
			# the finite extent of the segment
			minX = min(self.lowerEndpoint[0], self.upperEndpoint[0])
			maxX = max(self.lowerEndpoint[0], self.upperEndpoint[0])
			if minX < x < maxX:
				minY = min(self.lowerEndpoint[1], self.upperEndpoint[1])
				maxY = max(self.lowerEndpoint[1], self.upperEndpoint[1])
				if minY <  y < maxY:
					return True
		return False


	def __lt__(self, otherSegment):
		# first handle what happens if segments are equal
		if self == otherSegment:
			return False
		# first handle what happens if one or 
		# both of the segments are parallel to the sweep line
		elif type(self.sweepLineCrossing) == str and type(otherSegment.sweepLineCrossing) == np.ndarray:
			if self.sweepLineCrossing == 'Parallel':
				# then the parallel one is greater so self > other NOT self < other
				return False
		elif type(self.sweepLineCrossing) == np.ndarray and type(otherSegment.sweepLineCrossing) == str:
			if otherSegment.sweepLineCrossing == 'Parallel':
				# the parallel one is greater again so self < other
				return True
		elif type(self.sweepLineCrossing) == str and type(otherSegment.sweepLineCrossing) == str:
			if self.sweepLineCrossing == 'Parallel' and otherSegment.sweepLineCrossing == 'Parallel':
				# both lines are parallel to the sweep line and on sweepline. sort by upper endpoint
				if self.upperEndpoint[0] < otherSegment.upperEndpoint[0]:
					return True
				elif self.upperEndpoint[0] > otherSegment.upperEndpoint[0]:
					return False
				elif self.upperEndpoint[0] == otherSegment.upperEndpoint[0]:
					# if the upper endpoints are the same, compare lower endpoints
					if self.lowerEndpoint[0] < otherSegment.lowerEndpoint[0]:
						return True
					elif self.lowerEndpoint[0] > otherSegment.lowerEndpoint[0]:
						return False
					elif self.lowerEndpoint[0] == otherSegment.lowerEndpoint[0]:
						# if the lower endpoints are also the same, sort by segment number
						# this should never ever ever ever happen
						return self.segment < otherSegment
						# at this point any parallel segments are well ordered
		elif type(self.sweepLineCrossing) == np.ndarray and type(otherSegment.sweepLineCrossing) == np.ndarray:
			# if both segments have points for their line crossing, we sort by which point is farther left
			if self.sweepLineCrossing[0] < otherSegment.sweepLineCrossing[0]:
				return True
			elif self.sweepLineCrossing[0] > otherSegment.sweepLineCrossing[0]:
				return False
			else:
				# if the points are equal, sort by segment number
				return self.segment < otherSegment.segment
		else:
			# at this point, at least one of the segments does not cross the sweep line
			# this should never happen, so we raise an error
			raise Exception('Line segment does not cross sweep line. Something went wrong somewhere.')

	def __le__(self, otherSegment):
		if self < otherSegment or self == otherSegment:
			return True
		return False

	def __repr__(self):
		return ''.join(['ActiveSegment(seg = ', str(self.segment),
				', upper = ', str(self.upperEndpoint), 
				', lower = ', str(self.lowerEndpoint), 
				', sweepLineCrossing = ', str(self.sweepLineCrossing), ')'])

	# crude but it works
	def __eq__(self, other):
		return self.segment == other.segment

	def __hash__(self):
		return hash((self.segment,))

def _find_intersection(seg0, seg1):
	# the segments should be given as a list of tuples which are the endpoints of the lines
	tolerance = 10**(-10)
	x1, y1 = seg0[0]
	x2, y2 = seg0[1]

	x3, y3 = seg1[0]
	x4, y4 = seg1[1]

	# make the coefficients matrix
	A = np.array([ (-(y2 - y1), x2 - x1),
			(-(y4 - y3), x4 - x3)])
	b = np.array([ -(y2 - y1) * x1 + (x2 - x1) * y1,
			-(y4 - y3) * x3 + (x4 - x3) * y3])
	
	try:
		x, y = np.linalg.solve(A, b)
	except np.linalg.LinAlgError:
		# means lines are parallel
		# return empty tuple
		return 'Parallel'
	else:
		# check to make sure result is on finite extent of lines
		if min(x1,x2) - tolerance <= x <= max(x1,x2) + tolerance: 
			if min(y1,y2) - tolerance <= y <= max(y1,y2) + tolerance:
				if min(x3,x4) - tolerance <= x <= max(x3,x4) + tolerance:
					if min(y3, y4) - tolerance <= y <= max(y3, y4) + tolerance:
						return np.array((x,y))
		return 'Solution is off segments'
	






def find_intersections(dictOfEndpoints, xMin, xMax, yMax):
	# 1. initializing sorted list that will serve as event queue
	Q = EventPointQueue()
	for segment, endpoints in dictOfEndpoints.items():
		eventPoint0  = EventPoint(x = endpoints[0][0],
						y = endpoints[0][1])
		eventPoint1 = EventPoint(x = endpoints[1][0],
						y = endpoints[1][1])
		# adding the segment number to only top points
		if eventPoint0.coordinates < eventPoint1.coordinates:
			eventPoint0.upperSegments = [segment]
			eventPoint1.lowerSegments = [segment]
		elif eventPoint0.coordinates > eventPoint1.coordinates:
			eventPoint1.upperSegments = [segment]
			eventPoint0.lowerSegments = [segment]
		else:
			exceptionMsg = 'The segment ' + str(segment) + ' has identical endpoints.'
			raise Exception(exceptionMsg)

		# adding event points to our event queue
		Q.insert(eventPoint0)
		Q.insert(eventPoint1)
	# 2. initialize empty status structure T
	sweepLine = np.array([(xMin, yMax), (xMax, yMax)])
	T = StatusStructure(sweepLine)

	# initialize container for intersections
	dictOfIntersections = dict()

	# initialize troubleshooting object
	probDict = {'inserted' : sc.SortedSet(),
			'removed' : sc.SortedSet()}

	# 3. while Q is not empty
	while Q:
		# 4. determine next point p in Q and delete it
		keysList = list(Q.keys())
		keyToRemove = keysList[0]
		p = Q.pop(keyToRemove)
		
		# 5. handle event point
		Q, T, dictOfIntersections, probDict = handle_event_point(p, Q, T, dictOfIntersections, dictOfEndpoints, probDict = probDict)
		
	return dictOfIntersections

def main():
	sequence = [EventPoint(x = random.randint(0,10), y = random.randint(0,10), upperSegments = []) for n in range(10)]
	sd = sc.SortedDict()
	for item in sequence:
		sd[item.coordinates] = item.upperSegments


if  __name__ == '__main__':
	main()
