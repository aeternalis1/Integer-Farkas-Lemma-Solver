import numpy as np
from .models import Operation
from fractions import Fraction
from copy import deepcopy
import sys
from random import randint

sys.setrecursionlimit(1000000000)

def ConvertToHNF(M, bmark=False):
	'''
	Accepts an integer m by n matrix M of full row rank and returns
	the result of converting it to HNF, as well a list of operations
	and the intermediate matrices
	'''
	ops = []
	m, n = M.shape
	mats = [deepcopy(M)]
	maxInt = 0
	for i in range(m):
		# B is currently i x i matrix in top left
		# D is (m-i) x (n-i) matrix in bottom right
		for j in range(n):
			if bmark:
				maxInt = max(maxInt, abs(M[i,j]))
		for j in range(i, n):
			# convert all columns to have d[1][k] nonnegative, k > i
			if M[i][j] < 0:
				M[:,j] = M[:,j] * -1
				ops.append(Operation(2, j))
				if not bmark:
					mats.append(deepcopy(M))

		# take two nonzero columns, subtract from each other
		# repeatedly to get a 1, and then zero out everything else
		nonzero = []
		for j in range(i, n):
			if M[i][j] != 0:
				nonzero.append(j)

		# swap so d[1][1] is nonzero 
		if nonzero and nonzero[0] != i:
			M[:,[i, nonzero[0]]] = M[:,[nonzero[0], i]]
			nonzero[0] = i
			ops.append(Operation(1, i, nonzero[0]))
			if not bmark:
				mats.append(deepcopy(M))

		# zero out d[1][k] for k > i
		for k in range(1, len(nonzero)):
			while M[i][nonzero[k]]:
				if M[i][i] > M[i][nonzero[k]]:
					m = M[i][i] // M[i][nonzero[k]]
					M[:,i] = M[:,i] - M[:,nonzero[k]] * m
					ops.append(Operation(3, i, nonzero[k], -m))
					if not bmark:
						mats.append(deepcopy(M))
				else:
					m = M[i][nonzero[k]] // M[i][i]
					M[:,nonzero[k]] = M[:,nonzero[k]] - M[:,i] * m
					ops.append(Operation(3, nonzero[k], i, -m))
					if not bmark:
						mats.append(deepcopy(M))

				if not M[i][i]:
					M[:,[i, nonzero[k]]] = M[:,[nonzero[k], i]]
					ops.append(Operation(1, i, nonzero[k]))
					if not bmark:
						mats.append(deepcopy(M))

		# update columns in C to make M[i][i] unique maximum in row
		for k in range(i):
			if M[i][k] < 0:
				m = abs(M[i][k]) // M[i][i] + (1 if M[i][k]%M[i][i] else 0)
				M[:,k] = M[:,k] + M[:,i] * m
				ops.append(Operation(3, k, i, m))
				if not bmark:
					mats.append(deepcopy(M))
			elif M[i][k] >= M[i][i]:
				m = M[i][k] // M[i][i]
				M[:,k] = M[:,k] - M[:,i] * m
				ops.append(Operation(3, k, i, -m))
				if not bmark:
					mats.append(deepcopy(M))

	if bmark:
		return M, ops, maxInt

	return M, ops, mats


def ExtEuclid(a, b):
    x = 0
    y = 1
    lastx = 1
    lasty = 0
    while b != 0:
        quo = a // b
        a, b = b, a % b
        x, lastx = lastx - quo * x, x
        y, lasty = lasty - quo * y, y
    return (a, lastx, lasty)


def GetGCDVec(arr):
	n = len(arr)
	if n == 1:
		return arr[0], [1]
	elif n == 2:
		a, b, c = ExtEuclid(arr[0], arr[1])
		return a, [b, c]

	g = [0 for x in range(2*n)]
	y = [0 for x in range(2*n)]
	z = [0 for x in range(2*n)]
	for i in range(n, 2*n):
		g[i] = abs(arr[i-n])

	def bottomup(n, i):
		if i*2 >= 2*n:
			return
		bottomup(n, i*2)
		bottomup(n, i*2+1)
		if i*2+1 >= 2*n:
			g[i] = g[i*2]
			y[i*2] = 1
			return
		a, b, c = ExtEuclid(g[i*2], g[i*2+1])
		g[i] = a
		y[i*2] = b
		y[i*2+1] = c

	def topdown(n, i):
		if i*2 >= 2*n:
			return
		if i*2+1 >= 2*n:
			z[i*2] = z[i]
			return

		c = z[i] * g[i]
		neg = False
		if c < 0:
			neg = True

		a = y[2*i] * z[i]
		b = y[2*i+1] * z[i]
		if neg:
			a,b,c = -a,-b,-c

		j = (b // g[2*i])
		for k in range(-1, 2):
			nj = j + k

			z[i*2] = a + g[2*i+1] * nj
			z[i*2+1] = b - g[2*i] * nj
			if abs(z[i*2]) <= max(g[2*i+1]//2 + 1, abs(c)) and abs(z[i*2+1]) <= g[2*i]//2 + 1:
				if neg:
					z[i*2] *= -1
					z[i*2+1] *= -1
					a,b,c = -a,-b,-c
				topdown(n, i*2)
				topdown(n, i*2+1)
				return

	bottomup(n, 1)
	z[2] = y[2]
	z[3] = y[3]
	topdown(n, 2)
	topdown(n, 3)
	res = z[n:]
	for i in range(n):
		if arr[i] < 0:
			res[i] *= -1

	return g[1], res


def ConvertToHNF2(M, bmark=False):
	'''
	Accepts an integer m by n matrix M of full row rank and returns
	the result of converting it to HNF, as well a list of operations
	and the intermediate matrices
	'''
	ops = []
	m, n = M.shape
	mats = [deepcopy(M)]
	maxInt = 0
	for i in range(m):
		# B is currently i x i matrix in top left
		# D is (m-i) x (n-i) matrix in bottom right
		for j in range(n):
			if bmark:
				maxInt = max(maxInt, abs(M[i,j]))
		for j in range(i, n):
			# convert all columns to have d[1][k] nonnegative, k > i
			if M[i][j] < 0:
				M[:,j] = M[:,j] * -1
				ops.append(Operation(2, j))
				if not bmark:
					mats.append(deepcopy(M))

		# take two nonzero columns, subtract from each other
		# repeatedly to get a 1, and then zero out everything else
		nonzero = []
		for j in range(i, n):
			if M[i][j] != 0:
				nonzero.append(j)

		# swap so d[1][1] is nonzero 
		if nonzero and nonzero[0] != i:
			M[:,[i, nonzero[0]]] = M[:,[nonzero[0], i]]
			nonzero[0] = i
			ops.append(Operation(1, i, nonzero[0]))
			if not bmark:
				mats.append(deepcopy(M))

		# zero out d[1][k] for k > i
		gcd, vec = GetGCDVec([M[i, j] for j in nonzero])
		for j in range(len(nonzero)):
			if vec[j]:
				vec[j] -= 1
				if j:
					M[:,[i, nonzero[j]]] = M[:,[nonzero[j], i]]
					ops.append(Operation(1, i, nonzero[j]))
					if not bmark:
						mats.append(deepcopy(M))
					vec[0], vec[j] = vec[j], vec[0]
				break

		for j in range(len(nonzero)):
			if vec[j]:
				M[:,i] = M[:,i] + M[:,nonzero[j]] * vec[j]
				ops.append(Operation(3, i, nonzero[j], vec[j]))
				if not bmark:
					mats.append(deepcopy(M))

		for k in range(1, len(nonzero)):
			m = M[i][nonzero[k]] // M[i][i]
			M[:,nonzero[k]] = M[:,nonzero[k]] - M[:,i] * m
			ops.append(Operation(3, nonzero[k], i, -m))
			if not bmark:
				mats.append(deepcopy(M))

		# update columns in C to make M[i][i] unique maximum in row
		for k in range(i):
			if M[i][k] < 0:
				m = abs(M[i][k] // M[i][i])
				M[:,k] = M[:,k] + M[:,i] * m
				ops.append(Operation(3, k, i, m))
				if not bmark:
					mats.append(deepcopy(M))
			elif M[i][k] >= M[i][i]:
				m = M[i][k] // M[i][i]
				M[:,k] = M[:,k] - M[:,i] * m
				ops.append(Operation(3, k, i, -m))
				if not bmark:
					mats.append(deepcopy(M))

	if bmark:
		return M, ops, maxInt

	return M, ops, mats
