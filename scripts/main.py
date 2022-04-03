from .hermite import ConvertToHNF
from .inout import GetInput
from .utils import LCM, ConvertToRatMat
from .models import Result
from copy import deepcopy

import numpy as np
import sympy as sy

a = [
	[2, 1, 3, 4],
	[6, 7, 8, 9],
	[11, 12, 13, 14]
]

M = np.array(a)

b = [
	[2, 0, 0, 0],
	[6, 7, 0, 0],
	[11, 12, 13, 0]
]

M2 = np.array(b)


def Solve(M):
	m, n = M.shape

	# scale array to integers
	den = 1
	for i in range(m):
		for j in range(n):
			den = LCM(den, M[i][j].d)

	intM = np.array([[M[i][j].n * den // M[i][j].d for j in range(n)] for i in range(m)])
	A = ConvertToRatMat(intM[:, :-1])

	# convert scaled matrix to HNF
	M, MOps, MMats = ConvertToHNF(intM[:, :-1])

	b = ConvertToRatMat(intM[:, -1])
	B = M[:, :m]

	# candidate solution for HNF form is B^{-1}b
	invB = ConvertToRatMat(B).inv()
	print (invB)
	invBb = invB * b
	sol = sy.Matrix(list(invBb) + [sy.Rational('0') for x in range(n-m-1)])
	print (sol)
	# if any member of sol is not integer, then we identify which row
	# in B^{-1} gave y^T B integer and y^T b not integer
	for i in range(m):
		print (i, sol[i])
		if sol[i] != int(sol[i]):
			barY = invB.row(i)
			print ("No integer solution exists")
			print (barY)
			barYM = barY * A
			print (barYM, len(barYM))
			print (A, len(A))
			barYb = barY * b
			print (barYb)
			res = Result(
				1,
				MOps,
				MMats,
				A,
				B,
				invB,
				b,
				invBb,
				sol,
				sol,
				i,
				None,
				barY,
				barYM,
				barYb
			)
			print (res.ops, len(res.ops))
			return res

	# recover a solution using sol via back substitution
	ogSol = deepcopy(sol)
	sols = [deepcopy(sol)]
	for op in MOps[::-1]:
		sol = op.reverse(sol)
		sols.append(deepcopy(sol))
	sols = sols[::-1]

	print ("Solution exists")
	print (sol)
	print (np.matmul(intM[:, :-1], sol))
	res = Result(
		0,
		MOps,
		MMats,
		A,
		B,
		invB,
		b,
		invBb,
		ogSol,
		sol,
		None,
		sols,
		None,
		None,
		None
	)
	return res


if __name__ == '__main__':
	M = GetInput()
	Solve(M)

	'''
	print (IsInHNF(M))
	print (IsInHNF(M2))

	M, MOps = ConvertToHNF(M2)
	for i in M:
		print (i)
	'''