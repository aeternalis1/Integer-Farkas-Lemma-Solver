from .hermite import ConvertToHNF
from .inout import GetInput
from .utils import LCM, ConvertToRatMat
from .models import Result
from copy import deepcopy

import numpy as np
import sympy as sy


def Solve(M):
	m, n = M.shape

	# scale array to integers
	den = 1
	for i in range(m):
		for j in range(n):
			den = LCM(den, M[i][j].d)

	ogMat = deepcopy(M)
	intM = np.array([[M[i][j].n * den // M[i][j].d for j in range(n)] for i in range(m)])
	A = ConvertToRatMat(intM[:, :-1])

	# convert scaled matrix to HNF
	M, MOps, MMats = ConvertToHNF(intM[:, :-1])

	b = ConvertToRatMat(intM[:, -1])
	B = M[:, :m]

	# candidate solution for HNF form is B^{-1}b
	invB = ConvertToRatMat(B).inv()
	invBb = invB * b
	sol = sy.Matrix(list(invBb) + [sy.Rational('0') for x in range(n-m-1)])
	# if any member of sol is not integer, then we identify which row
	# in B^{-1} gave y^T B integer and y^T b not integer
	for i in range(m):
		if sol[i] != int(sol[i]):
			barY = invB.row(i)
			barYM = barY * A
			barYb = barY * b
			res = Result(
				1,
				MOps,
				ogMat,
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
			return res

	# recover a solution using sol via back substitution
	ogSol = deepcopy(sol)
	sols = [deepcopy(sol)]
	for op in MOps[::-1]:
		sol = op.reverse(sol)
		sols.append(deepcopy(sol))
	sols = sols[::-1]

	res = Result(
		0,
		MOps,
		ogMat,
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