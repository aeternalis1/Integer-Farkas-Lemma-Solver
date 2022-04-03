import math
import sympy as sy


def LCM(a, b):
	return a * b // math.gcd(a,b)


def ConvertToRatMat(M):
	# takes numpy array M and converts to sympy array of rationals
	try:
		m, n = M.shape
		return sy.Matrix([[sy.Rational(str(M[i][j])) for j in range(n)] for i in range(m)])
	except ValueError:
		m = M.shape[0]
		return sy.Matrix([sy.Rational(str(M[i])) for i in range(m)])