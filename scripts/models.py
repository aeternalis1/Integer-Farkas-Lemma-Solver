import sympy as sy

class Rational:
	def __init__(self, n, d):
		self.n = n
		self.d = d

	def __str__(self):
		if self.n == 0 or self.d == 1:
			return str(self.n)
		return "%d/%d" % (self.n, self.d)


class Operation:
	def __init__(self, t, c1, c2 = None, m = None):
		self.type = t
		self.c1 = c1
		self.c2 = c2
		self.m = m

	def reverse(self, sol):
		if self.type == 1:
			sol[self.c1], sol[self.c2] = sol[self.c2], sol[self.c1]
		elif self.type == 2:
			sol[self.c1] *= sy.Rational('-1')
		else:
			sol[self.c2] += sy.Rational(str(self.m)) * sol[self.c1]
		return sol


class Result:
	def __init__(
		self, t, ops, ogMat, mats, A, B, invB, b, invBb, ogSol,
		sol, nonIntInd=None, sols=None, barY=None, barYM=None, barYb=None
	):
		self.t = t
		self.ops = ops
		self.ogMat = ogMat
		self.mats = mats
		self.A = A
		self.B = B
		self.invB = invB
		self.b = b
		self.invBb = invBb
		self.ogSol = ogSol
		self.sol = sol 
		self.nonIntInd = nonIntInd
		self.sols = sols 
		self.barY = barY
		self.barYM = barYM
		self.barYb = barYb