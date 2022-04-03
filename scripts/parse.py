import numpy as np
from .models import Rational
from .utils import LCM


class InputError(Exception):
	# error raised when input is improperly formatted
	def __init__(self, msg):
		self.msg = msg
		super().__init__(self.msg)

	def __str__(self):
		return "" if not self.msg else self.msg


def ToRational(s):
	s = s.split('/')
	if len(s) == 1:
		try:
			return Rational(int(s[0]), 1)
		except:
			raise InputError('improperly formatted rational %s' % s[0])
	elif len(s) == 2:
		try:
			return Rational(int(s[0]), int(s[1]))
		except:
			raise InputError('improperly formatted rational %s' % s.join('/'))
	else:
		raise InputError('improperly formatted rational %s' % s.join('/'))


def CheckFullRowRank(M):
	m, n = M.shape
	den = 1
	for i in range(m):
		for j in range(n):
			den = LCM(den, M[i][j].d)

	intM = np.array([[M[i][j].n * den // M[i][j].d for j in range(n)] for i in range(m)])

	return np.linalg.matrix_rank(intM) == m


def ParseInput(m, n, s):
	# takes m and n, the dimensions of the matrix
	# and the m x (n + 1) rational matrix itself
	# where the last column corresponds to b
	mat = np.empty((m, n + 1), dtype=Rational)

	inp = s.split('\n')
	if len(inp) != m:
		raise InputError('input does not contain %i rows' % m)

	for i in range(m):
		row = list(map(ToRational, inp[i].split()))

		if len(row) != n + 1:
			raise InputError('row %i has %i =/= n+1 columns' % (i, len(row)))

		mat[i] = np.array(row)

	if not CheckFullRowRank(mat):
		raise InputError('matrix is not of full row rank')

	return mat