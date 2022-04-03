import numpy as np
from random import randint
from hermite import ConvertToHNF, ConvertToHNF2
import time

def genMat(m, n):
	M = np.array([[randint(-pow(2, 32), pow(2,32)-1) for x in range(n)] for x in range(m)])
	while np.linalg.matrix_rank(M) != m:
		M = np.array([[randint(-pow(2, 32), pow(2,32)-1) for x in range(n)] for x in range(m)])
	return np.array(M, dtype=object)


if __name__ == "__main__":
	sz = [[5, 5], [5, 10], [10, 10], [10, 15], [15, 15], [15, 25]]
	numTrials = 100

	for m, n in sz:
		avgTime = 0
		maxTime = 0
		avgNumOps = 0
		maxNumOps = 0
		avgMaxElement = 0
		maxMaxElement = 0
		for i in range(numTrials):
			M = genMat(m, n)
			start = time.time() * 1000
			M, ops, maxInt = ConvertToHNF2(M, True)
			curTime = time.time() * 1000 - start

			avgTime += curTime
			maxTime = max(maxTime, curTime)

			avgNumOps += len(ops)
			maxNumOps = max(maxNumOps, len(ops))

			avgMaxElement += maxInt
			maxMaxElement = max(maxMaxElement, maxInt)
		#print (m,n,len(str(avgMaxElement//numTrials)))
		print ("%d x %d: avgTime: %d, maxTime: %d, avgOps: %d, maxOps: %d, avgMax: %d, maxMax: %d" \
			% (m,n,avgTime//numTrials,maxTime,avgNumOps//numTrials,maxNumOps,len(str(avgMaxElement//numTrials)),len(str(maxMaxElement))))
