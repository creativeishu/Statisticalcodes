import numpy
import pylab

Principal = 1e5
MaxTime = 20
rate = 0.1

Principal += Principal * MaxTime * rate

def remaining(paid):
	global Principal, MaxTime, rate
	Principal -= paid
	return Principal

	