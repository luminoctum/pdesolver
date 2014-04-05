#! /usr/bin/env python2.7
from pylab import *
import sys

ua = []

for line in sys.stdin:
    ua.append(map(float, line.split()[2:]))
nx = len(ua[0])
xa = linspace(-5, 5, nx)

figure(figsize = (12, 10))
#plot(xa, ua[0, :], label = 'T = 0')
plot(xa, ua[0], label = 'T = 1')
#plot(xa, ue, '.', label = 'Exact solution')
#ylim(0, 4)
xlim(-5, 5)
legend(loc = 0)
#savefig('./HW3.1.200.png',bbox_inches = 'tight')
show()
