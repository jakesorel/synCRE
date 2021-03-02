import numpy as np
import matplotlib.pyplot as plt
from numba import jit
from scipy.optimize import minimize,Bounds

seq_size = 2000
hit_number = 20
weights = 1


def build_seq(seq_size,hit_number):
    seq = np.zeros(seq_size)
    for i in range(hit_number):
        seq[int(np.random.random()*seq_size)] = 1
    return seq

seq = build_seq(seq_size,hit_number)
b = np.nonzero(seq)[0]
w = np.ones_like(b)

n_frag = 4
x_init = np.concatenate(([0],np.random.uniform(0,seq_size-1,n_frag-1),[seq_size-1]))
x = x_init
x_space = np.arange(seq_size,dtype=np.float64)


# x.sort()
@jit(nopython=True)
def quad(x,x1,x2):
    # return (x2-x1) - 4*(x - (x1+x2)/2)**2/(x2-x1)
    return (x2-x1)**2 - 4*(x - (x1+x2)/2)**2


def get_E(y):
    x = np.concatenate(([0],y,[seq_size-1]))
    x.sort()
    E = np.zeros_like(x_space,dtype=np.float64)
    for i in range(x.size-1):
        E[int(x[i]):int(x[i+1])] = quad(np.arange(int(x[i]),int(x[i+1])),x[i],x[i+1])
    return E.sum()

y_init = x_init[1:-1]
y_init.sort()
LB = np.zeros_like(y_init)
UB = np.ones_like(y_init)*(seq_size-1)

ysol = minimize(get_E,y_init,method = "trust-constr",bounds = Bounds(LB,UB))

y = ysol.x
x = np.concatenate(([0], y, [seq_size - 1]))
x.sort()
E = np.zeros_like(x_space, dtype=np.float64)
for i in range(x.size - 1):
    E[int(x[i]):int(x[i + 1])] = quad(np.arange(int(x[i]), int(x[i + 1])), x[i], x[i + 1])
plt.plot(E)
plt.show()

n_trial = 1000
y_trial = np.random.uniform(0,seq_size-1,(n_trial,3))

E_trial = [get_E(y) for y in y_trial]

y = y_trial[np.where(E == min(E))[0][0]]
x = np.concatenate(([0], y, [seq_size - 1]))
x.sort()
E = np.zeros_like(x_space, dtype=np.float64)
for i in range(x.size - 1):
    E[int(x[i]):int(x[i + 1])] = quad(np.arange(int(x[i]), int(x[i + 1])), x[i], x[i + 1])
plt.plot(E)
plt.show()