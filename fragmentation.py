import numpy as np
import matplotlib.pyplot as plt
from numba import jit

seq_size = 200
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
# x_init = np.linspace(0,seq_size-1,n_frag+1)
x = x_init

@jit(nopython=True)
def get_disp(x,b):
    D = np.outer(x,np.ones(b.size)) - np.outer(np.ones(x.size),b)
    return D

@jit(nopython=True)
def get_b_mask(x,b):
    xm1 = x[:-2]
    xp1 = x[2:]
    b_mask = np.zeros((x.size,b.size),dtype=np.bool_)
    b_mask[1:-1] = (b > xm1.reshape(-1,1))*(b<xp1.reshape(-1,1))
    return b_mask

@jit(nopython=True)
def get_E(x,b,w):
    b_mask = get_b_mask(x,b)
    disp = get_disp(x,b)
    E = -np.sum(w*disp**2 * b_mask)
    return E



def get_F(x,b,w):
    b_mask = get_b_mask(x,b)
    disp = get_disp(x,b)
    F = np.sum(w*disp*b_mask,axis=1)
    return F

@jit(nopython=True)
def get_Ffrag(x,k):
    Ffrag = np.zeros_like(x)
    Ffrag[1:-1] = -k*(2*x[1:-1] - x[2:]-x[:-2])
    return Ffrag
#
@jit(nopython=True)
def get_Efrag(x,k):
    return 0.5*k*((2*x[1:-1] - x[2:]-x[:-2])**2).sum()
    # return (k*(x-x[2:]-x[:-2] -2*dist_opt)**2).sum()


def replace_x(x,xmin,xmax,nx):
    mask = ((x<xmin)+(x>xmax))
    if mask.any():
        # x = ~mask*x + np.random.uniform(xmin,xmax,nx)*mask
        # print("yes")
        x = ~mask*x + (xmax - np.random.uniform(0,(xmax-xmin)*0.1,nx))*(x>xmax) + (x<xmin)*np.random.uniform(0,(xmax-xmin)*0.1,nx)
    return x

def simulate(x_init,b,w,k = 10,dt = 0.001,tfin = 100):
    tspan = np.arange(0,tfin,dt)
    xmin,xmax,nx = x_init.min(),x_init.max(),x_init.size
    x = x_init.copy()
    x_save = np.zeros((tspan.size,x.size))
    E_save = np.zeros((tspan.size))
    for i in range(tspan.size):
        F = get_F(x,b,w)
        x[1:-1] = x[1:-1][x[1:-1].argsort()]
        Ffrag = get_Ffrag(x,k)
        x+= dt*(F+Ffrag)
        x = replace_x(x,xmin,xmax,nx)
        x_save[i] = x
        E_save[i] = get_E(x,b,w) + get_Efrag(x,k)

    return x_save,E_save

dt = 0.001
tfin = 20
x_save,E_save = simulate(x_init,b,w*0.13,1,dt = dt,tfin = tfin)

plt.plot(x_save,color="black")
x_save,E_save = simulate(x_init,b,w*0,1,dt = dt,tfin = tfin)

plt.plot(x_save,color="red")

plt.show()

fig, ax = plt.subplots()
ax.plot(np.arange(0,tfin,dt),E_save)
# ax.set(xscale="log")
fig.show()

fig, ax = plt.subplots()
ax.plot(seq)
for xi in x_save[-1]:
    ax.plot((xi,xi),(-0.3,1.3),linestyle="--",color="grey")
fig.show()

n_sample = 10000
x_sample = np.random.uniform(0,seq_size-1,(n_sample,4))
x_sample = np.column_stack((np.zeros(n_sample),x_sample,np.ones(n_sample)*(seq_size-1)))


@jit(nopython=True)
def get_fragment_sizes(x):
    x_sorted = x[x.argsort()]
    return x_sorted[1:] -x_sorted[:-1]

@jit(nopython=True)
def get_fragment_sizes_sample(x_sample):
    fragments = np.zeros((x_sample.shape[0],x_sample.shape[1]-1))
    for i, x in enumerate(x_sample):
        fragments[i] = get_fragment_sizes(x)
    return fragments

# fragment_sizes_sample = get_fragment_sizes_sample(x_sample)
# opt_frag_size = seq_size/n_frag
# dev = 0.3
# fragsize_lim = (opt_frag_size*(1-dev),opt_frag_size*(1+dev))
# x_sample = x_sample[(fragment_sizes_sample.min(axis=1)>fragsize_lim[0])*(fragment_sizes_sample.max(axis=1)<fragsize_lim[1])]

# @jit(nopython=True)
def get_E_sample(x_sample):
    E_sample,E_frag_sample = np.zeros(x_sample.shape[0]),np.zeros(x_sample.shape[0])
    for i, x in enumerate(x_sample):
        E_sample[i], E_frag_sample[i] = get_E(x,b,w),get_Efrag(x,1)
    return E_sample,E_frag_sample

E_sample,E_frag_sample = get_E_sample(x_sample)

fig, ax = plt.subplots()
ax.scatter(E_sample,E_frag_sample,s=2)
# ax.set(xscale="log",yscale="log")
fig.show()

k = 10


x_best = x_sample[np.argmin(E_sample + E_frag_sample*k)]
fig, ax = plt.subplots()
ax.plot(seq)
for xi in x_best:
    ax.plot((xi,xi),(-0.3,1.3),linestyle="--",color="grey")
fig.show()