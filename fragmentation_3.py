import numpy as np
import matplotlib.pyplot as plt
from numba import jit
from scipy.optimize import minimize,Bounds

seq_size = 2000
hit_number = 200
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

# b = np.linspace(0,seq_size-1,2*n_frag + 1)[1::2]
# b += np.random.normal(0,50,b.shape)

x_init = np.linspace(0,seq_size - 1,n_frag+1)
x_init = np.concatenate(([0],np.random.uniform(0,seq_size-1,n_frag-1),[seq_size-1]))
x = x_init
x_space = np.arange(seq_size,dtype=np.float64)

def tooth(x,x1,x2):
    m = (x1+x2)/2
    h = (x2 - x1)/2
    return (x<=m)*(x-x1)/h + (x>m)*(x2-x)/h

@jit(nopython=True)
def get_F(x,b):
    ##1. Classify bases into 2*fragment classes
    x.sort()
    m = (x[1:] + x[:-1])/2
    bins = np.zeros((x.size + m.size))
    bins[::2] = x
    bins[1::2] = m
    classif = np.digitize(b,bins)-1
    fragment_id = (classif/2).astype(np.int64)
    side_id = (np.ceil(classif/2) - np.floor(classif/2)).astype(np.int64)

    ##2. Calculate differentials
    xi = x[fragment_id]
    xip1 = x[fragment_id+1]

    dydxi = (1-2*side_id)*2*(b-xip1)/(xi-xip1)**2
    dydxip1 = -(1-2*side_id)*2*(b-xi)/(xi-xip1)**2

    dydx = np.zeros_like(x)
    for i, id in enumerate(fragment_id):
        dydx[id] += dydxi[i]
        dydx[id+1] += dydxip1[i]
    dydx[0],dydx[-1] = 0,0

    return dydx

@jit(nopython=True)
def get_Ffrag(x,k):
    Ffrag = np.zeros_like(x)
    Ffrag[1:-1] = -k*(2*x[1:-1] - x[2:]-x[:-2])
    return Ffrag
#
@jit(nopython=True)
def get_Efrag(x,k):
    return -0.5*k*((2*x[1:-1] - x[2:]-x[:-2])**2).sum()
    # return (k*(x-x[2:]-x[:-2] -2*dist_opt)**2).sum()



@jit(nopython=True)
def simulate(x_init,b,dt=0.01,tfin = 10000,k = 5e-3):
    x = x_init.copy()
    tspan = np.arange(0,tfin,dt)
    nt = tspan.size
    x_save = np.ones((nt,x_init.shape[0]))*np.nan
    i = 0
    Fmax = 100
    while (i<nt):
        x.sort()
        F = dt*(get_F(x,b)+get_Ffrag(x,k))
        x = x + F
        x_save[i] = x
        i+=1
        # Fmax = np.abs(F).max()
        # if np.mod(i,1000)==0:
        #     print(Fmax)
    # print("Final av dx/dt = ",get_F(x,b)[1:-1].mean()*dt)
    return x_save

def get_E(y):
    x = np.concatenate(([0],y,[seq_size-1]))
    x.sort()
    E = np.zeros_like(x_space,dtype=np.float64)
    for i in range(x.size-1):
        E[int(x[i]):int(x[i+1])] = tooth(np.arange(int(x[i]),int(x[i+1])),x[i],x[i+1])
    return (E*seq).sum()

x_save = simulate(x_init,b,dt=10,tfin=100000,k=5e-3)

def run_optim(nrun = 100,k = 5e-3):
    x_runs = np.zeros((nrun,n_frag+1))
    E_runs = np.zeros(nrun)
    for i in range(nrun):
        x_init = np.concatenate(([0], np.random.uniform(0, seq_size - 1, n_frag - 1), [seq_size - 1]))
        x_save = simulate(x_init, b, dt=10, tfin=100000,k=k)
        x_runs[i] = x_save[-1]
        try:
            E_runs[i] = get_E(x_runs[i][1:-1]) + get_Efrag(x,k)
        except:
            E_runs[i] = np.nan
    return x_runs,E_runs

x_runs,E_runs = run_optim(nrun = 20,k = 1e-4)
#
#
# plt.plot(x_save)
# plt.show()


# x = np.concatenate(([0], y, [seq_size - 1]))
x = x_runs[np.nanargmax(E_runs)]
E = np.zeros_like(x_space, dtype=np.float64)
for i in range(x.size - 1):
    E[int(x[i]):int(x[i + 1])] = tooth(np.arange(int(x[i]), int(x[i + 1])), x[i], x[i + 1])
import seaborn as sb
fig, ax = plt.subplots()
sb.kdeplot(b,bw_adjust=0.1,ax=ax)
# plt.hist(b,bins = 30)
# plt.scatter(b,np.ones_like(b)*0.5)
ax.plot(E*0.002)
ax.set(xlim = (0,seq_size))
fig.show()




y_init = x_init[1:-1]
y_init.sort()
LB = np.zeros_like(y_init)
UB = np.ones_like(y_init)*(seq_size-1)

def run_optim(n_run):
    y_out = np.zeros((n_run,n_frag-1))
    E_out = np.zeros(n_run)
    for i in range(n_run):
        y_init = np.random.uniform(0, seq_size - 1, n_frag - 1)
        y_init.sort()
        try:
            ysol = minimize(get_E,y_init,method = "trust-constr",bounds = Bounds(LB,UB))
            y_out[i] = ysol.x
            E_out[i] = ysol.fun
        except:
            y_out[i] = np.nan
            E_out[i] = np.nan
    return y_out,E_out

y_out, E_out = run_optim(10)

y = y_out[np.nanargmin(E_out)]
x = np.concatenate(([0], y, [seq_size - 1]))
x.sort()
E = np.zeros_like(x_space, dtype=np.float64)
for i in range(x.size - 1):
    E[int(x[i]):int(x[i + 1])] = tooth(np.arange(int(x[i]), int(x[i + 1])), x[i], x[i + 1])
plt.plot(E)
plt.plot(seq*0.5)
plt.show()

