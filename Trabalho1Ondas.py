from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import numpy as np

K = 200 # steps
N = 2000 # steps
dz = 0.1 # m

c = 300000000
l = 1e11
u = 0.9 * c

R = 0
L = 5e-3
G = 0
C = 2e-6

Zo = 50
# Rs = 75
Rs = 75

Rl = "infinito"
# Rl = 0
# Rl = 1e2

if Rl == "infinito":
    RoL = 1
else:
    RoL = ((Rl / Zo) - 1) / ((Rl / Zo) + 1)

RoS = ((Rs / Zo) - 1) / ((Rs / Zo) + 1)


tfim = 10 * l / u

vp = 1 / np.sqrt(L * C)

dt = dz / vp

c1 = - (2 * dt) / (dt * dz * R + 2 * dz * L)
c2 = (2 * L - dt * R) / (2 * L + dt * R)
c3 = - (2 * dt) / (dt * dz * G + 2 * dz * C)
c4 = (2 * C - dt * G) / (2 * C + dt * G)

Vs = np.ones(N) * 2
# Vs = np.array([1 if i <= (l / (100 * u)) else 0 for i in range(N)])

i = np.array([np.zeros(K)] * N)
v = np.array([np.zeros(K)] * N)

x = [i for i in range(N)]
y1 = [0]
y2 = [0]
y3 = [0]
y4 = [0]
xx = [i for i in range(K)]

Vs_t = np.array([i * Zo / (Zo + Rs) for i in Vs])

for n in range(1, N):
    i[n] = c1 * (v[n - 1] - np.concatenate((Vs_t[n - 1:n], v[n - 1][:-1]))) + c2 * i[n - 1]

    if Rl == "infinito":
        v[n] = c3 * (np.concatenate((i[n][1:], [0])) - i[n]) + c4 * v[n - 1]
    elif Rl == 0:
        v[n] = np.concatenate((c3 * (i[n][1:] - i[n][:-1]), [0])) + c4 * v[n - 1]
    else:
        v[n] = c3 * (np.concatenate((i[n][1:], v[n - 1][K - 1:K] / Rl)) - i[n]) + c4 * v[n - 1]

    y1.append(v[n][0])
    y2.append(v[n][K - 2])
    # y2.append(v[n][(K - 1) * 1 // 2])
    y3.append(v[n][K - 1])



f = 1

def func(y):
    global f
    ydata = y
    line.set_data(xdata, ydata)
    plt.title("N = " + str(f))
    f += 1
    return line,

def func_incial():
    global f, Vs
    f = 0
    if Vs[-1] == 0:
        ax.set_ylim(-1, 1)
    else:
        ax.set_ylim(-2, 2)
    ax.set_xlim(-50, K + 50)
    line.set_data(xdata, ydata)
    return line,

fig, ax = plt.subplots()
line, = ax.plot([], [])
ax.grid()
xdata = xx
ydata = []

ani = FuncAnimation(fig, func, frames = v, blit = False,
                    interval = 10, repeat = True, init_func = func_incial)

plt.show()


fig, ax = plt.subplots(3, 2, figsize = (10, 6))
fig.suptitle("Resistencia = " + str(Rl))

plt.subplot(321)
plt.plot(x, y1, label = "Tensão em z = -l no tempo")
plt.legend()

plt.subplot(322)
plt.plot(xx, v[599], label = "Tensão no tempo = 600")
plt.legend()

plt.subplot(323)
plt.plot(x, y2, label = "Tensão em z = 1 no tempo")
plt.legend()

plt.subplot(324)
plt.plot(xx, v[1199], label = "Tensão no tempo = 1200")
plt.legend()

plt.subplot(325)
plt.plot(x, y3, label = "Tensão em z = 0 no tempo")
plt.legend()

plt.subplot(326)
plt.plot(xx, v[-1], label = "Tensão no tempo final")
plt.legend()

plt.show()
