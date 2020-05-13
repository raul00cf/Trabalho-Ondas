from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import numpy as np

escolha1 = int(input("Digite o numero da escolha da Carga:"
                 "\n[1] Rl = infinito"
                 "\n[2] Rl = 0"
                 "\n[3] Rl = 100"
                 "\n"))

Rl = 0

if escolha1 == 1:
    Rl = "infinito"
elif escolha1 == 2:
    Rl = 0
elif escolha1 == 3:
    Rl = 1e2
else:
    print("Escolha Impossível!")
    quit()

escolha2 = int(input("Digite o número da escolha da fonte:"
                 "\n[1] Vs = 2u(t)"
                 "\n[2] Vs = u(t) - u(t - t')"
                 "\n"))

exp = 6

K = 200
N = 10 * K

c = 3e8
l = 10 ** exp
u = 0.9 * c

R = 0
L = 1.85185e-07
G= 0
C = 7.40741e-11

Zo = 50
Rs = 75

if Rl == "infinito":
    RoL = 1
else:
    RoL = ((Rl / Zo) - 1) / ((Rl / Zo) + 1)

RoS = ((Rs / Zo) - 1) / ((Rs / Zo) + 1)

tfim = 10 * l / u

vp = 1 / np.sqrt(L * C)

t = l / u
deltat = 10 * t / N
deltaz = l / K

tt = l / (10 * u)

nn = tt / deltat

Vs = np.zeros(N)

c1 = - (2 * deltat) / (deltat * deltaz * R + 2 * deltaz * L)
c2 = (2 * L - deltat * R) / (2 * L + deltat * R)
c3 = - (2 * deltat) / (deltat * deltaz * G + 2 * deltaz * C)
c4 = (2 * C - deltat * G) / (2 * C + deltat * G)

if escolha2 == 1:
    Vs = np.ones(N) * 2
elif escolha2 == 2:
    Vs = np.array([1 if i <= nn else 0 for i in range(N)])
else:
    print("Escolha Impossível!")
    quit()

print('Dados:')
print('\tN = {} K = {}'.format(N, K))
print('\tR = {} Ohm/m\t L = {:.2e} H/m'.format(R, L))
print('\tG = {} S/m\t C = {:.2e} F/m'.format(G, C))
print('\tZo = {} Rs = {}'.format(Zo, Rs))
print('\tRoS = {:.2f} RoL = {:.2f}'.format(RoS, RoL))
print('\tComprimento da linha = {:.1f} km'.format(l / 1e3))
print('\tVelocidade de propagação = {:.1e} m/s'.format(u))
print('\tTempo de ida = {:.2f} ms'.format(t * 1e3))
print('\tTempo estacionário = {:.2f} ms'.format(tfim * 1e3))
print('\tdt = {:.3f} ms'.format(deltat * 1e3))
print('\tdz = {:.1f} km'.format(deltaz / 1e3))
print('\tt\' = {:.2f} ms'.format(tt * 1e3))
print('\tC1 = {:.2f}'.format(c1))
print('\tC2 = {:.2f}'.format(c2))
print('\tC3 = {:.2f}'.format(c3))
print('\tC4 = {:.2f}'.format(c4))

i = np.array([np.zeros(K)] * N)
v = np.array([np.zeros(K)] * N)

y1 = [0]
y2 = [0]
y3 = [0]
y4 = [0]
xii = [i * 1e3 * deltat for i in range(N)]
xxii = [i * deltaz / 1e3 for i in range(K)]

for n in range(1, N):
    i[n] = np.concatenate(((Vs[n - 1:n] - v[n - 1][0:1]) / Rs, c1 * (v[n - 1][1:] - v[n - 1][:-1]) + c2 * i[n - 1][1:]))

    if Rl == "infinito":
        v[n] = c3 * (np.concatenate((i[n][1:], [0])) - i[n]) + c4 * v[n - 1]
    elif Rl == 0:
        v[n] = np.concatenate((c3 * (i[n][1:] - i[n][:-1]), [0])) + c4 * v[n - 1]
    else:
        v[n] = c3 * (np.concatenate((i[n][1:], v[n - 1][K - 1:K] / Rl)) - i[n]) + c4 * v[n - 1]

    y1.append(v[n][0])
    y3.append(v[n][K - 1])
    y2.append(i[n][0])
    y4.append(i[n][K - 1])

f = 1
fim = Vs[-1]

v_ylim = (0, 0)
i_ylim = (0, 0)

if fim == 0:
    if Rl == "infinito":
        v_ylim = (-0.5, 1.5)
        i_ylim = (-0.02, 0.02)
    elif Rl == 0:
        v_ylim = (-1, 1)
        i_ylim = (-0.01, 0.02)
    else:
        v_ylim = (-0.5, 1)
        i_ylim = (-0.01, 0.02)
else:
    if Rl == "infinito":
        v_ylim = (-0.5, 2.5)
        i_ylim = (-0.02, 0.03)
    elif Rl == 0:
        v_ylim = (-1, 2)
        i_ylim = (-0.01, 0.05)
    else:
        v_ylim = (-0.5, 2)
        i_ylim = (-0.01, 0.04)


def func(y):
    global f, fim
    xdata = xxii
    ydata = y[0]
    ydata2 = y[1]
    line.set_data(xdata, ydata)
    line2.set_data(xdata, ydata2)
    plt.suptitle(
        "$R_L = {} \Omega$ | $V_s(t) = {}$ | $t = {:.2f} ms$".format(Rl if Rl != "infinito" else r"\infty",
                                                                             "2u(t)" if fim != 0 else r"u(t) − u(t-{:.2f})".format(tt * 1e3),
                                                                             f * deltat * 1e3))
    f += 1
    return line, line2


def func_inicial():
    global f, fim
    f = 0

    ax.set_ylim(v_ylim)
    ay.set_ylim(i_ylim)

    ax.set_xlim(-5, (l / 1e3) + 5)
    ay.set_xlim(-5, (l / 1e3) + 5)
    line.set_data(xdata, ydata)
    line2.set_data(xdata, ydata2)

    return line, line2


fig = plt.figure()
ax = fig.add_subplot(211)
ay = fig.add_subplot(212)

plt.suptitle(
        "$R_L = {} \Omega$ | $V_s(t) = {}$".format(Rl if Rl != "infinito" else r"\infty",
                                                                             "2u(t)" if fim != 0 else r"u(t) − u(t-{:.2f})".format(tt)))

fig.subplots_adjust(hspace=0.4)

line, = ax.plot([], [])
line2, = ay.plot([], [])

ax.grid()
ax.set_xlabel('z (km)')
ax.set_ylabel('v(z,t) (V)')

ay.grid()
ay.set_xlabel('z (km)')
ay.set_ylabel('i(z,t) (A)')

xdata = []
ydata = []
ydata2 = []

kk = zip(v, i)

ani = FuncAnimation(fig, func, frames=kk, blit=False,
                    interval=20, repeat=False, init_func=func_inicial)

plt.show()



fig, axs = plt.subplots(4, 2, constrained_layout=True, figsize=(10, 8))
fig.suptitle("$R_L = {} \Omega$ | $V(t) = {}$".format(Rl if Rl != "infinito" else r"\infty",
                                                      "2u(t)" if fim != 0 else r"u(t) − u(t-{:.2f})".format(tt * 1e3)))

axs[0][0].plot(xii, y1)
axs[0][0].set_title('Tensão em z = -{:.0f} km durante o tempo'.format(l / 1e3))
axs[0][0].set_xlabel('t (ms)')
axs[0][0].set_ylabel('v(-{:.0f},t) (V)'.format(l / 1e3))

axs[1][0].plot(xii, y3)
axs[1][0].set_title('Tensão em z = 0 km durante o tempo')
axs[1][0].set_xlabel('t (ms)')
axs[1][0].set_ylabel('v(0,t) (V)')

axs[2][0].plot(xxii, v[100])
axs[2][0].set_title('Tensão ao longo da linha em t = {:.2f} ms'.format(100 * deltat * 1e3))
axs[2][0].set_xlabel('z (km)')
axs[2][0].set_ylabel('v(z,{:.2f}) (V)'.format(100 * deltat * 1e3))

axs[3][0].plot(xxii, v[1900])
axs[3][0].set_title('Tensão ao longo da linha em t = {:.2f} ms'.format(1900 * deltat * 1e3))
axs[3][0].set_xlabel('z (km)')
axs[3][0].set_ylabel('v(z,{:.2f}) (V)'.format(1900 * deltat * 1e3))

for zz in range(4):
    axs[zz][0].set_ylim(v_ylim)

axs[0][1].plot(xii, y2)
axs[0][1].set_title('Corrente em z = -{:.0f} km durante o tempo'.format(l / 1e3))
axs[0][1].set_xlabel('t (ms)')
axs[0][1].set_ylabel('i(-{:.0f},t) (A)'.format(l / 1e3))

axs[1][1].plot(xii, y4)
axs[1][1].set_title('Corrente em z = 0 km durante o tempo')
axs[1][1].set_xlabel('t (ms)')
axs[1][1].set_ylabel('i(0,t) (A)')

axs[2][1].plot(xxii, i[100])
axs[2][1].set_title('Corrente ao longo da linha em t = {:.2f} ms'.format(100 * deltat * 1e3))
axs[2][1].set_xlabel('z (km)')
axs[2][1].set_ylabel('i(z,{:.2f}) (A)'.format(100 * deltat * 1e3))

axs[3][1].plot(xxii, i[1900])
axs[3][1].set_title('Corrente ao longo da linha em t = {:.2f} ms'.format(1900 * deltat * 1e3))
axs[3][1].set_xlabel('z (km)')
axs[3][1].set_ylabel('i(z,{:.2f}) (A)'.format(1900 * deltat * 1e3))

for zz in range(4):
    axs[zz][1].set_ylim(i_ylim)

plt.show()
