# Bibliotecas utilizadas

from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import numpy as np

#######################################################################################################################

# Parametros determinados pelo PDF do trabalho:
# apresenta R = G = 0, uma vez que a linha não apresenta perdas, e L = 0,185185 micro H e C = 74,0741 pico F, uma vez
# que determinam a impedância caracteristica Zo = 50 Ohms
# Determina 200 pontos para a determinação dos passos em Z e 2000 pontos para a determinação dos passos em t
# de acordo com o método FDTD
# Determina a resistencia interna Rs = 75 Ohm

K = 200
N = 2000

c = 3e8
l = 10 ** 6
u = 0.9 * c

R = 0
L = 1.85185e-07
G = 0
C = 7.40741e-11

Zo = (L / C) ** 0.5
Rs = 75

#######################################################################################################################

# Recebe a escolha das cargas a serem consideradas, sendo Rl infinito, 0 ou 100 Ohm
# Com os valores, calcula alguns parametros necessários, sendo os valores dos coeficientes de reflexão das cargas
# Recebe também a escolha da tensão possíveis, sendo Vs = 2u(t) e Vs = u(t)-u(t-t'), sendo t' o tempo de um décimo
# do tempo de propagação da onda, determinado no PDF

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

c1 = - (2 * deltat) / (deltat * deltaz * R + 2 * deltaz * L)
c2 = (2 * L - deltat * R) / (2 * L + deltat * R)
c3 = - (2 * deltat) / (deltat * deltaz * G + 2 * deltaz * C)
c4 = (2 * C - deltat * G) / (2 * C + deltat * G)

Vs = np.zeros(N)

if escolha2 == 1:
    Vs = np.ones(N) * 2
elif escolha2 == 2:
    Vs = np.array([1 if i <= nn else 0 for i in range(N)])
else:
    print("Escolha Impossível!")
    quit()

#######################################################################################################################

# Apresenta dados importantes para o usuário

print('Dados:')
print('\tN = {} K = {}'.format(N, K))
print('\tR = {} Ohm/m\t L = {:.2e} H/m'.format(R, L))
print('\tG = {} S/m\t C = {:.2e} F/m'.format(G, C))
print('\tZo = {:.2f} Rs = {}'.format(Zo, Rs))
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

#######################################################################################################################

# Cria matriz de corrente e tensão que será utilizada no método FDTD e determina como tensão inicial em z = -l como a
# tensão do gerador após a queda de tensão da resistencia interna
# Cria 4 vetores que serão salvos os dados para os eixos y e 2 vetores para os eixos x nos gráficos

i = np.array([np.zeros(K)] * N)
v = np.array([np.zeros(K)] * N)
v[0][0] = Vs[0] * (Zo / (Rs + Zo))

y1 = [0]
y2 = [0]
y3 = [0]
y4 = [0]
xii = [i * 1e3 * deltat for i in range(N)]
xxii = [i * deltaz / 1e3 for i in range(K)]

#######################################################################################################################

# Realiza o método FDTD, sendo em cada iteração de n calcula a tensão e a corrente ao longo da linha de acordo com as
# equações discretas do método FDTD
# Para diminuir o ruido, quando há uma variação na tensão do gerador, determina a tensão em z = -l como a tensão do
# gerador após a queda de tensão em Rs
# Quando Rl for 0, ou seja, apresentar um curto circuito na carga Rl, a tensão em z = 0 é fixada como 0
# Quando Rl for infinito, ou seja, apresentar um ramo aberto na carga Rl, a corrente em z = 0 é fixada como 0
# Adiciona, após a iteração, a tensão em z = -l e z = 0 e a corrente em z = -l e z = 0

for n in range(1, N):
    pular = False
    if (Vs[n] != Vs[n - 1]):
        pular = True

    i[n] = np.concatenate(((Vs[n - 1:n] - v[n - 1][0:1]) / Rs, c1 * (v[n - 1][1:] - v[n - 1][:-1]) + c2 * i[n - 1][1:]))

    if Rl == "infinito":
        v[n] = c3 * (np.concatenate((i[n][1:], [0])) - i[n]) + c4 * v[n - 1]
    elif Rl == 0:
        v[n] = np.concatenate((c3 * (i[n][1:] - i[n][:-1]) + c4 * v[n - 1][:-1], [0]))
    else:
        v[n] = c3 * (np.concatenate((i[n][1:], v[n - 1][K - 1:K] / Rl)) - i[n]) + c4 * v[n - 1]

    if pular:
        v[n][0] = Vs[n] * (Zo / (Rs + Zo))

    y1.append(v[n][0])
    y3.append(v[n][K - 1])
    y2.append(i[n][0])
    y4.append(i[n][K - 1])

#######################################################################################################################

# Determina os parâmetros para os gráficos e a animação

f = 1
fim = Vs[-1]

v_ylim = (0, 0)
i_ylim = (0, 0)

if fim == 0:
    if Rl == "infinito":
        v_ylim = (-0.5, 1.4)
        i_ylim = (-0.02, 0.02)
    elif Rl == 0:
        v_ylim = (-0.9, 1)
        i_ylim = (-0.01, 0.03)
    else:
        v_ylim = (-0.3, 1)
        i_ylim = (-0.01, 0.02)
else:
    if Rl == "infinito":
        v_ylim = (-0.5, 2.7)
        i_ylim = (-0.01, 0.03)
    elif Rl == 0:
        v_ylim = (-0.5, 1.5)
        i_ylim = (-0.01, 0.04)
    else:
        v_ylim = (-0.5, 2)
        i_ylim = (-0.01, 0.03)

#######################################################################################################################

# Cria e apresenta 8 gráficos, sendo 4 deles sobre a tensão:
# Tensão em z = -l ao longo do tempo
# Tensão em z = 0 ao longo do tempo
# Tensão ao longo da linha no tempo de meia propagação da onda
# Tensão ao longo da linha no tempo de 9 propagações e meia da onda
# Os outros 4 gráficos são da corrente, sendo os mesmos parametros dos gráficos de tensão

fig, axs = plt.subplots(4, 2, figsize=(10, 8))

fig.subplots_adjust(hspace=0.7, wspace=0.25)

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

#######################################################################################################################

# Cria e apresenta uma animação da tensão e corrente ao longo da linha, variando o tempo, de t = 0 ao tempo de 10
# propagações da onda, no tempo estacionário determinado no PDF

def func(y):
    global f, fim
    xdata = xxii
    ydata = y[0]
    ydata2 = y[1]
    line.set_data(xdata, ydata)
    line2.set_data(xdata, ydata2)
    title.set_text("$t = {:.2f} ms$".format(f * deltat * 1e3))
    f += 1
    return line, line2, title


def func_inicial():
    global f, fim
    f = 0

    xdata = []
    ydata = []
    ydata2 = []

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

plt.suptitle("$R_L = {} \Omega$ | $V_s(t) = {}$".format(Rl if Rl != "infinito" else r"\infty",
                                                        "2u(t)" if fim != 0 else r"u(t) − u(t-{:.2f})".format(
                                                            tt * 1e3)), ha="center")

title = ax.text(0.5, 0.85, "", transform=ax.transAxes, ha="center", fontsize=11)

fig.subplots_adjust(hspace=0.4, left=0.15)

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

ani = FuncAnimation(fig, func, frames=kk, blit=True,
                    interval=50, repeat=False, init_func=func_inicial)

plt.show()

#######################################################################################################################
