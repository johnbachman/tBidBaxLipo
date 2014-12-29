import pickle
import triangle
from matplotlib import pyplot as plt

plt.ion()

with open('pt_140318_nbd_2_conf_1.pck') as f:
    (gf, sample) = pickle.load(f)

chain = sample.flatchain
print chain.shape
fig = triangle.corner(chain[0])
fig.savefig('triangle_low.png')
fig = triangle.corner(chain[-1])
fig.savefig('triangle_high.png')

plt.figure('Chain convergence')
plt.subplot(2, 2, 1)
plt.plot(sample._lnprob[0,:,:].T, alpha=0.1)
plt.title('0th chain')
plt.subplot(2, 2, 2)
plt.plot(sample._lnprob[1,:,:].T, alpha=0.1)
plt.title('1th chain')
plt.subplot(2, 2, 3)
plt.plot(sample._lnprob[2,:,:].T, alpha=0.1)
plt.title('2th chain')
plt.subplot(2, 2, 4)
plt.plot(sample._lnprob[3,:,:].T, alpha=0.1)
plt.title('3th chain')
plt.savefig('chain_convergence.png')
