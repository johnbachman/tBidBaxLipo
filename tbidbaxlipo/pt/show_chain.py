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

