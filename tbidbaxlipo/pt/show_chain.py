import pickle
import triangle
from matplotlib import pyplot as plt

plt.ion()

with open('140318fit_pt.pck') as f:
    (gf, chain) = pickle.load(f)

print chain.shape
fig = triangle.corner(chain[0])
fig.savefig('triangle_low.png')
fig = triangle.corner(chain[-1])
fig.savefig('triangle_high.png')

