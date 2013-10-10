from pylab import *
from tbidbaxlipo.util.xkcdify import XKCDify

bax_concs = linspace(0, 1000, 100)

saturating = bax_concs / (bax_concs + 200)
saturating /= max(saturating)
#cooperative = bax_concs**10 / (bax_concs**10 + 500**10)
cooperative = bax_concs**4
cooperative /= max(cooperative)
linear = 0.001 * bax_concs
f = figure()
plot(bax_concs, linear, color='r')
plot(bax_concs, saturating, color='g')
plot(bax_concs, cooperative, color='b')
xlabel('[Bax]')
ylabel('Permeabilization')
title('What is the Bax dose-response?')
text(50, 0.72, 'Saturating?')
text(680, 0.1, 'Cooperative?')
text(400, 0.3, 'Linear?')
#XKCDify(gca())
show()
