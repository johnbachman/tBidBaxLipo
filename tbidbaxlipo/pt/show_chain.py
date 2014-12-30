"""A script that shows the triangle plots of the parameters and the posterior
probability at different steps of the chain for inspection of convergence.
"""
import pickle
import triangle
from matplotlib import pyplot as plt
import sys

def triangle_plots(sampler):
    """Triangle plots of lowest and highest temperature chains."""
    chain = sampler.flatchain
    fig = triangle.corner(chain[0])
    fig.suptitle('Triangle plot, lowest temp')
    #fig.savefig('triangle_low.png')
    fig = triangle.corner(chain[-1])
    fig.suptitle('Triangle plot, highest temp')
    #fig.savefig('triangle_high.png')

def plot_chain_convergence(sampler):
    """Four-panel plot showing convergence of four lowest-temperature chains.

    Useful for determining if the temperature spacing at the lowest part of
    the chain is adequate.
    """
    plt.figure('Chain convergence')
    plt.subplot(2, 2, 1)
    plt.plot(sampler._lnprob[0,:,:].T, alpha=0.1)
    plt.title('0th chain')
    plt.subplot(2, 2, 2)
    plt.plot(sampler._lnprob[1,:,:].T, alpha=0.1)
    plt.title('1st chain')
    plt.subplot(2, 2, 3)
    plt.plot(sampler._lnprob[2,:,:].T, alpha=0.1)
    plt.title('2nd chain')
    plt.subplot(2, 2, 4)
    plt.plot(sampler._lnprob[3,:,:].T, alpha=0.1)
    plt.title('3rd chain')
    #plt.savefig('chain_convergence.png')

if __name__ == '__main__':

    usage_msg =  "Usage:\n"
    usage_msg += " python show_chain.py chain_filename.mcmc\n"

    if len(sys.argv) < 2:
        print usage_msg
        sys.exit()

    chain_filename = sys.argv[1]

    # Unpickle the chain
    with open(chain_filename) as f:
        (gf, sampler) = pickle.load(f)

    # Show plots
    plt.ion()
    triangle_plots(sampler)
    plot_chain_convergence(sampler)

