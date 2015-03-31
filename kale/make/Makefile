CODEDIR := ../../tbidbaxlipo
PLOTSDIR := $(CODEDIR)/plots
DATADIR := $(CODEDIR)/data
NBD_DATA := $(DATADIR)/Compiled\ Release\ percentages\ and\ NBD\ F-F0.xlsx \
            $(DATADIR)/parse_bid_bim_nbd_release.py
NBD_FRET_DATA := $(DATADIR)/2015-02-27\ -\ Bid-Bim\ FRET\ with\ Bax\ NBD\ mutants\ in\ Tb-DPA\ Liposomes\ compiled\ reps+WT\ Bax.xlsx \
                 $(DATADIR)/parse_bid_bim_fret_nbd_release.py
PQUEUE := sorger_par_unlimited

all: nbd_release fret_nbd_release

clean:
	rm -f ./*.pdf ./*.png

deploy:
	cp ./*.pdf ~/Dropbox/Bax_insertion_paper/figure_panels_JAB

# This bit here allows us to only rebuild when the hash of a file changes.
# See blog post describing the approach here:
# http://blog.jgc.org/2006/04/rebuilding-when-hash-has-changed-not.html
to-md5 = $(patsubst %,%.md5,$1)
from-md5 = $(patsubst %.md5,%,$1)

# The .md5 file is updated only if the hash has changed
%.md5: FORCE
	@$(if $(filter-out $(shell cat $@ 2>/dev/null), $(shell md5sum $*)),md5sum $* > $@)

# Dummy target so the .md5 step is always run
FORCE:

.PHONY: clean deploy nbd_release fret_nbd_release

# Temporary target for testing of MCMC fitting of NBD data
test.mcmc:
	bsub -q $(PQUEUE) -n 20 -W 12:00 -a openmpi mpirun.lsf python -m tbidbaxlipo.plots.mcmc_54_script

nbd_release: nbd_release_endpts nbd_release_initial_rates
nbd_release_endpts: data1_release_endpt_norm.pdf data1_release_endpt_no_norm.pdf data1_nbd_endpt.pdf
nbd_release_initial_rates: data1_init_nbd_bar.pdf

fret_nbd_release: fret_nbd_release_endpts fret_nbd_release_initial_rates
fret_nbd_release_endpts: data2_release_endpt_norm.pdf data2_release_endpt_no_norm.pdf data2_nbd_endpt.pdf
fret_nbd_release_initial_rates: data2_init_nbd_bar.pdf

# --- Normalized release endpoint data for release, dataset 1 ---
data1_release_endpt_norm.pdf: \
        $(PLOTSDIR)/nbd_bax_analysis.py \
        $(PLOTSDIR)/bid_bim_nbd_release/plots.py \
        $(NBD_DATA)
	python -m tbidbaxlipo.plots.bid_bim_nbd_release.plots release_endpoint_norm $(basename $@)

# --- Raw release endpoint data for release, dataset 1 ---
data1_release_endpt_no_norm.pdf: \
        $(PLOTSDIR)/nbd_bax_analysis.py \
        $(PLOTSDIR)/bid_bim_nbd_release/plots.py \
        $(NBD_DATA)
	python -m tbidbaxlipo.plots.bid_bim_nbd_release.plots release_endpoint_no_norm $(basename $@)

# --- NBD endpoints, dataset 1 ---
data1_nbd_endpt.pdf: \
        $(PLOTSDIR)/nbd_bax_analysis.py \
        $(PLOTSDIR)/bid_bim_nbd_release/plots.py \
        $(NBD_DATA)
	python -m tbidbaxlipo.plots.bid_bim_nbd_release.plots nbd_endpoint $(basename $@)

# NBD/release, initial rates
# Building this target creates four files:
#     data1_init_nbd_bar.pdf,
#     data1_init_nbd_bar.pdf,
#     data1_init_scatter_norm_Bid.pdf,
#     data1_init_scatter_norm_Bim.pdf
data1_init_nbd_bar.pdf: \
        $(PLOTSDIR)/nbd_bax_analysis.py \
        $(PLOTSDIR)/bid_bim_nbd_release/plots.py \
        $(NBD_DATA)
	python -m tbidbaxlipo.plots.bid_bim_nbd_release.plots initial_rate_samples data1

# --- Normalized release endpoint data for release, dataset 2 ---
data2_release_endpt_norm.pdf: \
        $(PLOTSDIR)/nbd_bax_analysis.py \
        $(PLOTSDIR)/bid_bim_fret_nbd_release/plots.py \
        $(NBD_FRET_DATA)
	python -m tbidbaxlipo.plots.bid_bim_fret_nbd_release.plots release_endpoint_norm $(basename $@)

# --- Raw release endpoint data for release, dataset 2 ---
data2_release_endpt_no_norm.pdf: \
        $(PLOTSDIR)/nbd_bax_analysis.py \
        $(PLOTSDIR)/bid_bim_fret_nbd_release/plots.py \
        $(NBD_FRET_DATA)
	python -m tbidbaxlipo.plots.bid_bim_fret_nbd_release.plots release_endpoint_no_norm $(basename $@)

# --- NBD endpoints, dataset 2 ---
data2_nbd_endpt.pdf: \
        $(PLOTSDIR)/nbd_bax_analysis.py \
        $(PLOTSDIR)/bid_bim_fret_nbd_release/plots.py \
        $(NBD_FRET_DATA)
	python -m tbidbaxlipo.plots.bid_bim_fret_nbd_release.plots nbd_endpoint $(basename $@)

# FRET/NBD/release, initial rates
# Building this target creates four files:
#     data2_init_nbd_bar.pdf,
#     data2_init_nbd_bar.pdf,
#     data2_init_scatter_norm_Bid.pdf,
#     data2_init_scatter_norm_Bim.pdf
data2_init_nbd_bar.pdf: \
        $(PLOTSDIR)/nbd_bax_analysis.py \
        $(PLOTSDIR)/bid_bim_fret_nbd_release/plots.py \
        $(NBD_FRET_DATA)
	python -m tbidbaxlipo.plots.bid_bim_fret_nbd_release.plots initial_rate_samples data2



