CODEDIR := ../../tbidbaxlipo
PLOTSDIR := $(CODEDIR)/plots
DATADIR := $(CODEDIR)/data
NBD_MCMC_DIR := $(PLOTSDIR)/bid_bim_nbd_release/mcmc
NBD_FRET_MCMC_DIR := $(PLOTSDIR)/bid_bim_fret_nbd_release/mcmc
NBD_FRET_DUAL_MCMC_DIR := $(PLOTSDIR)/bid_bim_fret_nbd_release/fret_mcmc
BAX_NBD_FRET_DUAL_MCMC_DIR := $(PLOTSDIR)/bax_bax_fret_nbd_release/fret_mcmc
NBD_DATA := $(DATADIR)/Compiled\ Release\ percentages\ and\ NBD\ F-F0.xlsx \
            $(DATADIR)/parse_bid_bim_nbd_release.py
NBD_FRET_DATA := $(DATADIR)/2015-02-27\ -\ Bid-Bim\ FRET\ with\ Bax\ NBD\ mutants\ in\ Tb-DPA\ Liposomes\ compiled\ reps+WT\ Bax.xlsx \
                 $(DATADIR)/parse_bid_bim_fret_nbd_release.py
BAX_FRET_DATA := $(DATADIR)/2016\ -\ Bax-Bax\ FRET\ dataset.xlsx \
                 $(DATADIR)/parse_bax_bax_fret_nbd_release.py
x140430 := $(CODEDIR)/plots/x140430_Bid_Bax_FRET/mcmc
#PQUEUE := sorger_par_unlimited
PQUEUE := parallel
# md5 sum program to use (md5 on mac, md5sum on Linux)
MD5 := $(MD5_NAME)

all: nbd_release fret_nbd_release bax_bax_fret nbd_labeling_ratios_table

deploy:
	rsync -av ./ ~/Dropbox/Bachman-Sorger\ Talks\ and\ Papers/Kale-Bachman\ Bax\ mutants/figures/panels

clean:
	rm -f ./*.pdf ./*.png ./*.csv

clean_pt_data1:
	cd $(NBD_MCMC_DIR); rm -f *.md5 *.fit *.deps.txt *.ensemble

clean_pt_data2:
	cd $(NBD_FRET_MCMC_DIR); rm -f *.md5 *.fit *.deps.txt *.ensemble

clean_pt_data3_fret:
	cd $(BAX_NBD_FRET_MCMC_DIR); rm -f *.md5 *.fit *.deps.txt *.ensemble

# This bit here allows us to only rebuild when the hash of a file changes.
# See blog post describing the approach here:
# http://blog.jgc.org/2006/04/rebuilding-when-hash-has-changed-not.html
to-md5 = $(patsubst %,%.md5,$1)
from-md5 = $(patsubst %.md5,%,$1)

# The .md5 file is updated only if the hash has changed
%.md5: FORCE
	@$(if $(filter-out $(shell cat $@ 2>/dev/null), $(shell $(MD5) $*)),$(MD5) $* > $@)

# Dummy target so the .md5 step is always run
FORCE:

.PHONY: clean clean_pt deploy nbd_release fret_nbd_release bax_bax_fret

# Running this script generates both the dependency list and the .yaml files
# specifying the fit parameters for each individual model
%.deps.txt: $(call to-md5, %.fit.ensemble)
	python -m tbidbaxlipo.pt.generate_model_ensemble_fit_files $(call from-md5, $<)

# In this case, we know that the .fit file exists, since it has been
# regenerated in the step that regenerated the dependency file. What
# we care about is whether the hash is any different.
%.mcmc: $(call to-md5, %.fit)
	#python -m tbidbaxlipo.pt.compile_models $(call from-md5, $<)
	#bsub -q $(PQUEUE) -n 10 -W 100:00 -N -o $(call from-md5, $<).out -a openmpi mpirun.lsf python -m tbidbaxlipo.pt.run_pt $(call from-md5, $<) 1 $(call from-md5, $<).pos
	#qsub -b y -cwd -V -o $(call from-md5, $<).out -e $(call from-md5, $<).err -pe orte 16 mpirun python $(CODEDIR)/pt/run_pt.py $(call from-md5, $<) 1 $(call from-md5, $<).pos
	qsub -b y -cwd -V -o $(call from-md5, $<).out -e $(call from-md5, $<).err python $(CODEDIR)/pt/run_pt.py $(call from-md5, $<) 1 $(call from-md5, $<).pos nompi

# Don't delete the .md5 files! Without this rule they are treated as
# intermediates and deleted.
.PRECIOUS: %.md5

# -- MCMC Fitting for 140430 Bid/Bax FRET pilot data -------------------
pt_140430_36C: $(x140430)/pt_140430_36C.deps.txt
-include $(x140430)/pt_140430_36C.deps.txt

pt_140430_68C: $(x140430)/pt_140430_68C.deps.txt
-include $(x140430)/pt_140430_68C.deps.txt

pt_140430_126C: $(x140430)/pt_140430_126C.deps.txt
-include $(x140430)/pt_140430_126C.deps.txt

# -- MCMC Fitting for Dataset 1, complete set of NBD-Bax trajectories --
pt_data1: pt_data1_deps

# The pt_data1_deps target allows for the dependencies list to be created
# without running all of the jobs that pt_data1 itself depends on. For example,
# > make pt_data1_deps
# then
# > make pt_data1_Bid_NBD_126_r1
# or
# > make pt_data1_Bid_NBD_126_r1_2confs.mcmc
# ...allows selected reps/models to be run explicitly.
pt_data1_deps: $(NBD_MCMC_DIR)/pt_data1.deps.txt
-include $(NBD_MCMC_DIR)/pt_data1.deps.txt

$(NBD_MCMC_DIR)/pt_data1.deps.txt: \
    $(NBD_MCMC_DIR)/generate_fit_files_for_nbd_curves.py
	python $(NBD_MCMC_DIR)/generate_fit_files_for_nbd_curves.py $(NBD_MCMC_DIR)

# -- MCMC Fitting for Dataset 2, complete set of NBD-Bax trajectories with FRET --
pt_data2: pt_data2_deps

pt_data2_deps: $(NBD_FRET_MCMC_DIR)/pt_data2.deps.txt
-include $(NBD_FRET_MCMC_DIR)/pt_data2.deps.txt

$(NBD_FRET_MCMC_DIR)/pt_data2.deps.txt: \
    $(NBD_FRET_MCMC_DIR)/generate_fit_files_for_nbd_curves.py
	python $(NBD_FRET_MCMC_DIR)/generate_fit_files_for_nbd_curves.py $(NBD_FRET_MCMC_DIR)

# -- MCMC fitting for NBD + FRET, dataset 2
pt_data2_fret: pt_data2_fret_deps

pt_data2_fret_deps: $(NBD_FRET_DUAL_MCMC_DIR)/pt_data2_fret.deps.txt
-include $(NBD_FRET_DUAL_MCMC_DIR)/pt_data2_fret.deps.txt

$(NBD_FRET_DUAL_MCMC_DIR)/pt_data2_fret.deps.txt: \
    $(NBD_FRET_DUAL_MCMC_DIR)/generate_fit_files_for_nbd_curves.py
	python $(NBD_FRET_DUAL_MCMC_DIR)/generate_fit_files_for_nbd_curves.py $(NBD_FRET_DUAL_MCMC_DIR)


# -- MCMC fitting for NBD + FRET, dataset 3
pt_data3_fret: pt_data3_fret_deps

pt_data3_fret_deps: $(BAX_NBD_FRET_DUAL_MCMC_DIR)/pt_data3_fret.deps.txt
-include $(BAX_NBD_FRET_DUAL_MCMC_DIR)/pt_data3_fret.deps.txt

$(BAX_NBD_FRET_DUAL_MCMC_DIR)/pt_data3_fret.deps.txt: \
    $(BAX_NBD_FRET_DUAL_MCMC_DIR)/generate_fit_files_for_nbd_curves.py
	python $(BAX_NBD_FRET_DUAL_MCMC_DIR)/generate_fit_files_for_nbd_curves.py $(BAX_NBD_FRET_DUAL_MCMC_DIR)


# -- Dataset 1, NBD + Tb release, summary plots ------------------------
nbd_release: nbd_release_endpts nbd_release_initial_rates nbd_example_curves \
             data1_release_peak_times.csv data1_derivatives_Bid_15_r1.pdf \
             data1_evidence_2confs.pdf
nbd_release_endpts: \
    data1_release_endpt_norm.pdf data1_release_endpt_no_norm.pdf \
    data1_nbd_endpt_norm.pdf data1_nbd_endpt_no_norm.pdf
nbd_release_initial_rates: data1_init_nbd_bar.pdf
nbd_example_curves: data1_nbd_example_curves.pdf
nbd_param_densities: data1_Bid_k1k2_densities.pdf data1_Bid_c1c2_densities.pdf\
                     data1_Bim_k1k2_densities.pdf data1_Bim_c1c2_densities.pdf

# -- Dataset 2, Bid FRET + NBD + Tb release, summary plots -------------
fret_nbd_release: fret_nbd_release_endpts fret_nbd_release_initial_rates \
    data2_fret_endpt.pdf data2_example_curves.pdf data2_delta_fret.pdf \
    data2_example_nbd_fret_fits.pdf
fret_nbd_release_endpts: \
    data2_release_endpt_norm.pdf data2_release_endpt_no_norm.pdf \
    data2_nbd_endpt_norm.pdf data2_nbd_endpt_no_norm.pdf
fret_nbd_release_initial_rates: data2_init_nbd_bar.pdf

# -- Bax-Bax FRET data (141203) ------------------------------------------
bax_bax_fret: 141203_example_nbd_fret_fits.pdf

# -- Dataset 3, Bax-Bax FRET + NBD + Tb release, summary plots -------------
bax_fret_nbd_release: bax_fret_nbd_release_endpts bax_fret_nbd_release_initial_rates \
    data3_fret_endpt.pdf data3_example_curves.pdf data3_delta_fret.pdf \
    data3_example_nbd_fret_fits.pdf
bax_fret_nbd_release_endpts: \
    data3_release_endpt_norm.pdf data3_release_endpt_no_norm.pdf \
    data3_nbd_endpt_norm.pdf data3_nbd_endpt_no_norm.pdf
bax_fret_nbd_release_initial_rates: data2_init_nbd_bar.pdf


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

# --- NBD endpoints, dataset 1, normalized to F/F0 ---
data1_nbd_endpt_norm.pdf: \
        $(PLOTSDIR)/nbd_bax_analysis.py \
        $(PLOTSDIR)/bid_bim_nbd_release/plots.py \
        $(NBD_DATA)
	python -m tbidbaxlipo.plots.bid_bim_nbd_release.plots nbd_endpoint_norm $(basename $@)

# --- NBD endpoints, dataset 1, no normalization (raw values) ---
data1_nbd_endpt_no_norm.pdf: \
        $(PLOTSDIR)/nbd_bax_analysis.py \
        $(PLOTSDIR)/bid_bim_nbd_release/plots.py \
        $(NBD_DATA)
	python -m tbidbaxlipo.plots.bid_bim_nbd_release.plots nbd_endpoint_no_norm $(basename $@)

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

data1_nbd_example_curves.pdf: \
        $(PLOTSDIR)/bid_bim_nbd_release/plot_example_curves.py \
        $(NBD_DATA)
	python -m tbidbaxlipo.plots.bid_bim_nbd_release.plot_example_curves

data1_release_peak_times.csv: \
        $(PLOTSDIR)/nbd_bax_analysis.py \
        $(PLOTSDIR)/bid_bim_nbd_release/plots.py \
        $(NBD_DATA)
	python -m tbidbaxlipo.plots.bid_bim_nbd_release.plots calc_release_peaks

data1_derivatives_Bid_15_r1.pdf: \
        $(PLOTSDIR)/nbd_bax_analysis.py \
        $(PLOTSDIR)/bid_bim_nbd_release/plots.py \
        $(NBD_DATA)
	python -m tbidbaxlipo.plots.bid_bim_nbd_release.plots example_derivatives

data1_evidence_2confs.pdf: \
        $(PLOTSDIR)/nbd_bax_analysis.py \
        $(PLOTSDIR)/bid_bim_nbd_release/plots.py \
        $(PLOTSDIR)/bid_bim_nbd_release/plot_bf_values.py \
        $(NBD_DATA)
	python -m tbidbaxlipo.plots.bid_bim_nbd_release.plots evidence

data1_Bid_k1k2_densities.pdf: \
        $(PLOTSDIR)/bid_bim_nbd_release/plot_k1_k2_dists.py \
        $(PLOTSDIR)/bid_bim_nbd_release/process_3conf_mcmc.py \
        $(NBD_MCMC_DIR)/Bid_k1_density_mx.txt \
        $(NBD_MCMC_DIR)/Bid_k2_density_mx.txt
	python -m tbidbaxlipo.plots.bid_bim_nbd_release.plot_k1_k2_dists plot k1k2 \
		$(NBD_MCMC_DIR)/Bid_k1_density_mx.txt \
		$(NBD_MCMC_DIR)/Bid_k2_density_mx.txt \
		data1_Bid_k1k2_densities

data1_Bid_c1c2_densities.pdf: \
        $(PLOTSDIR)/bid_bim_nbd_release/plot_k1_k2_dists.py \
        $(PLOTSDIR)/bid_bim_nbd_release/process_3conf_mcmc.py \
        $(NBD_MCMC_DIR)/Bid_c1_density_mx.txt \
        $(NBD_MCMC_DIR)/Bid_c2_density_mx.txt
	python -m tbidbaxlipo.plots.bid_bim_nbd_release.plot_k1_k2_dists plot c1c2 \
		$(NBD_MCMC_DIR)/Bid_c1_density_mx.txt \
		$(NBD_MCMC_DIR)/Bid_c2_density_mx.txt \
		data1_Bid_c1c2_densities

data1_Bim_k1k2_densities.pdf: \
        $(PLOTSDIR)/bid_bim_nbd_release/plot_k1_k2_dists.py \
        $(PLOTSDIR)/bid_bim_nbd_release/process_3conf_mcmc.py \
        $(NBD_MCMC_DIR)/Bim_k1_density_mx.txt \
        $(NBD_MCMC_DIR)/Bim_k2_density_mx.txt
	python -m tbidbaxlipo.plots.bid_bim_nbd_release.plot_k1_k2_dists plot k1k2 \
		$(NBD_MCMC_DIR)/Bim_k1_density_mx.txt \
		$(NBD_MCMC_DIR)/Bim_k2_density_mx.txt \
		data1_Bim_k1k2_densities

data1_Bim_c1c2_densities.pdf: \
        $(PLOTSDIR)/bid_bim_nbd_release/plot_k1_k2_dists.py \
        $(PLOTSDIR)/bid_bim_nbd_release/process_3conf_mcmc.py \
        $(NBD_MCMC_DIR)/Bim_c1_density_mx.txt \
        $(NBD_MCMC_DIR)/Bim_c2_density_mx.txt
	python -m tbidbaxlipo.plots.bid_bim_nbd_release.plot_k1_k2_dists plot c1c2 \
		$(NBD_MCMC_DIR)/Bim_c1_density_mx.txt \
		$(NBD_MCMC_DIR)/Bim_c2_density_mx.txt \
		data1_Bim_c1c2_densities


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

# --- NBD endpoints, dataset 2, normalized to F/F0 ---
data2_nbd_endpt_norm.pdf: \
        $(PLOTSDIR)/nbd_bax_analysis.py \
        $(PLOTSDIR)/bid_bim_fret_nbd_release/plots.py \
        $(NBD_FRET_DATA)
	python -m tbidbaxlipo.plots.bid_bim_fret_nbd_release.plots nbd_endpoint_norm $(basename $@)

# --- NBD endpoints, dataset 2, no normalization (raw values) ---
data2_nbd_endpt_no_norm.pdf: \
        $(PLOTSDIR)/nbd_bax_analysis.py \
        $(PLOTSDIR)/bid_bim_fret_nbd_release/plots.py \
        $(NBD_FRET_DATA)
	python -m tbidbaxlipo.plots.bid_bim_fret_nbd_release.plots nbd_endpoint_no_norm $(basename $@)

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

data2_fret_endpt.pdf: \
        $(PLOTSDIR)/nbd_bax_analysis.py \
        $(PLOTSDIR)/bid_bim_fret_nbd_release/plots.py \
        $(PLOTSDIR)/bid_bim_fret_nbd_release/preprocess_data.py \
        $(NBD_FRET_DATA)
	python -m tbidbaxlipo.plots.bid_bim_fret_nbd_release.plots fret_endpoint $(basename $@)

data2_example_curves.pdf: \
        $(PLOTSDIR)/bid_bim_fret_nbd_release/plot_example_fret_curves.py \
        $(NBD_FRET_DATA)
	python -m tbidbaxlipo.plots.bid_bim_fret_nbd_release.plot_example_fret_curves

data2_example_nbd_fret_fits.pdf: \
        $(PLOTSDIR)/bid_bim_fret_nbd_release/plot_example_fret_fits.py \
        $(NBD_FRET_DATA)
	python -m tbidbaxlipo.plots.bid_bim_fret_nbd_release.plot_example_fret_fits

data2_delta_fret.pdf: \
        $(PLOTSDIR)/nbd_bax_analysis.py \
        $(PLOTSDIR)/bid_bim_fret_nbd_release/plot_delta_fret.py \
        $(PLOTSDIR)/bid_bim_fret_nbd_release/preprocess_data.py \
        $(NBD_FRET_DATA)
	python -m tbidbaxlipo.plots.bid_bim_fret_nbd_release.plot_delta_fret $(basename $@)


# ---------------------------------------------------------------

nbd_labeling_ratios_table: \
        $(PLOTSDIR)/bid_bim_nbd_release/nbd_labeling_ratios_table.py \
        $(NBD_DATA)
	python -m tbidbaxlipo.plots.bid_bim_nbd_release.nbd_labeling_ratios_table \
                      nbd_labeling_ratios_table.csv

141203_example_nbd_fret_fits.pdf: \
    $(PLOTSDIR)/x141203_Bax_Bax_FRET/plot_example_fret_fits.py
	python -m tbidbaxlipo.plots.x141203_Bax_Bax_FRET.plot_example_fret_fits


# ---------------------------------------------------------------

# --- Normalized release endpoint data for release, dataset 3 ---
# Not done, no WT control
#data3_release_endpt_norm.pdf: \
#        $(PLOTSDIR)/nbd_bax_analysis.py \
#        $(PLOTSDIR)/bax_bax_fret_nbd_release/plots.py \
#        $(BAX_FRET_DATA)
#	python -m tbidbaxlipo.plots.bax_bax_fret_nbd_release.plots release_endpoint_norm $(basename $@)

# --- Raw release endpoint data for release, dataset 3 ---
data3_release_endpt_no_norm.pdf: \
        $(PLOTSDIR)/nbd_bax_analysis.py \
        $(PLOTSDIR)/bax_bax_fret_nbd_release/plots.py \
        $(BAX_FRET_DATA)
	python -m tbidbaxlipo.plots.bax_bax_fret_nbd_release.plots release_endpoint_no_norm $(basename $@)

# --- NBD endpoints, dataset 3, normalized to F/F0 ---
data3_nbd_endpt_norm.pdf: \
        $(PLOTSDIR)/nbd_bax_analysis.py \
        $(PLOTSDIR)/bax_bax_fret_nbd_release/plots.py \
        $(BAX_FRET_DATA)
	python -m tbidbaxlipo.plots.bax_bax_fret_nbd_release.plots nbd_endpoint_norm $(basename $@)

# --- NBD endpoints, dataset 3, no normalization (raw values) ---
# Not done, because data is already normalized in the spreadsheet
#data3_nbd_endpt_no_norm.pdf: \
#        $(PLOTSDIR)/nbd_bax_analysis.py \
#        $(PLOTSDIR)/bax_bax_fret_nbd_release/plots.py \
#        $(BAX_FRET_DATA)
#	python -m tbidbaxlipo.plots.bax_bax_fret_nbd_release.plots nbd_endpoint_no_norm $(basename $@)

data3_fret_endpt.pdf: \
        $(PLOTSDIR)/nbd_bax_analysis.py \
        $(PLOTSDIR)/bax_bax_fret_nbd_release/plots.py \
        $(PLOTSDIR)/bax_bax_fret_nbd_release/preprocess_data.py \
        $(BAX_FRET_DATA)
	python -m tbidbaxlipo.plots.bax_bax_fret_nbd_release.plots fret_endpoint $(basename $@)

data3_example_curves.pdf: \
        $(PLOTSDIR)/bax_bax_fret_nbd_release/plot_example_fret_curves.py \
        $(BAX_FRET_DATA)
	python -m tbidbaxlipo.plots.bax_bax_fret_nbd_release.plot_example_fret_curves

# FRET/NBD/release, initial rates
# Building this target creates four files:
#     data3_init_nbd_bar.pdf,
#     data3_init_nbd_bar.pdf,
#     data3_init_scatter_norm_Bid.pdf,
#     data3_init_scatter_norm_Bim.pdf
#data3_init_nbd_bar.pdf: \
#        $(PLOTSDIR)/nbd_bax_analysis.py \
#        $(PLOTSDIR)/bax_bax_fret_nbd_release/plots.py \
#        $(BAX_FRET_DATA)
#	python -m tbidbaxlipo.plots.bax_bax_fret_nbd_release.plots initial_rate_samples data3


#data3_example_nbd_fret_fits.pdf: \
#        $(PLOTSDIR)/bax_bax_fret_nbd_release/plot_example_fret_fits.py \
#        $(BAX_FRET_DATA)
#	python -m tbidbaxlipo.plots.bax_bax_fret_nbd_release.plot_example_fret_fits
#
#data3_delta_fret.pdf: \
#        $(PLOTSDIR)/nbd_bax_analysis.py \
#        $(PLOTSDIR)/bax_bax_fret_nbd_release/plot_delta_fret.py \
#        $(PLOTSDIR)/bax_bax_fret_nbd_release/preprocess_data.py \
#        $(BAX_FRET_DATA)
#	python -m tbidbaxlipo.plots.bax_bax_fret_nbd_release.plot_delta_fret $(basename $@)


