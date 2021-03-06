SHELL=/bin/bash -O extglob -c

FIGDIR := results/figures/panels
DATADIR := tbidbaxlipo/data
CODEDIR := tbidbaxlipo
MCMCDIR := results/mcmc
PQUEUE := sorger_par_unlimited

x140320 := $(CODEDIR)/plots/x140320_NBD_Bax_BimBH3_unlab_Bax_titration/mcmc
x140318 := $(CODEDIR)/plots/x140318_Bax_liposome_titration/mcmc
x141119 := $(CODEDIR)/plots/x141119_Bax_Bid_saturation/mcmc
x141203 := $(CODEDIR)/plots/x141203_Bax_Bax_FRET/mcmc

figures: \
		$(FIGDIR)/141119_Bid_20nm_timecourses.pdf \
		$(FIGDIR)/fig_141016_1.pdf \
		$(FIGDIR)/fig_fmax_fit_comparison.pdf \
		$(FIGDIR)/poisson_bax_fmax.pdf \
		$(FIGDIR)/slice_bax_fixed.pdf \
		$(FIGDIR)/140320_exp_fits.pdf \
		$(FIGDIR)/140318_exp_fits_lstsq_fmax_var.pdf \
		$(FIGDIR)/140429_exact_comp_bind_fit.pdf \
		$(FIGDIR)/140429_gouy_chap_fit.pdf \
		$(FIGDIR)/model_predictions_bax_titration.pdf \
		$(FIGDIR)/model_predictions_lipo_titration.pdf \
		$(FIGDIR)/140724_requench_bid_2.pdf \
		$(FIGDIR)/140710_requench_bax.pdf \
		$(FIGDIR)/requenching_examples.pdf \
		$(FIGDIR)/140318_evidence_barplot1.pdf \
		$(FIGDIR)/140320_evidence_barplot.pdf
#		$(FIGDIR)/pt_140318_1c_Baxtr1Activ1Rever1Nbd1.mcmc.ibax_reverse.pdf

deploy:
	rsync -av results/figures/ ~/Dropbox/Bachman-Sorger\ Talks\ and\ Papers/Bachman-Kale\ Bax\ kinetics/figures

mcmc_figures: $(FIGDIR)/pt_140318_nbd_2_conf_fits.pdf

pt_140318_figures:
	for f in $(x140318)/*.mcmc; \
	do \
		OUTPUT=$(FIGDIR)/pt_140318/$$(basename -s .mcmc $$f) ;\
		mkdir -p $$OUTPUT ;\
		qsub -b y -cwd -V python $(CODEDIR)/pt/show_chain.py $$f $$OUTPUT ;\
	done

pt_140320_figures:
	for f in $(x140320)/*.mcmc; \
	do \
		OUTPUT=$(FIGDIR)/pt_140320/$$(basename -s .mcmc $$f) ;\
		mkdir -p $$OUTPUT ;\
		qsub -b y -cwd -V python $(CODEDIR)/pt/show_chain.py $$f $$OUTPUT ;\
	done

#python $(CODEDIR)/pt/show_chain.py $$f $$OUTPUT ;\


clean:
	cd $(FIGDIR); rm -f *.pdf

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

.PHONY: clean pt_140320 pt_140318 pt_141119 pt_141203_54C pt_141203_126C

# Keying the MCMC execution on the dependency file ensures that a re-run is
# performed if the model ensemble specification (or fit parameters) has
# changed. Also ensures that the dependency file exists before we try to
# include it. As a convention, the target used here should correspond to the
# basename of the .yaml file specifying the fit for the model ensemble
# (e.g., pt_140320.yaml --> pt_140320)
pt_140320: $(x140320)/pt_140320.deps.txt
# This file specifies the list of *.mcmc files that pt_140320 depends on
-include $(x140320)/pt_140320.deps.txt

# See comments for pt_140320, above
pt_140318: $(x140318)/pt_140318.deps.txt
-include $(x140318)/pt_140318.deps.txt

# See comments for pt_140320, above
pt_141119: $(x141119)/pt_141119.deps.txt
-include $(x141119)/pt_141119.deps.txt

# See comments for pt_140320, above
pt_141203_54C: $(x141203)/pt_141203_54C.deps.txt
-include $(x141203)/pt_141203_54C.deps.txt

# See comments for pt_140320, above
pt_141203_126C: $(x141203)/pt_141203_126C.deps.txt
-include $(x141203)/pt_141203_126C.deps.txt

# Running this script generates both the dependency list and the .yaml files
# specifying the fit parameters for each individual model
%.deps.txt: %.fit.ensemble
	python -m tbidbaxlipo.pt.generate_model_ensemble_fit_files $<

# In this case, we know that the .fit file exists, since it has been
# regenerated in the step that regenerated the dependency file. What
# we care about is whether the hash is any different.
%.mcmc: $(call to-md5, %.fit)
	python -m tbidbaxlipo.pt.compile_models $(call from-md5, $<)
	#bsub -q $(PQUEUE) -n 50 -W 12:00 -a openmpi mpirun.lsf python -m tbidbaxlipo.pt.run_pt $(call from-md5, $<) 1
	qsub -b y -cwd -V -o $(call from-md5, $<).out -e $(call from-md5, $<).err -pe orte 128 mpirun python $(CODEDIR)/pt/run_pt.py $(call from-md5, $<) 1 $(call from-md5, $<).pos

# Don't delete the .md5 files! Without this rule they are treated as
# intermediates and deleted.
.PRECIOUS: %.md5

# ==== FIGURES ==============================================

$(FIGDIR)/141119_Bid_20nm_timecourses.pdf: \
        $(CODEDIR)/plots/x141119_Bax_Bid_saturation/plot_exp_fits.py
	python $(CODEDIR)/plots/x141119_Bax_Bid_saturation/plot_exp_fits.py
	mv *.pdf $(FIGDIR)
	mv *.png $(FIGDIR)


# --- Bax depletion figures ----
# fig_141016_1.pdf, fig_141016_2.pdf
$(FIGDIR)/fig_141016_1.pdf: \
		$(DATADIR)/141016_Bax_depletion_preincubation.txt \
		$(DATADIR)/141016_Bax_depletion_timecourse.txt \
		$(DATADIR)/141016_Bax_depletion_triton.txt \
		$(DATADIR)/141016_Bax_depletion_added_ANTS_EF.txt \
		$(CODEDIR)/util/plate_assay.py \
		$(CODEDIR)/util/__init__.py \
		$(CODEDIR)/plots/layout_141016.py
	python $(CODEDIR)/plots/layout_141016.py
	mv *.pdf $(FIGDIR)

# --- Figs generated by hazard_rate.py ----
# fig_fmax_fit_comparison.pdf, fig_hazard_rate.pdf
$(FIGDIR)/fig_fmax_fit_comparison.pdf: \
		$(CODEDIR)/plots/layout_140311.py \
		$(CODEDIR)/plots/hazard_rate.py \
		$(CODEDIR)/util/__init__.py \
		$(CODEDIR)/util/fitting.py
	python $(CODEDIR)/plots/hazard_rate.py
	mv *.pdf $(FIGDIR)

# --- Figures showing the relationship between Bax distribution and fraction
#     permeabilized ----
# poisson_bax_fmax.pdf, poisson_bax_fmax_fit.pdf
$(FIGDIR)/poisson_bax_fmax.pdf: \
		$(CODEDIR)/plots/poisson_bax_fmax.py
	python $(CODEDIR)/plots/poisson_bax_fmax.py
	mv *.pdf $(FIGDIR)

# --- Fits of models to 140318 Bax-NBD/liposome titration data by MCMC ----
# Makes plots pt..._fits.pdf, pt..._tri.pdf
# (for 2_conf, 2_conf_rev, and 3_conf models)
$(FIGDIR)/pt_140318_nbd_2_conf_fits.pdf: \
		$(CODEDIR)/plots/x140318_Bax_liposome_titration/pt_140318_nbd_2_conf_5.mcmc \
		$(CODEDIR)/plots/x140318_Bax_liposome_titration/pt_140318_nbd_2_conf_rev_4.mcmc \
		$(CODEDIR)/plots/x140318_Bax_liposome_titration/pt_140318_nbd_3_conf_3.mcmc \
		$(CODEDIR)/plots/x140318_Bax_liposome_titration/plot_mcmc_model_fits.py
	python $(CODEDIR)/plots/x140318_Bax_liposome_titration/plot_mcmc_model_fits.py \
			$(CODEDIR)/plots/x140318_Bax_liposome_titration/pt_140318_nbd_2_conf_5.mcmc pt_140318_nbd_2_conf
	python $(CODEDIR)/plots/x140318_Bax_liposome_titration/plot_mcmc_model_fits.py \
			$(CODEDIR)/plots/x140318_Bax_liposome_titration/pt_140318_nbd_2_conf_rev_4.mcmc pt_140318_nbd_2_conf_rev
	python $(CODEDIR)/plots/x140318_Bax_liposome_titration/plot_mcmc_model_fits.py \
			$(CODEDIR)/plots/x140318_Bax_liposome_titration/pt_140318_nbd_3_conf_3.mcmc pt_140318_nbd_3_conf
	mv *.pdf $(FIGDIR)

# --- Slice diagrams ---
$(FIGDIR)/slice_bax_fixed.pdf: \
		$(CODEDIR)/plots/slice_diagrams.py \
		$(CODEDIR)/util/__init__.py
	python $(CODEDIR)/plots/slice_diagrams.py
	mv *.pdf $(FIGDIR)

# --- Exponential fits to Bax titration, 140320
#     140320_exp_fit_curves.pdf, 140320_exp_fits.pdf ---
$(FIGDIR)/140320_exp_fits.pdf: \
		$(CODEDIR)/plots/x140320_NBD_Bax_BimBH3_unlab_Bax_titration/preprocess_data.py \
		$(CODEDIR)/plots/x140320_NBD_Bax_BimBH3_unlab_Bax_titration/exp_fits_lstsq.py \
		$(CODEDIR)/data/140320_NBD_Bax_BimBH3_unlab_Bax_titration.txt \
		$(CODEDIR)/util/fitting.py \
		$(CODEDIR)/util/plate_assay.py \
		$(CODEDIR)/util/__init__.py
	python $(CODEDIR)/plots/x140320_NBD_Bax_BimBH3_unlab_Bax_titration/exp_fits_lstsq.py
	mv *.pdf $(FIGDIR)

# --- 140318_exp_fits_lstsq_fmax_var.pdf,
#     140318_exp_fits_lstsq_curves_fmax_var.pdf,
#     140318_exp_fits_lstsq_fmax_fixed.pdf,
#     140318_exp_fits_lstsq_curves_fmax_fixed.pdf ---
$(FIGDIR)/140318_exp_fits_lstsq_fmax_var.pdf: \
		$(CODEDIR)/plots/x140318_Bax_liposome_titration/exp_fits_lstsq.py \
		$(CODEDIR)/plots/x140318_Bax_liposome_titration/preprocess_data.py \
		$(CODEDIR)/util/plate_assay.py \
		$(CODEDIR)/data/140318_NBD_Bax_BimBH3_lipo_titration.txt \
		$(CODEDIR)/util/__init__.py
	python $(CODEDIR)/plots/x140318_Bax_liposome_titration/exp_fits_lstsq.py
	mv *.pdf $(FIGDIR)


# --- Fits of 140429 Bid FRET competition experiment by MCMC to the exact
#     competition binding model ----
# Files: 140429_exact_comp_bind_fit.pdf, 140429_exact_comp_bind_marginals.pdf
$(FIGDIR)/140429_exact_comp_bind_fit.pdf: \
		$(CODEDIR)/plots/x140429_Bid_membrane_FRET/exact_comp_bind_mcmc.py \
		$(CODEDIR)/plots/x140429_Bid_membrane_FRET/calculate_fret.py \
		$(CODEDIR)/data/parse_140429_Bid_membrane_FRET.py \
		$(CODEDIR)/data/140429_Bid_membrane_FRET.xlsx \
		$(CODEDIR)/util/__init__.py \
		$(CODEDIR)/plots/x140429_Bid_membrane_FRET/140429_exact_comp_bind.mcmc
	python $(CODEDIR)/plots/x140429_Bid_membrane_FRET/exact_comp_bind_mcmc.py plot \
		   $(CODEDIR)/plots/x140429_Bid_membrane_FRET/140429_exact_comp_bind.mcmc
	mv *.pdf $(FIGDIR)

$(CODEDIR)/plots/x140429_Bid_membrane_FRET/140429_exact_comp_bind.mcmc: \
		$(CODEDIR)/plots/x140429_Bid_membrane_FRET/exact_comp_bind_mcmc.py \
		$(CODEDIR)/plots/x140429_Bid_membrane_FRET/calculate_fret.py \
		$(CODEDIR)/data/parse_140429_Bid_membrane_FRET.py \
		$(CODEDIR)/data/140429_Bid_membrane_FRET.xlsx
	python $(CODEDIR)/plots/x140429_Bid_membrane_FRET/exact_comp_bind_mcmc.py sample $(CODEDIR)/plots/x140429_Bid_membrane_FRET/140429_exact_comp_bind.mcmc

# Files: 140429_gouy_chap_fit.pdf, 140429_gouy_chap_marginals.pdf
$(FIGDIR)/140429_gouy_chap_fit.pdf: \
		$(CODEDIR)/plots/x140429_Bid_membrane_FRET/gouy_chap_mcmc.py \
		$(CODEDIR)/plots/x140429_Bid_membrane_FRET/calculate_fret.py \
		$(CODEDIR)/data/parse_140429_Bid_membrane_FRET.py \
		$(CODEDIR)/data/140429_Bid_membrane_FRET.xlsx \
		$(CODEDIR)/util/__init__.py \
		$(CODEDIR)/plots/x140429_Bid_membrane_FRET/140429_gouy_chap.mcmc
	python $(CODEDIR)/plots/x140429_Bid_membrane_FRET/gouy_chap_mcmc.py plot \
		   $(CODEDIR)/plots/x140429_Bid_membrane_FRET/140429_gouy_chap.mcmc
	mv *.pdf $(FIGDIR)

$(CODEDIR)/plots/x140429_Bid_membrane_FRET/140429_gouy_chap.mcmc: \
		$(CODEDIR)/plots/x140429_Bid_membrane_FRET/gouy_chap_mcmc.py \
		$(CODEDIR)/plots/x140429_Bid_membrane_FRET/calculate_fret.py \
		$(CODEDIR)/data/parse_140429_Bid_membrane_FRET.py \
		$(CODEDIR)/data/140429_Bid_membrane_FRET.xlsx
	python $(CODEDIR)/plots/x140429_Bid_membrane_FRET/gouy_chap_mcmc.py sample $(CODEDIR)/plots/x140429_Bid_membrane_FRET/140429_gouy_chap.mcmc

$(FIGDIR)/model_predictions_bax_titration.pdf: \
        $(CODEDIR)/plots/model_predictions_bax_titration.py
	python $(CODEDIR)/plots/model_predictions_bax_titration.py
	mv *.pdf $(FIGDIR)
	mv *.png $(FIGDIR)

$(FIGDIR)/model_predictions_lipo_titration.pdf: \
        $(CODEDIR)/plots/model_predictions_lipo_titration.py
	python $(CODEDIR)/plots/model_predictions_lipo_titration.py
	mv *.pdf $(FIGDIR)
	mv *.png $(FIGDIR)

$(FIGDIR)/140724_requench_bid_2.pdf: \
        $(CODEDIR)/plots/layout_140724.py \
        $(CODEDIR)/util/dpx_assay.py
	python $(CODEDIR)/plots/layout_140724.py
	mv *.pdf $(FIGDIR)
	mv *.png $(FIGDIR)

$(FIGDIR)/140710_requench_bax.pdf: \
        $(CODEDIR)/plots/layout_140710.py \
        $(CODEDIR)/util/dpx_assay.py
	python $(CODEDIR)/plots/layout_140710.py
	mv *.pdf $(FIGDIR)
	mv *.png $(FIGDIR)

$(FIGDIR)/requenching_examples.pdf: $(CODEDIR)/util/dpx_assay.py
	python $(CODEDIR)/util/dpx_assay.py
	mv *.pdf $(FIGDIR)
	mv *.png $(FIGDIR)

$(FIGDIR)/140318_evidence_barplot1.pdf: \
        $(CODEDIR)/plots/x140318_Bax_liposome_titration/evidence_list.pck \
        $(CODEDIR)/plots/x140318_Bax_liposome_titration/plot_evidence_barplot.py
	python $(CODEDIR)/plots/x140318_Bax_liposome_titration/plot_evidence_barplot.py $(CODEDIR)/plots/x140318_Bax_liposome_titration/evidence_list.pck
	mv *.pdf $(FIGDIR)
	mv *.png $(FIGDIR)

$(FIGDIR)/140320_evidence_barplot.pdf: \
        $(CODEDIR)/plots/x140320_NBD_Bax_BimBH3_unlab_Bax_titration/evidence_list.pck \
        $(CODEDIR)/plots/x140320_NBD_Bax_BimBH3_unlab_Bax_titration/plot_evidence_barplot.py
	python $(CODEDIR)/plots/x140320_NBD_Bax_BimBH3_unlab_Bax_titration/plot_evidence_barplot.py $(CODEDIR)/plots/x140320_NBD_Bax_BimBH3_unlab_Bax_titration/evidence_list.pck
	mv *.pdf $(FIGDIR)
	mv *.png $(FIGDIR)

$(FIGDIR)/pt_140318_1c_Baxtr1Activ1Rever1Nbd1.mcmc.ibax_reverse.pdf: \
        $(CODEDIR)/plots/x140318_Bax_liposome_titration/plot_m2_m3_posteriors.pdf
	python $(CODEDIR)/plots/x140318_Bax_liposome_titration/plot_m2_m3_posteriors.pdf
