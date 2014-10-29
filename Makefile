FIGDIR := results/figures
DATADIR := tbidbaxlipo/data
CODEDIR := tbidbaxlipo

#figures: ~/matplotlib/.matplotlibrc

$(FIGDIR)/fig_141016_1.pdf: \
          $(DATADIR)/141016_Bax_depletion_preincubation.txt \
          $(DATADIR)/141016_Bax_depletion_timecourse.txt \
          $(DATADIR)/141016_Bax_depletion_triton.txt \
          $(DATADIR)/141016_Bax_depletion_added_ANTS_EF.txt \
          $(CODEDIR)/util/plate_assay.py \
          $(CODEDIR)/util/__init__.py \
          $(CODEDIR)/plots/layout_141016.py
	python $(CODEDIR)/plots/layout_141016.py
	mv fig_141016_1.pdf $(FIGDIR)

$(MCMCDIR)/141016_2conf.pck:
# the code
# the data
# the models
