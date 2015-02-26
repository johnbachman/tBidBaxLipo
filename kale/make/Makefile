CODEDIR := ../../tbidbaxlipo

all: nbd_release fret_nbd_release

clean:
	rm -f ./*.pdf ./*.png

deploy:
	cp ./*.pdf ~/Dropbox/Bax_insertion_paper/figure_panels_JAB

.PHONY: clean deploy nbd_release fret_nbd_release

nbd_release: nbd_release_endpts
nbd_release_endpts: data1_release_endpt.pdf data1_nbd_endpt.pdf

fret_nbd_release: fret_nbd_release_endpts
fret_nbd_release_endpts: data2_release_endpt.pdf data2_nbd_endpt.pdf

# --- Raw data for NBD/release, dataset 1 ---
# data1_release_endpt.pdf, data1_nbd_endpt.pdf
data1_release_endpt.pdf:
	python -m tbidbaxlipo.plots.bid_bim_nbd_release.plots release_endpoint $(basename $@)

data1_nbd_endpt.pdf:
	python -m tbidbaxlipo.plots.bid_bim_nbd_release.plots nbd_endpoint $(basename $@)

# --- Raw data for NBD/release, dataset 2 ---
# data2_release_endpt.pdf, data2_nbd_endpt.pdf
data2_release_endpt.pdf:
	python -m tbidbaxlipo.plots.bid_bim_fret_nbd_release.plots release_endpoint $(basename $@)

data2_nbd_endpt.pdf:
	python -m tbidbaxlipo.plots.bid_bim_fret_nbd_release.plots nbd_endpoint $(basename $@)


