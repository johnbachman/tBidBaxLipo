CODEDIR := ../../tbidbaxlipo

all: nbd_release_endpts

clean:
	rm -f ./*.pdf ./*.png

.PHONY: clean

nbd_release: nbd_release_endpts

nbd_release_endpts: data1_release_endpt.pdf data1_nbd_endpt.pdf

# --- Raw data for NBD/release, dataset 1 ---
# data1_release_endpt.pdf, data1_nbd_endpt.pdf
data1_release_endpt.pdf:
	python -m tbidbaxlipo.plots.bid_bim_nbd_release.plots release_endpoint data1_release_endpt

data1_nbd_endpt.pdf:
	python -m tbidbaxlipo.plots.bid_bim_nbd_release.plots nbd_endpoint data1_nbd_endpt


