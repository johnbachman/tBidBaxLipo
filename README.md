# tBidBaxLipo: Analysis of Bax Interactions from in vitro Measurements

## Make instructions for Kale NBD/FRET datasets

Run the following make commands from `kale/make`.

### Dataset 1 (NBD-Bax)

* `make pt_data1`: Run MCMC fitting (requires Starcluster setup). Files are
  stored in `tbidbaxlipo/plots/bid_bim_nbd_release/mcmc`

* `make pt_data1_norm`: Run MCMC fitting with a Gaussian prior distribution for
  the NBD fluorescence scaling parameters (requires Starcluster setup).
  Files are stored in `tbidbaxlipo/plots/bid_bim_nbd_release/mcmc_norm_prior`

* `make nbd_labeling_ratios_table`

### Dataset 2 (NBD-Bax + Bid/Bax FRET)

* `make pt_data2_fret`: Run MCMC fitting with coupled NBD and FRET dynamics
  (requires Starcluster setup). Files are stored in
  `tbidbaxlipo/plots/bid_bim_fret_nbd_release/fret_mcmc`. (54 fits)

* `make pt_data2_fret_norm`: Run MCMC fitting with coupled NBD and FRET
  dynamics with a Gaussian prior distribution for the NBD fluorescence
  scaling parameters (requires Starcluster setup). Files are stored in
  `tbidbaxlipo/plots/bid_bim_fret_nbd_release/fret_mcmc_norm`. (54 fits)

### Dataset 3 (NBD-Bax + Bax/Bax FRET)

* `make pt_data3_fret`: Run MCMC fitting with coupled NBD and FRET dynamics
  (requires Starcluster setup). Files are stored in
  `tbidbaxlipo/plots/bax_bax_fret_nbd_release/fret_mcmc`. (27 fits)

* `make pt_data3_fret_norm`: Run MCMC fitting with coupled NBD and FRET
  dynamics with a Gaussian prior distribution for the NBD fluorescence
  scaling parameters (requires Starcluster setup). Files are stored in
  `tbidbaxlipo/plots/bax_bax_fret_nbd_release/fret_mcmc_norm`. (27 fits)

## All others

* `make all`: Equivalent to building all of the following targets:

  * `make bax_bax_fret`: Bax-Bax FRET Pilot Dataset (141203)

  * `make nbd_release`: Generates the following plots:
    * Raw and normalized NBD fluorescence endpoint bar plots
    * Raw and WT-normalized permeabilization endpoint bar plots
    * Others...

  * `make fret_nbd_release`: Generates the following plots:
    * Raw and normalized NBD fluorescence endpoint bar plots
    * Raw and WT-normalized permeabilization endpoint bar plots
    * Others...

