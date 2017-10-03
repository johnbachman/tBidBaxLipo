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

* Building "stoplight" plot (do for both `fret_mcmc` and `fret_mcmc_norm`):
  * `python process_3conf_mcmc_submit.py fret_mcmc/*`. Run from
    `bax_bax_fret_nbd_release` directory to generate density files after fitting.
  * (From `fret_mcmc` directory): `python ../plot_k1_k2_dists.py assemble ...`
    * `... k1 *.k1_hist`
    * `... k2 *.k2_hist`
    * `... c1 *.c1_scaling`
    * `... c2 *.c2_scaling`
    * `... f1 *.fret1_scaling`
    * `... f2 *.fret2_scaling`
  * (From `fret_mcmc` directory): `python ../plot_k1_k2_dists.py plot ...`
    * `... k1k2 fret_Bid_k1_density_mx.txt fret_Bid_k2_density_mx.txt pt_data3_k1k2`
    * `... c1c2 fret_Bid_c1_density_mx.txt fret_Bid_c2_density_mx.txt pt_data3_c1c2`
    * `... f1f2 fret_Bid_f1_density_mx.txt fret_Bid_f2_density_mx.txt pt_data3_f1f2`

```
python ../plot_k1_k2_dists.py assemble k1 *.k1_hist
python ../plot_k1_k2_dists.py assemble k2 *.k2_hist
python ../plot_k1_k2_dists.py assemble c1 *.c1_scaling
python ../plot_k1_k2_dists.py assemble c2 *.c2_scaling
python ../plot_k1_k2_dists.py assemble f1 *.fret1_scaling
python ../plot_k1_k2_dists.py assemble f2 *.fret2_scaling

# For fret_mcmc:
python ../plot_k1_k2_dists.py plot k1k2 fret_Bid_k1_density_mx.txt fret_Bid_k2_density_mx.txt pt_data3_k1k2
python ../plot_k1_k2_dists.py plot c1c2 fret_Bid_c1_density_mx.txt fret_Bid_c2_density_mx.txt pt_data3_c1c2
python ../plot_k1_k2_dists.py plot f1f2 fret_Bid_f1_density_mx.txt fret_Bid_f2_density_mx.txt pt_data3_f1f2

# For fret_mcmc_norm:
python ../plot_k1_k2_dists.py plot k1k2 fret_norm_Bid_k1_density_mx.txt fret_norm_Bid_k2_density_mx.txt pt_data3_norm_k1k2
python ../plot_k1_k2_dists.py plot c1c2 fret_norm_Bid_c1_density_mx.txt fret_norm_Bid_c2_density_mx.txt pt_data3_norm_c1c2
python ../plot_k1_k2_dists.py plot f1f2 fret_norm_Bid_f1_density_mx.txt fret_norm_Bid_f2_density_mx.txt pt_data3_norm_f1f2
```

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

