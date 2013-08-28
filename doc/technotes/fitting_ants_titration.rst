Fitting dye release kinetics by MCMC
====================================

Local Fitting
-------------

1. After the dye release experiment, save the files as `.csv`. If the kinetics
   are by plate (typical for slower kinetics, like Bax assays), then it is
   necessary to sort the spreadsheet by the `Well` column, then save as `.csv`.
   Save the `.csv` as a Windows `.csv`. Do this for both the timecourse and the
   Triton Excel files. If the data is kinetics by well, then the timecourses
   don't need to be sorted first.

2. Make the data a reusable, importable, and plottable resource by creating a
   layout file. The layout file can live in :py:mod:`tbidbaxlipo.plots`. It
   should have a ``plot()`` function to run all plots, and it should create
   top-level package variables creating the data in various forms (normalized,
   averaged, reset to zero, etc.) For fitting by MCMC, create a top-level
   variable for a pandas dataframe (with the same format as the one created by
   the function :py:func:`tbidbaxlipo.util.plate_assay.to_dataframe`).

3. Tweak the main function of :py:mod:`tbidbaxlipo.mcmc.pore_mcmc` to use the
   desired dataset.  Set the step size, hessian scale or other parameters as
   desired.  The `get_likelihood_function` may need tweaking, or it may not.

4. Run `pore_mcmc`, creating the .mcmc file. For example::

    mkdir pore_fits
    cd pore_fits
    python -m tbidbaxlipo.mcmc.pore_mcmc random_seed=0 model=bax_schwarz \
        cpt_type=one_cpt nsteps=10000

5. Modify :py:mod:`tbidbaxlipo.run_report` to contain the list of reports you
   wish to run on the MCMC chain. Run with::

    python -m tbidbaxlipo.run_report *.mcmc

   Note that this will dump out `.html`, `.png`, and other files, so this
   should be run in the directory that will contain the report.

6. View the report by opening the generated `index.html` file.

Fitting on Orchestra
--------------------

Perform steps 1-3 as above.

4. Tweak the script :py:mod:`tbidbaxlipo.mcmc.pore_mcmc_jobs` to use the
   appropriate number of steps, number of chains per model, and desired models
   and compartmentalization types. Then go to a directory on the server where the `.mcmc`
   files are to be created and run::

    cd /hms/scratch1/jab69
    mkdir some_pore_fits
    cd some_pore_fits
    python -m tbidbaxlipo.mcmc.pore_mcmc_jobs

5. When the jobs are done running, modify :py:mod:`tbidbaxlipo.run_report` to contain
   the list of reports you wish to run on the chains. Switch to a web accessible
   directory and run (note that since it may take a long time to run all reports,
   it's best to submit this as a job)::

    cd /www/sorger.med.harvard.edu/docroot/data/bachman
    mkdir some_pore_fits
    cd some_pore_fits
    bsub -q short -W 12:00 python -m tbidaxlipo.run_report /hms/scratch1/jab69/some_pore_fits/*.mcmc

6. View the report by opening a browser and navigating to the report, e.g., at
   ``http://sorger.med.harvard.edu/data/bachman/some_pore_fits/index.html``.

