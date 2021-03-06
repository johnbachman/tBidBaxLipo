#nosetests tbidbaxlipo.plots.nbd_analysis tbidbaxlipo.plots.grid_analysis
#nosetests tbidbaxlipo.tests.test_one_cpt_models
#nosetests tbidbaxlipo.tests.test_multi_cpt_models
#nosetests --with-doctest tbidbaxlipo.util.plate_assay

# Test the code
nosetests tbidbaxlipo/plots/x140318_Bax_liposome_titration
nosetests tbidbaxlipo/plots/x140320_NBD_Bax_BimBH3_unlab_Bax_titration
nosetests tbidbaxlipo/models
# Test main figures
make figures
export BASEDIR=`pwd`
# Test kale paper main text figures
cd $BASEDIR/kale/make
make
# Test kale paper supplemental figures
cd $BASEDIR/kale/doc
make html
cd $BASEDIR
