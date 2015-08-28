#nosetests tbidbaxlipo.plots.nbd_analysis tbidbaxlipo.plots.grid_analysis
#nosetests tbidbaxlipo.tests.test_one_cpt_models
#nosetests tbidbaxlipo.tests.test_multi_cpt_models
#nosetests --with-doctest tbidbaxlipo.util.plate_assay
export BASEDIR=`pwd`
cd $BASEDIR/kale/make
make
cd $BASEDIR/kale/doc
make html
cd $BASEDIR
nosetests tbidbaxlipo/plots/x140318_Bax_liposome_titration
nosetests tbidbaxlipo/plots/x140320_NBD_Bax_BimBH3_unlab_Bax_titration
nosetests tbidbaxlipo/models
make figures
