from tbidbaxlipo.util.report import Report
from nbd_analysis import plot_raw, plot_normalized, plot_fit

rep = Report()
plot_raw(report=rep)
rep.writeReport('nbd_raw')

rep = Report()
plot_normalized(report=rep)
rep.writeReport('nbd_normalized_by_fit')

rep = Report()
plot_fit(report=rep)
rep.writeReport('nbd_double_exp')

