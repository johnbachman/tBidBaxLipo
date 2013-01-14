from pysb import *
import nbd_model_shared

Model()

nbd_model_shared.declare_shared_components()

Parameter('c3_insertion_rate', 0.004)
Parameter('c62_insertion_rate', 0.001)
Parameter('c120_insertion_rate', 0.0017)
Parameter('c122_insertion_rate', 0.0012)
Parameter('c126_insertion_rate', 0.002)

Rule('c3_insertion', Bax(c3='s') >> Bax(c3='m'), c3_insertion_rate)
Rule('c62_insertion', Bax(c62='s') >> Bax(c62='m'), c62_insertion_rate)
Rule('c120_insertion', Bax(c120='s') >> Bax(c120='m'), c120_insertion_rate)
Rule('c122_insertion', Bax(c122='s') >> Bax(c122='m'), c122_insertion_rate)
Rule('c126_insertion', Bax(c126='s') >> Bax(c126='m'), c126_insertion_rate)

