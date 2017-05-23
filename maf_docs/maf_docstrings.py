from __future__ import print_function
import lsst.sims.maf.metrics as metrics
from inspect import cleandoc
import pandas as pd

metricNames = [metric for metric in dir(metrics) if ('Metric' in metric) & ('Metrics' not in metric)]
metricList = [getattr(metrics, metric) for metric in metricNames]

output = []
for name, metric in zip(metricNames, metricList):
    doc = metric.__doc__
    if doc is None:
        doc = 'None'
    # Take just the first line
    doc = doc.splitlines()[0]
    output.append([name, cleandoc(doc)])

df = pd.DataFrame(output)
latex_table = df.to_latex(index=False)


