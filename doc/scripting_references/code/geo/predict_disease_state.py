""" Documentation script """
from Orange.classification import LogisticRegressionLearner
from Orange.evaluation.testing import CrossValidation
from Orange.evaluation.scoring import AUC

from orangecontrib.bioinformatics.geo.dataset import GDS

gds = GDS("GDS2960")
data = gds.get_data(sample_type="disease state", transpose=True, report_genes=True)
print("Samples: %d, Genes: %d" % (len(data), len(data.domain.attributes)))

learners = [LogisticRegressionLearner()]
results = CrossValidation(data, learners, k=10)

print("AUC = %.3f" % AUC(results)[0])
