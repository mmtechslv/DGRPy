import numpy as np
from biom.table import Table
data = np.arange(40).reshape(10, 4)
sample_ids = ['S%d' % i for i in range(4)]
observ_ids = ['O%d' % i for i in range(10)]
sample_metadata = [{'environment': 'A'}, {'environment': 'B'},
{'environment': 'A'}, {'environment': 'B'}]
observ_metadata = [{'taxonomy': ['Bacteria', 'Firmicutes']},{'taxonomy': ['Bacteria', 'Firmicutes']},{'taxonomy': ['Bacteria', 'Proteobacteria']},
                    {'taxonomy': ['Bacteria', 'Proteobacteria']},
                    {'taxonomy': ['Bacteria', 'Proteobacteria']},
                    {'taxonomy': ['Bacteria', 'Bacteroidetes']},
                    {'taxonomy': ['Bacteria', 'Bacteroidetes']},
                    {'taxonomy': ['Bacteria', 'Firmicutes']},
                    {'taxonomy': ['Bacteria', 'Firmicutes']},
                    {'taxonomy': ['Bacteria', 'Firmicutes']}]
table = Table(data, observ_ids, sample_ids, observ_metadata,sample_metadata, table_id='Example Table')