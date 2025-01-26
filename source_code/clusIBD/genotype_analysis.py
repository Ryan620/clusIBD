import numpy as np

def convert_genotypes(block, genotype_mapping):
    converted = [genotype_mapping.get(item, "NA") if not np.isnan(item) else "NA" for item in block.ravel()]
    return np.array(converted).reshape(block.shape)
def compare_snps_dask(g1, g2, types):
    if types == "IBD1":
        return (g1 - g2) ** 2 == 4
    else:
        valid_comparison = ~(np.isnan(g1) | np.isnan(g2))
        comparison_result = (g1 != g2) & valid_comparison
        return comparison_result