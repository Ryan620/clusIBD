import time
import pandas_plink
import numpy as np
import pandas as pd
import random
import math
from IBD_data_structures import initialize_data_structures
def read_and_process_files(file_prefix, file_prefix_2=None, bin_size=-1, random_size=10, fpr=0.001, minQ1=0.15, minQ2=0.05):

    bim, fam, bed = pandas_plink.read_plink(file_prefix)
    all_samples_data = bed[:, :len(fam)].compute()
    if (bim.iloc[:, 2] != 0).all():
        bim['pos'] = bim.iloc[:, 2]
    else:
        bim['pos'] = bim.iloc[:, 3] / 1000000
    if file_prefix_2:
        bim_2, fam_2, bed_2 = pandas_plink.read_plink(file_prefix_2)
        all_samples_data_2 = bed_2[:, :len(fam_2)].compute()
    else:
        bim_2, fam_2, bed_2, all_samples_data_2 = None, None, None, None

    call_rate_array = np.ones(fam.shape[0])
    for i in range(fam.shape[0]):
        genotypes = all_samples_data[:, i]
        nan_count = len(np.where(np.isnan(genotypes))[0])
        call_rate_array[i] = 1 - nan_count / bim.shape[0]

    # init_bin
    #chromosomes = bim['chrom'].unique()
    shifted_pos = bim['pos'].shift(-1)
    condition = bim['pos'] > shifted_pos
    condition = condition.fillna(False)
    indices = condition[condition].index
    x = bim.loc[indices, 'pos']
    genome_len = x.sum()
    ##to identify IBD segments over 10 
    if bin_size < 0: bin_size = math.floor(9 * bim.shape[0] / (genome_len * 8))
    if bin_size < 150: bin_size = 150


    groups_idx = pd.DataFrame(columns=['chr', 'groups', 'start_idx', 'end_idx'])
    total_n = 0
    group_start = [0] * 23
    for i in range(1, 23):
        chrom_data = bim[bim['chrom'] == str(i)]
        n = len(chrom_data)
        m, m2 = divmod(n, bin_size)
        chr_groups = np.arange(1, m + (1 if m2 > 0 else 0) + 1)
        start_indices = total_n + np.arange(0, len(chr_groups) * bin_size, bin_size)
        end_indices = np.append(start_indices[1:], total_n + n) - 1
        chrom = i
        chr_df = pd.DataFrame({
            'chr': chrom,
            'groups': chr_groups,
            'start_idx': start_indices,
            'end_idx': end_indices
        })
        groups_idx = pd.concat([groups_idx, chr_df], ignore_index=True)
        group_start[i] = len(chr_groups)
        total_n += n

    for i in range(1, 23):
        group_start[i] += group_start[i - 1]

    groups_idx_dict = {}
    for _, row in groups_idx.iterrows():
        chr_num = int(row['chr'])
        group_num = row['groups']
        start_idx = row['start_idx']
        end_idx = row['end_idx']
        if chr_num not in groups_idx_dict:
            groups_idx_dict[chr_num] = {}
        groups_idx_dict[chr_num][group_num] = {'start_idx': start_idx, 'end_idx': end_idx}

    pos_list = bim['pos'].tolist()
    start_indices_list, end_indices_list, chr_nums_list, group_nums_list, chr_group_data, rate_array, state_array = initialize_data_structures(groups_idx, groups_idx_dict)
    ##set the seed so that consistent rescults can be obtained for every sampling
    random.seed(random_size)
    num_samples = min(random_size, fam.shape[0])
    sample_random = random.sample(range(fam.shape[0]), num_samples)
    het_array = np.zeros((groups_idx.shape[0], len(sample_random)))

    snp_data_samples = all_samples_data[:, sample_random]
    combined_ibs_states_samples = snp_data_samples == 1.0  # homozygote

    for j in range(groups_idx.shape[0]):
        idx1, idx2 = groups_idx.iloc[j, 2:4]
        states = combined_ibs_states_samples[idx1:idx2, :]
        het_array[j, :] = np.where(np.count_nonzero(~np.isnan(states), axis=0), np.nanmean(states, axis=0), 0)

    het_vec = het_array.flatten()
    het_vec = het_vec[het_vec > 0]

    minQ01 = np.quantile(het_vec, q=minQ1, interpolation='linear')
    minQ02 = np.quantile(het_vec, q=minQ2, interpolation='linear')
    Q25 = np.quantile(het_vec, q=0.25, interpolation='linear')
    t1 = minQ01 ** 2 / 2
    t2 = minQ02
    total_windows = len(groups_idx)
    Q50 = np.quantile(het_vec, q=0.50, interpolation='linear')
    min_num1 = math.ceil(math.log10(fpr / total_windows) / math.log10(minQ1))
    min_num2 = math.ceil(math.log10(fpr / total_windows) / math.log10(minQ2))

    return (bim, fam, bed, all_samples_data, bim_2, fam_2, bed_2, all_samples_data_2, call_rate_array,
            groups_idx, group_start, groups_idx_dict, pos_list, het_vec, minQ01, minQ02, Q25, t1, t2,bin_size, total_windows,
            Q50, min_num1, min_num2, start_indices_list, end_indices_list, chr_nums_list, group_nums_list,
            chr_group_data, rate_array, state_array)
