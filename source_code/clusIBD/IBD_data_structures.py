import numpy as np

start_indices_list, end_indices_list, chr_nums_list, group_nums_list = None, None, None, None
chr_group_data, rate_array, state_array = None, None, None


def initialize_data_structures(groups_idx, groups_idx_dict):
    # global start_indices_list, end_indices_list, chr_nums_list, group_nums_list
    # global chr_group_data, rate_array, state_array

    start_indices_list = groups_idx['start_idx'].tolist()
    end_indices_list = groups_idx['end_idx'].tolist()
    chr_nums_list = groups_idx['chr'].tolist()
    group_nums_list = groups_idx['groups'].tolist()

    chr_group_data = {chr_num: {group: [] for group in groups_idx_dict[chr_num]} for chr_num in groups_idx_dict}
    rate_array = np.zeros(len(start_indices_list))
    state_array = np.zeros(len(start_indices_list), dtype=bool)
    return start_indices_list, end_indices_list, chr_nums_list, group_nums_list, chr_group_data, rate_array, state_array
