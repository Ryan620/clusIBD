import numpy as np
from scipy.stats.mstats import mquantiles
from scipy.stats import gaussian_kde


def estimate_threshold(x, Q25, het, min_percent=0.1, min_threshold=0.02, min_Rate_num=100):
    x = np.array(x)
    x = x[~np.isnan(x)]  # Remove NaN values
    if len(x) < min_Rate_num:
        return max(np.max([min_threshold, mquantiles(x, prob=min_percent)]))

    x1 = np.min(x)
    x2 = mquantiles(x, prob=0.95)
    delt_x = (x2 - x1) / 100
    if delt_x == 0:
        return min_threshold
    x_ranges = np.arange(start=x1, stop=x2, step=delt_x)
    kde = gaussian_kde(x, bw_method="silverman")
    probs = kde(x_ranges)
    probs_dif = np.diff(probs)

    lowest_pos = np.where(
        (probs_dif[:-5] < 0) & (probs_dif[1:-4] < 0) & (probs_dif[2:-3] < 0) & (probs_dif[3:-2] < 0) & (
                    probs_dif[4:-1] < 0) & (probs_dif[5:] > 0))[0] + 5
    highest_pos = np.where(
        (probs_dif[:-5] > 0) & (probs_dif[1:-4] > 0) & (probs_dif[2:-3] > 0) & (probs_dif[3:-2] > 0) & (
                    probs_dif[4:-1] > 0) & (probs_dif[5:] < 0))[0] + 5

    if len(lowest_pos) > 0:
        #if there are two peaks and bottom prob is high than 0.9 times the peak prob, it is hard to differentiate them
        #if the lowest prob is higher than the het value, it is the bottom between IBD0 and IBD1 for IBD2 estimate, not the expected 
        #between IBD2 and IBD1. 
        if np.max(x_ranges[lowest_pos]) > het or np.max(probs[lowest_pos]) > np.max(probs) * 0.9:
            lowest_pos = []

    t1 = x_ranges[lowest_pos] if len(lowest_pos) > 0 else np.array([])
    '''
    if len(lowest_pos) == 0 and len(highest_pos) > 0:
        if np.mean(x) < Q25:
            num = np.count_nonzero(x > x_ranges[highest_pos[0]])
            if num / len(x) < 2.5 / 4:  # PC
                t1 = mquantiles(x, prob=0.99)
            else:  # FS
                t1 = mquantiles(x, prob=0.75 + 0.25 * min_percent)
    '''
    t2 = mquantiles(x, prob=min_percent)
    if len(t1) > 0:
        t1 = t1[0]
    threshold = np.max([t1, t2, min_threshold])

    return threshold


def relocate_position(x, threshold):
    new_x = np.where(~x)[0]
    for i, j in enumerate(new_x):
        rate = (i + 1) / (j + 1)
        if rate > threshold:
            return j - 1




def binary_search_adjusted_rate(chr_group_data, chr_num, current_group, previous_group, type_id):
    current_data = chr_group_data[chr_num][current_group]
    previous_data = chr_group_data[chr_num][previous_group]

    if type_id == 1:
        combined_data = np.concatenate((previous_data, current_data[:len(current_data) // 2]))
    elif type_id == 2:
        combined_data = np.concatenate((previous_data[len(previous_data) // 2:], current_data))

    if type_id == 1:
        for i in range(len(combined_data)-1, -1, -1):
            if combined_data[i]:
                return i
    elif type_id == 2:
        for i in range(len(combined_data)):
            if combined_data[i]:
                return i

    return 0


def get_position_by_group_element(bim, groups_idx, chr_num, group_num, element_index):
    group_info = groups_idx[(groups_idx['chr'] == chr_num) & (groups_idx['groups'] == group_num)]
    if not group_info.empty:
        start_idx = group_info.iloc[0]['start_idx']
        if start_idx + element_index < len(bim):
            return bim.at[start_idx + element_index, 'pos']
    return None
