import re
import os

def read_sample_pairs(file_path, fam, fam2=None):
    sample_pairs = []

    # mapping the sample names
    sample_name_to_index = {name: index for index, name in enumerate(fam.iloc[:, 1])}

    if fam2 is not None:
        # another dataset
        sample_name_to_index_2 = {name: index for index, name in enumerate(fam2.iloc[:, 1])}

    with open(file_path, 'r') as file:
        for line in file:
            sample_names = line.strip().split('\t')
            if len(sample_names) == 2:
                sample1, sample2 = sample_names[0], sample_names[1]

                if fam2 is None:
                    # within match
                    if sample1 in sample_name_to_index and sample2 in sample_name_to_index:
                        index1 = sample_name_to_index[sample1]
                        index2 = sample_name_to_index[sample2]
                        pair = (index1, index2)
                        sample_pairs.append(pair)
                else:
                    # between-match
                    if sample1 in sample_name_to_index and sample2 in sample_name_to_index_2:
                        index1 = sample_name_to_index[sample1]
                        index2 = sample_name_to_index_2[sample2]
                        pair = (index1, index2)
                        sample_pairs.append(pair)

    return sample_pairs

def validate_lengths(results, results_details):
    # save the sum lengths for each pair
    summary_lengths = {}
    for line in results:
        parts = line.split('\t')
        pair_key = (parts[0], parts[1])  # 
        length = float(parts[3])  # 
        if pair_key in summary_lengths:
            summary_lengths[pair_key] += length
        else:
            summary_lengths[pair_key] = length

    # save the IBD details for each pair
    details_lengths = {}
    for _, row in results_details.iterrows():
        pair_key = (row['family1'], row['family2'])  # 
        length = row['lengths']
        if pair_key in details_lengths:
            details_lengths[pair_key] += length
        else:
            details_lengths[pair_key] = length

    # validate summary and details
    for pair_key in summary_lengths:
        summary_total_length = summary_lengths[pair_key]
        details_total_length = details_lengths.get(pair_key, 0)
        if summary_total_length != details_total_length:
            print(f"Warning: Lengths for pair {pair_key} do not match! Summary: {summary_total_length}, Details: {details_total_length}")

    # 
    summary_total_length = sum(summary_lengths.values())
    details_total_length = sum(details_lengths.values())
    return summary_total_length, details_total_length
