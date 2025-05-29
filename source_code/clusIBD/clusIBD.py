import datetime
from itertools import combinations
import os
import sys
import itertools
import numpy as np
import math
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
from identify_continuous_win import identify_continuous_win
import time
import argparse
from file_operations import read_sample_pairs,validate_lengths
from genotype_analysis import compare_snps_dask, convert_genotypes
from position_and_threshold_calculation import get_position_by_group_element, binary_search_adjusted_rate, \
    relocate_position, estimate_threshold
from collections import defaultdict
from IBD_preprocess import read_and_process_files
from IBD_data_structures import (start_indices_list, end_indices_list, chr_nums_list, group_nums_list,
                                 chr_group_data, rate_array, state_array)


def process_ibd2(pair):
    future_ibd2 = executor.submit(process_pair, pair, 'IBD2')
    return future_ibd2.result()
# 20240321

parser = argparse.ArgumentParser(description='Detection of IBD segments with clusIBD using unphased genetic data.')

#parser = argparse.ArgumentParser(description='Detection of IBD segments using clusIBD')
parser.add_argument('-f', '--file_prefix', type=str, required=True, help='The prefix of bed file. This is required.')
parser.add_argument('-F', '--file_prefix_2', type=str,help='The prefix of another bed file.')
parser.add_argument('-n', '--bin_size', default=-1, type=int,
                    help='The number of SNPs per window. The default is -1.')
parser.add_argument('-q', '--minQ1', default=0.15, type=float,
                    help='The minimal percentage for IBD1. The default is 0.15.')
parser.add_argument('-Q', '--minQ2', default=0.05, type=float,
                    help='The minimal percentage for IBD2. The default is 0.05.')
parser.add_argument('-R', '--fpr', default=0.001, type=float,
                    help='The maximal false positive rate. The default is 0.001.')
parser.add_argument('-s', '--size', default=5, type=int,
                    help='The random size for parameter estimation. The default is 5.')
parser.add_argument('-L', '--IBD2_length', default=500, type=int,
                    help='The minimal length of IBD1 for IBD2 estimation. The default is 500.')
parser.add_argument('-l', '--min_length', default=7, type=int,
                    help='The minimal length for an IBD segment to be considered true. The default is 7.')
parser.add_argument('-c', '--cpu', default=5, type=int, help='The number of CPU cores to be used. The default is 5.')
parser.add_argument('-p', '--pairs_file', type=str,help='File containing sample pairs, each line represents a pair.')
parser.add_argument('-o', '--out', default='out', type=str, metavar='',
                    help='The prefix of output file.If it is not specified,"out" is used')

# Check if no arguments were passed
if len(sys.argv) == 1:
    parser.print_help(sys.stderr)  # Print the help message
    sys.stderr.write("\nPlease provide arguments. For detailed information, see the options above.\n")
    sys.exit(1)

# Parse the command line arguments
timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
args = parser.parse_args()
bin_size = args.bin_size
IBD2_length = args.IBD2_length
min_length = args.min_length
minQ1 = args.minQ1
minQ2 = args.minQ2
num_cpu = args.cpu
out_file = args.out

#if not os.path.exists(out_file):
#    os.makedirs(out_file)
#out_file = out_file+f"/{args.file_prefix.split('/')[-1]}_{timestamp}"
# call preprocess
start_time = time.time()
(bim, fam, bed, all_samples_data, bim_2, fam_2, bed_2, all_samples_data_2, call_rate_array,
 groups_idx, group_start, groups_idx_dict, pos_list, het_vec, minQ01, minQ02, Q25, t1, t2,bin_size, total_windows,
 Q50, min_num1, min_num2, start_indices_list, end_indices_list, chr_nums_list, group_nums_list,
 chr_group_data, rate_array, state_array) = read_and_process_files(args.file_prefix, args.file_prefix_2, args.bin_size, args.size, args.fpr, args.minQ1, args.minQ2)
print(f"Init processing time: {time.time() - start_time} seconds")

print("Estimated parameters:\n Q25=", Q25, ",\n Q50=", Q50, ",\n t1=", t1, ",\n t2=", t2,"\n SNPs per window = ",bin_size, ",\n Total windows = ",
       total_windows, ",\n min N for IBD1 = ", min_num1, ",\n min N for IBD2 =", min_num2)

num_rows = len(groups_idx)
results = {}

num_samples_to_load = 10
compare_two_files = False
# read sample pairs that will be processed
if args.pairs_file:
    sample_pairs = read_sample_pairs(args.pairs_file, fam, fam2=fam_2)
    print(f"Successfully read {len(sample_pairs)} sample pairs from {args.pairs_file}")
else:
    if args.file_prefix_2:
        compare_two_files = True      ##byliran 20240829
        sample_pairs = list(itertools.product(range(len(fam)), range(len(fam_2))))
    else:
        sample_pairs = list(itertools.combinations(range(len(fam)), 2))
# if not args.pairs_file:
#     print("no samples")
invalid_pairs = []

##function for detecting IBD segments for each sample pair
def process_pair(sample_indices, type_id):
    if compare_two_files:
        family1 = fam.iloc[sample_indices[0]]['iid']
        family2 = fam_2.iloc[sample_indices[1]]['iid']
    else:
        family1 = fam.iloc[sample_indices[0]]['iid']
        family2 = fam.iloc[sample_indices[1]]['iid']

    try:
        if compare_two_files:
            snp_data_sample1 = all_samples_data[:, sample_indices[0]]
            snp_data_sample2 = all_samples_data_2[:, sample_indices[1]]
        else:
            snp_data_sample1 = all_samples_data[:, sample_indices[0]]
            snp_data_sample2 = all_samples_data[:, sample_indices[1]]
        combined_results = compare_snps_dask(snp_data_sample1, snp_data_sample2, type_id)
        # calculate rate
        for i in range(num_rows):
            start = start_indices_list[i]
            end = end_indices_list[i]
            chr_num = chr_nums_list[i]
            group = group_nums_list[i]
            group_data = combined_results[start:end+1]
            chr_group_data[chr_num][group] = group_data
            effective_SNPnum = len(group_data) * call_rate_array[sample_indices[0]] * call_rate_array[sample_indices[1]]
            true_values = np.count_nonzero(group_data)
            rate_array[i] = true_values / effective_SNPnum if effective_SNPnum > 0 else np.nan
        threshold = estimate_threshold(rate_array, Q25 ** 2 / 2, Q50, minQ1 if type_id == "IBD1" else minQ2,
                                       t1 if type_id == "IBD1" else t2, 100)
        #threshold = t1 if type_id == "IBD1" else t2
        state_array = rate_array < threshold
        '''
        #save the rates
        oph_rate = np.append(rate_array,threshold)
        #if type_id == "IBD2": np.savetxt(f"./results/clusIBD_0_{sample_indices[0]}_{sample_indices[1]}.rate",oph_rate)
        if sample_indices[0] == 13 and sample_indices[1] == 29:
            np.savetxt(f"./results/clusIBD_0_{sample_indices[0]}_{sample_indices[1]}.rate",oph_rate)
            print(f"threshold is {threshold}")
        '''
        ibd_results_list = []
        #min_num1 = 6
        for chr_num in range(1, 23):
            start_idx = group_start[chr_num - 1]
            end_idx = group_start[chr_num]
            state_data = state_array[start_idx:end_idx]
            chr_length = len(state_data)
            idx_list = identify_continuous_win(state_data, min_num1 if type_id == "IBD1" else min_num2,
                                               min_num1 // 2 if type_id == "IBD1" else min_num2 // 2, 1)

            if len(idx_list) == 0:
                continue

            for window_num, idx in enumerate(idx_list, start=1):
                if idx[0] > 1:
                    previous_group = idx[0] - 1 + 1
                    current_group = idx[0] + 1
                    try:
                        low = binary_search_adjusted_rate(chr_group_data, chr_num, current_group, previous_group,
                                                          1)
                        start_idx = groups_idx_dict[chr_num][previous_group]['start_idx']
                        start_pos = pos_list[start_idx + low - 1]
                    except (IndexError, KeyError) as e:
                        print(f"Error in binary_search_adjusted_rate or accessing pos_list for start position: {e}")
                        start_pos = pos_list[start_idx]
                else:
                    start_idx = groups_idx_dict[chr_num][1]['start_idx']
                    start_pos = pos_list[start_idx]

                if idx[-1] < chr_length - 1:
                    previous_group = idx[-1] + 1
                    current_group = idx[-1] + 1 + 1
                    try:
                        high = binary_search_adjusted_rate(chr_group_data, chr_num, current_group, previous_group,
                                                           2)
                        end_idx = groups_idx_dict[chr_num][previous_group]['start_idx'] + math.floor(bin_size / 2)
                        end_pos = pos_list[end_idx + high - 1]
                    except (IndexError, KeyError) as e:
                        end_pos = pos_list[end_idx]
                else:
                    end_idx = groups_idx_dict[chr_num][chr_length]['start_idx']
                    end_pos = pos_list[end_idx - 1]


                if not pd.isna(start_pos):
                #exclude short segments
                    if (end_pos - start_pos) >= min_length:
                        ibd_results_list.append(
                        {'family1': family1, 'family2': family2, 'chr': chr_num, 'start_pos': start_pos,
                         'end_pos': end_pos, 'lengths': end_pos - start_pos})

        num_ibd_segments = len(ibd_results_list)
        total_ibd_length = sum(row['lengths'] for row in ibd_results_list)
        return {
            'family1': family1,
            'family2': family2,
            'num_ibd_segments': num_ibd_segments,
            'total_ibd_length': total_ibd_length,
            'threshold': threshold,
            'detail': ibd_results_list
        }
        return ibd_df
    except Exception as e:
        print(f"Error processing pair {sample_indices}: {e}")
        return None

start_time = time.time()
out_summary = pd.DataFrame(index=range(len(sample_pairs)), columns=['sample1', 'sample2', 'number_segments', 'lengths'])
n_finished = 0
p_finished = 0

results = []
results_details = []

#from collections import defaultdict

start_time = time.time()

ibd2_pairs = []
summary_dict = defaultdict(lambda: {'num_ibd_segments': 0, 'total_ibd_length': 0})

#detect IBD1 segments
with ProcessPoolExecutor(num_cpu) as executor:
    future_to_pair = {executor.submit(process_pair, pair, 'IBD1'): pair for pair in sample_pairs}
    for future in tqdm(as_completed(future_to_pair), total=len(sample_pairs)):
        pair = future_to_pair[future]
        try:
            result = future.result()

            if result is None:
                print(f"No result returned for pair: {pair}")
                continue

            if 'detail' not in result:
                print(f"No 'detail' in result for pair: {pair}")
                continue

            #  'detail' valid results
            for detail in result['detail']:
                detail['type'] = 'IBD1'
            results_details.extend(result['detail'])

            family_pair = (result['family1'], result['family2'])

            # refresh summary_dict to IBD1
            summary_dict[family_pair]['num_ibd_segments'] += result['num_ibd_segments']
            summary_dict[family_pair]['total_ibd_length'] += result['total_ibd_length']

            # if IBD will be analyzed, do not refresh summary
            if result['total_ibd_length'] > IBD2_length:
                ibd2_pairs.append(pair)  #save the sampe pairs that will be processed for IBD2 detection

        except Exception as e:
            print(f"Error when processing pair {pair}: {e}")

ibd1_end_time = time.time()
print(f"IBD1 process time: {ibd1_end_time - start_time} seconds")

# detect IBD2 segments
with ProcessPoolExecutor(num_cpu) as executor:
    future_to_pair = {executor.submit(process_pair, pair, 'IBD2'): pair for pair in ibd2_pairs}
    for future in tqdm(as_completed(future_to_pair), total=len(ibd2_pairs)):
        pair = future_to_pair[future]
        try:
            result_ibd2 = future.result()
            if result_ibd2 is None:
                continue

            family_pair = (result_ibd2['family1'], result_ibd2['family2'])

            # refresh summary_dict to IBD2
            summary_dict[family_pair]['num_ibd_segments'] += result_ibd2['num_ibd_segments']
            summary_dict[family_pair]['total_ibd_length'] += result_ibd2['total_ibd_length']

            # add the IBD2 results to details
            for detail in result_ibd2['detail']:
                detail['type'] = 'IBD2'
            results_details.extend(result_ibd2['detail'])

        except Exception as e:
            print(f"Error when processing pair {pair}: {e}")

ibd2_end_time = time.time()
print(f"IBD2 process time: {ibd2_end_time - ibd1_end_time} seconds")
print(f"Total process time: {ibd2_end_time - start_time} seconds")

# save summary
results = [f"{family1}\t{family2}\t{data['num_ibd_segments']}\t{data['total_ibd_length']}\n"
           for (family1, family2), data in summary_dict.items()]

with open(f"{out_file}.IBD.summary", 'w') as file:
    file.writelines(results)

# save details
df_details = pd.DataFrame(results_details)
df_details.to_csv(f"{out_file}.IBD.details", sep='\t', index=False, header=False)

print(f"Summary file generated: {out_file}.IBD.summary")
print(f"Details file generated: {out_file}.IBD.details")
