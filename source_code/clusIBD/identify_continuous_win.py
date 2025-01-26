import numpy as np
def identify_continuous_win(x, min_num=5, gap_between=5, gap_flank=1):
    def find_continuous_segments(data, min_len):
        segments = []
        start = 0
        for i in range(1, len(data)):
            if data[i] - data[i - 1] != 1:
                if i - start >= min_len:
                    segments.append(data[start:i])
                start = i
        if len(data) - start >= min_len:
            segments.append(data[start:])
        return segments

    while True:
        changed = False
        idx0 = np.where(x != 0)[0]
        if len(idx0) <= 1:
            return []
        idx_list = find_continuous_segments(idx0, min_num)
        new_idx_list = []
        i = 0
        while i < len(idx_list):
            if i < len(idx_list) - 1 and idx_list[i + 1][0] - idx_list[i][-1] <= gap_between + 1:
                x[idx_list[i][-1]:idx_list[i + 1][0] + 1] = 1
                new_segment = np.concatenate((idx_list[i], idx_list[i + 1]))
                new_idx_list.append(new_segment)
                i += 1
                changed = True
            else:
                new_idx_list.append(idx_list[i])
                i += 1

        idx_list = new_idx_list
        if changed:
            continue
        # ...
            # Handling gap_flank
        for idxes in idx_list:
            start_idx = max(idxes[0] - gap_flank - 1, 0)
            end_idx = min(idxes[-1] + gap_flank + 1, len(x))
            x[start_idx:end_idx] = np.where(np.isnan(x[start_idx:end_idx]), True, x[start_idx:end_idx])
            if np.any(np.isnan(x[start_idx:end_idx])):
                x[start_idx:end_idx] = np.where(np.isnan(x[start_idx:end_idx]), True, x[start_idx:end_idx])
                changed = True

        if not changed:
            break

    return idx_list