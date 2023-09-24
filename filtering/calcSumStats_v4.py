import math

def calculate_statistics(file_path):
    count = 0
    mean1 = mean2 = M2_1 = M2_2 = 0
    min1 = min2 = float('inf')
    max1 = max2 = float('-inf')
    
    with open(file_path, 'r') as f:
        # Get the column names
        header = next(f).strip().split('\t')
        col1_name, col2_name = header[0], header[1]
        
        for line in f:
            col1, col2 = line.strip().split('\t')
            
            # Skip entries that are '.'
            if col1 == '.' or col2 == '.':
                continue
            
            col1, col2 = float(col1), float(col2)
            
            # Update count, mean, and variance for column 1
            count += 1
            delta1 = col1 - mean1
            mean1 += delta1 / count
            delta2_1 = col1 - mean1
            M2_1 += delta1 * delta2_1
            
            # Update count, mean, and variance for column 2
            delta2 = col2 - mean2
            mean2 += delta2 / count
            delta2_2 = col2 - mean2
            M2_2 += delta2 * delta2_2
            
            # Update min and max for column 1
            if col1 < min1:
                min1 = col1
            if col1 > max1:
                max1 = col1
            
            # Update min and max for column 2
            if col2 < min2:
                min2 = col2
            if col2 > max2:
                max2 = col2
    
    # Calculate standard deviation for each column
    stdev1 = math.sqrt(M2_1 / count)
    stdev2 = math.sqrt(M2_2 / count)
    
    # Print the column names and calculated statistics
    print(f'Input file: {f.name}')
    print(f'{col1_name}: mean={mean1:.3f}, stdev={stdev1:.3f}, min={min1:.3f}, max={max1:.3f}')
    print(f'{col2_name}: mean={mean2:.3f}, stdev={stdev2:.3f}, min={min2:.3f}, max={max2:.3f}')
    
    # Return the calculated statistics as a dictionary
    return {
        col1_name: {
            'mean': mean1,
            'stdev': stdev1,
            'min': min1,
            'max': max1
        },
        col2_name: {
            'mean': mean2,
            'stdev': stdev2,
            'min': min2,
            'max': max2
        }
    }

#filter missing
calculate_statistics('filtered_canFam3.1_gvcfstats.DP_NMISS.out')
calculate_statistics('filtered_arcticfox_gvcfstats.DP_NMISS.out')
calculate_statistics('filtered_canFam4_gvcfstats.DP_NMISS.out')
calculate_statistics('filtered_grayfox_gvcfstats.DP_NMISS.out')

#unfiltered
#calculate_statistics('canFam3.1_gvcfstats.DP_NMISS.out')
#calculate_statistics('arcticfox_gvcfstats.DP_NMISS.out')
#calculate_statistics('canFam4_gvcfstats.DP_NMISS.out')
#calculate_statistics('grayfox_gvcfstats.DP_NMISS.out')

