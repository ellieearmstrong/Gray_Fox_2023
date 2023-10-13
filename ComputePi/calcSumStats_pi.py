import math
import sys

#input file has three columns chrom, pos, value of pi
if len(sys.argv) != 2:
    print("Usage: python3 calcSumStats_pi.py file.pi")
    sys.exit(1)

def calculate_statistics(file_path):
    count = 0
    mean3 = M3_3 = 0
    min3 = float('inf')
    max3 = float('-inf')
    
    with open(file_path, 'rt') as f:
        # Get the column names
        header = next(f).strip().split('\t')
        col3_name = header[2]
        
        for line in f:
            data = line.strip().split('\t')
            col3 = data[2]
            print(col3)

            # Skip entries that are '.'
            if col3 == '-nan':
                continue
            
            col3 = float(col3)
            
            # Update count, mean, and variance for column 3
            count += 1
            delta3 = col3 - mean3
            mean3 += delta3 / count
            delta3_3 = col3 - mean3
            M3_3 += delta3 * delta3_3
            
            
            # Update min and max for column 3
            if col3 < min3:
                min3 = col3
            if col3 > max3:
                max3 = col3
    
    # Calculate standard deviation for each column
    stdev3 = math.sqrt(M3_3 / count)
    
    # Print the column names and calculated statistics
    print(f''Input file: {f.name}', {col3_name}: mean={mean3:.3f}, stdev={stdev3:.3f}, min={min3:.3f}, max={max3:.3f}')
    
    # Return the calculated statistics as a dictionary
    return {
        
        col3_name: {
            'mean': mean3,
            'stdev': stdev3,
            'min': min3,
            'max': max3
        }
    }

#filter missing
calculate_statistics(sys.argv[1])


