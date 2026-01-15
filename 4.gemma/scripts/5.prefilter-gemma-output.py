import argparse
import os
import glob

# Function to filter rows based on condition and write to new file
def filter_and_write(input_file, output_dir):
    output_file = os.path.join(output_dir, "filtered_" + os.path.basename(input_file))
    with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
        # Keep the first row (header)
        header = next(f_in)
        f_out.write(header)
        
        # Filter rows based on condition (final column > 0.3)
        for line in f_in:
            if float(line.strip().split()[-1]) > 0.3:
                continue  # Skip this row
            f_out.write(line)

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Filter rows based on a condition in files with *.assoc.txt extension.")
    parser.add_argument("-d", "--directory", required=True, help="Directory to search for *.assoc.txt files and output the filtered files.")
    args = parser.parse_args()

    # Get a list of all files with the extension .assoc.txt
    files = glob.glob(os.path.join(args.directory, "*.assoc.txt"))

    # Loop through each file and filter rows
    for file in files:
        filter_and_write(file, args.directory)

if __name__ == "__main__":
    main()
