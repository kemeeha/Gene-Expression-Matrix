# Import libraries
import os
import pandas as pd
import numpy as np
import argparse

def parse_gene_expression_files(input_dir, output_file, gene_column_index=1, 
                               expression_column_index=5, skip_rows=6):
    """
    Parse gene expression TSV files and create a combined gene expression matrix.
    
    Parameters:
    -----------
    input_dir : str
        Directory containing the gene expression files or subdirectories with files
    output_file : str
        Path to save the output gene expression matrix
    gene_column_index : int
        Index of the column containing gene names (0-based, default: 1 for column 2)
    expression_column_index : int
        Index of the column containing expression values (0-based, default: 5 for column 6)
    skip_rows : int
        Number of header rows to skip (default: 6)
    """
    print(f"Searching for TSV files in: {input_dir}")
    
    # Dictionary to store gene expression data
    gene_data = {}
    sample_ids = []
    
    # Find all TSV files in the directory and subdirectories
    tsv_files = []
    for root, dirs, files in os.walk(input_dir):
        for file in files:
            if file.endswith('.tsv'):
                tsv_files.append(os.path.join(root, file))
    
    print(f"Found {len(tsv_files)} TSV files")
    
    # Process each TSV file
    for file_path in tsv_files:
        # Extract sample ID from filename
        sample_id = os.path.basename(file_path).split(".tsv")[0]
        if sample_id not in sample_ids:
            sample_ids.append(sample_id)
        
        print(f"Processing file: {file_path} (Sample ID: {sample_id})")
        
        try:
            # Read the file line by line, skipping header
            with open(file_path, 'r') as file:
                # Skip first N rows as specified
                for _ in range(skip_rows):
                    file.readline()
                
                # Process each line
                for line in file:
                    parts = line.strip().split('\t')
                    
                    # Skip malformed lines
                    if len(parts) <= max(gene_column_index, expression_column_index):
                        continue
                    
                    gene_name = parts[gene_column_index]
                    expression_value = parts[expression_column_index]
                    
                    # Initialize gene entry if not exists
                    if gene_name not in gene_data:
                        gene_data[gene_name] = {}
                    
                    # Store expression value
                    gene_data[gene_name][sample_id] = expression_value
        
        except Exception as e:
            print(f"Error processing file {file_path}: {str(e)}")
    
    print(f"Processed {len(sample_ids)} samples with {len(gene_data)} genes")
    
    # Create DataFrame from the collected data
    df_output = pd.DataFrame(index=gene_data.keys(), columns=sample_ids)
    
    # Fill the DataFrame with expression values
    for gene_name, samples in gene_data.items():
        for sample_id, value in samples.items():
            df_output.at[gene_name, sample_id] = value
    
    # Replace missing values with NA
    df_output = df_output.fillna("NA")
    
    # Set index name
    df_output.index.name = "Gene"
    
    # Save to file
    df_output.to_csv(output_file, sep='\t')
    print(f"Gene expression matrix saved to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Parse gene expression TSV files and create a combined expression matrix.')
    parser.add_argument('--input_dir', default='gdc-blca', help='Directory containing gene expression files')
    parser.add_argument('--output_file', default='gene_expression_matrix.tsv', help='Path to save the output gene expression matrix')
    parser.add_argument('--gene_column_index', type=int, default=1, help='Index of column for gene identifiers (0-based)')
    parser.add_argument('--expression_column_index', type=int, default=5, help='Index of column for expression values (0-based)')
    parser.add_argument('--skip_rows', type=int, default=6, help='Number of header rows to skip')
    
    args = parser.parse_args()
    
    parse_gene_expression_files(
        args.input_dir,
        args.output_file,
        args.gene_column_index,
        args.expression_column_index,
        args.skip_rows
    )
