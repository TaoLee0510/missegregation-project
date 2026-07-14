# README: data/raw/salehi/

**Last Updated:** May 7, 2025

## Overview

This directory stores the raw input data required for the Salehi dataset lineage analysis pipeline. The R analysis scripts rely on the precise structure and file contents described herein.

## Directory Contents

1.  **`metadata.csv`**
    * **Description:** Contains essential metadata for all samples and their hierarchical relationships.
    * **Required Columns:**
        * `uid`: Unique identifier for each sample/entry.
        * `parent`: The `uid` of the parent sample, defining lineage connections.
        * `datasetname`: Identifier for the dataset or patient group.
        * `timepoint`: Specifies the collection timepoint for the sample.
        * `PDX_id`: Patient-Derived Xenograft identifier. This influences time delta (`dt`) calculations for specific IDs (e.g., "SA039", "SA906").
        * `library_ids`: Semicolon-separated string of library identifiers. These map metadata entries to copy number data.

2.  **`arm_loci.Rds`**
    * **Description:** An R Data Serialization (RDS) file defining chromosome arm locations. It is used to aggregate binned copy number data.
    * **Required Data Structure:** An R object (e.g., data frame) with at least the following columns:
        * `chrom`: Chromosome identifier (e.g., "1", "X").
        * `arm`: Chromosome arm identifier ("p" or "q").
        * `start`: Genomic start coordinate of the arm.

3.  **`raw_post_jump/` (Subdirectory)**
    * **Description:** Contains binned copy number (CN) data. This directory is organized into further subdirectories, each corresponding to a specific sample or library batch.
    * **Structure:**
        * Individual subdirectories (e.g., `SA501X3P1`).
            * **Note:** The subdirectory `SA004` is explicitly excluded from processing by the `extract_cn_profiles` function.
        * Within each subdirectory, a file named **`named_mat.csv`** must be present.
            * **`named_mat.csv` Description:** A CSV file containing a matrix of copy number values. Rows represent genomic bins, and columns represent individual cells or samples.
            * **`named_mat.csv` Required Columns:**
                * `bin_name`: An identifier for each genomic bin, formatted as "chromosome_startposition_..." (e.g., "1_10000_details"). This column is parsed to extract genomic coordinates.
                * Additional columns: Each subsequent column contains copy number values for a specific cell/sample corresponding to the `bin_name`.

## Usage Notes

The analysis scripts depend on this exact directory structure and file naming conventions. Ensure all specified files and directories are present and correctly formatted before executing the pipeline.
