# Moore’s Rayleigh Analysis Python Script

A Python script for calculating Moore’s Modified Rayleigh test and generating publication-style circular plots of directional vector data.

The script reads input data from an Excel file, calculates the mean direction, rank-weighted vector length, 95% bootstrap confidence interval, Monte Carlo p-value, and creates a circular plot with individual vectors and statistical results.

The generated plots are automatically saved in three formats: PNG, PDF, and SVG.

## Repository Files

This repository contains:

* Moores_Rayleigh_analysis_Python_script.py — the main Python script;
* Data_format_example.xlsx — an example Excel file showing the required input data format.

## Input Data

The repository includes the example Excel file Data_format_example.xlsx, which shows how the input data should be organized.

The script uses the first two columns of the Excel file:

* the first column must contain directions in degrees;
* the second column must contain vector length or orientation strength.

Column names can be different, because the script always reads the first two columns of the Excel file.

## Requirements

The script requires Python 3 and the following Python packages:

pip install numpy pandas matplotlib openpyxl

## Usage

Run the script:

python Moores_Rayleigh_analysis_Python_script.py

After running the script, a file selection window will open. Select the Excel file with your data.

You can also provide the Excel file path directly:

python Moores_Rayleigh_analysis_Python_script.py data.xlsx

## Output

After processing the input data, the script automatically creates three output files:

* PNG image;
* PDF file;
* SVG vector file.

The output files are saved next to the input Excel file unless another output path is specified.
