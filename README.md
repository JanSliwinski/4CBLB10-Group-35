# Project: Data Evaluation and Visualization

## Overview
This project processes datasets stored in the `AdjustedData` folder to evaluate and visualize data. The workflow includes loading the data, performing calculations, and creating visualizations. Key scripts and folder structures are outlined below for ease of use and configuration.

---

## Folder Structure
```
ProjectRoot/
|
|-- ExampleDataSet/
|   |-- AdjustedData/          # Contains renamed experiment data TXT files
|   |-- CompiledEmissions.csv  # Additional CSV file for supplementary data
|
|-- MainScript.m               # Main script for data loading and analysis
|-- DataLoading.m              # Handles loading of TXT and CSV data
|-- PlottingScript.m           # Dedicated script for detailed plotting
|-- RedundantCode/             # Contains unused or archived project code
```

---

## Usage

### Data Loading
The dataset is loaded using the `DataLoading.m` script. Ensure the `ExampleDataSet` folder is on your current MATLAB path or update the file paths in the script to match your local machine:

```matlab
% Specify the folder containing the renamed experiment data TXT files
dataFolder = 'ExampleDataSet/AdjustedData'; % <-- Replace with your actual folder path

% Specify the path to the additional data CSV file
additionalCSVPath = 'ExampleDataSet/CompiledEmissions.csv'; % <-- Replace with your actual CSV file path
```

If the `ExampleDataSet` folder is not on your MATLAB path, you must:
1. Add the folder to your MATLAB path, or
2. Replace the paths in the code with the absolute paths on your machine.

### Main Workflow
The **`MainScript.m`** file serves as the central script to:
- Call most of the functions
- Perform the majority of data analysis and calculations

### Advanced Plotting
The **`PlottingScript.m`** is a separate script dedicated to creating precise and detailed visualizations. Use this script if the visualizations from the main script require further refinement or customization.

### Redundant Code
The **`RedundantCode/`** folder contains scripts that are not currently used in the project. These may include experimental or obsolete code for reference.

---

## Requirements
- MATLAB R2020a or later
- The dataset folder `ExampleDataSet/` should be organized as described above.

---

## Getting Started
1. Clone the repository or download the project files.
2. Add the project root folder to your MATLAB path.
3. Run `MainScript.m` to load and analyze the data.
4. For detailed plotting, execute `PlottingScript.m`.

---

## Notes
- Modify paths in `DataLoading.m` to reflect the location of your dataset if different from the default structure.
- Ensure all required files (e.g., TXT and CSV) are present in the specified directories before running the scripts.

---