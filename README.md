# Coordination Analysis

A Fortran program for analyzing atomic coordination numbers from molecular dynamics trajectories with high performance.

## Overview

This program efficiently calculates coordination numbers between different atom types in molecular dynamics simulations using an optimized cell-based spatial partitioning approach. It's designed for large-scale simulations with millions of atoms and supports parallel execution via OpenMP for excellent performance.

## Key Features

- **Fast Coordination Analysis**: Uses cell lists for O(N) performance instead of the naive O(N²) approach
- **Multi-core Processing**: Parallel processing with OpenMP for faster calculations
- **Flexible Input Options**: 
  - Supports LAMMPS trajectory files (.dump, .dmp) and data files
  - Dynamic box dimensions support (can read changing box sizes from trajectory)
- **Advanced Selection Capabilities**: 
  - Analyze specific atom types or ranges (e.g., "1-100,200,300-400")
  - Include/exclude specific atom types
- **Detailed Output**: 
  - Coordination numbers between all atom types
  - Optional neighbor tracking for detailed analysis
  - Parsable output format for post-processing
- **Progress Tracking**: Real-time progress display with estimated completion time

## Installation and Compilation

### Compilation

```bash
# Clone the repository
git clone repository
cd coordination-analysis

# Compile optimized version
make release

# Compile debug version
make debug

# Compile benchmark version
make benchmark
```

## Usage

1. Create a `setup.txt` configuration file (see [Configuration](#configuration) section)
2. Place your LAMMPS data file and trajectory file in the working directory
3. Run the analysis:

```bash
./coordination_analysis
```

## Configuration

The program is configured through a `setup.txt` file with the following sections:

### System Configuration

```
# System settings
Cores=4                 # Number of OpenMP threads to use
ATOM_TYPES=3            # Number of atom types

# Define atom types
TYPE_1=Li               # Type names
TYPE_2=Li  
TYPE_3=O

# Define cutoff distances for each pair
PAIR_CUTOFF_1_1=3.0     # Cutoffs in Angstroms
PAIR_CUTOFF_1_2=3.2
PAIR_CUTOFF_1_3=3.1
PAIR_CUTOFF_2_2=3.3
PAIR_CUTOFF_2_3=3.2
PAIR_CUTOFF_3_3=3.4
```

### Input/Output Configuration

```
# Input/output files
INPUT_DATA=lammps.data
INPUT_TRAJECTORY=trajectory.dmp
OUTPUT_FILE=coordination_numbers.dat

# Box dimensions source (optional)
USE_TRAJECTORY_BOX=yes  # Read box dimensions from trajectory instead of data file
```

### Analysis Options

```
# Frame selection
START_FRAME=0           # Starting frame (0 = beginning)
END_FRAME=all           # Ending frame (all = process all frames)

# Atom selection
SELECTED_ATOMS=all      # Optional: atom ID selection (e.g., "1-100,200,300-400")
SHOW_NEIGHBORS=yes      # Track neighbor atom IDs (required for visualization)

# Cell list options
CELL_UPDATE_FREQ=20     # How often to update cell lists
SELECTIVE_REBUILD=yes   # Use selective rebuilding for better performance
```

### Custom Trajectory Parsing

```
# Custom trajectory format (only needed for non-standard formats)
TRAJECTORY_HEADER=ITEM: ATOMS id type x y z
```

## Understanding Box Dimensions

The program offers two modes for handling simulation box dimensions:

1. **Data File Mode (default)**: Reads box dimensions from the LAMMPS data file and keeps them fixed throughout the analysis
2. **Trajectory Mode**: Reads box dimensions from each frame of the trajectory file, allowing for variable box sizes

Set the mode using the `USE_TRAJECTORY_BOX` option in `setup.txt`:

```
USE_TRAJECTORY_BOX=yes  # Use dimensions from trajectory file
```

This is especially useful for:
- NPT simulations where the box size changes
- Deformation simulations
- Sequential analysis of multiple different systems

When using trajectory box dimensions:
- Cell lists are automatically rebuilt when box dimensions change
- Performance statistics track these rebuilds properly
- All coordination calculations use the current frame's box dimensions

## Output Format

The program generates a coordination analysis file with the following format:

```
# Coordination Number Analysis
# Column format:
#   1. Atom ID
#   2. Atom Type
#   3. CN(Type 1-1)
#   4. CN(Type 1-2)
#   5. CN(Type 1-3)
#   ...
# Frame: 1
     1             Li      5(2,6754,18802,6754,18802)      0      4(17,12065,6750,6750)      0      0      0
     2             Li      6(1,33,12081,3,12051,25)      0      5(17,12065,30,18,30)      0      0      0
```

Each line shows:
- Atom ID
- Atom type/element name
- Series of coordination numbers for each pair type
- Optional neighbor atom IDs in parentheses (when `SHOW_NEIGHBORS=yes`)

The output format is designed to be easily parsed by other programs:
- Fixed-width columns for atom ID and type
- Consistent formatting for coordination numbers
- Neighbor IDs in compact parenthesized lists with no extraneous spaces

## Performance Considerations

### Cell List Optimization

The program uses a cell-based spatial partitioning system for efficient neighbor finding:

- **Cell Size**: Automatically optimized based on cutoff distances
- **Selective Rebuilding**: Only rebuilds cells that need updating
- **Update Frequency**: Configurable with `CELL_UPDATE_FREQ`

### Parallel Processing

OpenMP parallelization provides significant speedup on multicore systems:

- Set desired thread count with `Cores=N` in setup.txt
- Default uses all available cores if not specified
- Best performance typically at 4-16 cores depending on system size

## Advanced Usage

### Selective Atom Analysis

Analyze only specific atoms with the `SELECTED_ATOMS` option:

```
SELECTED_ATOMS=1-100,200,300-400
```

This will only process atoms with IDs 1-100, 200, and 300-400, significantly improving performance for large systems when you're only interested in specific regions.

### Custom Trajectory Formats

For non-standard trajectory formats, specify the column mapping:

```
TRAJECTORY_HEADER=ITEM: ATOMS id type x y z
```

### Benchmarking

A separate benchmark program is included to test performance:

```bash
./benchmark_coordination
```

This will test different thread counts and report performance statistics.

## Troubleshooting

### Common Issues

1. **Missing atoms in output**
   - Verify atom selection is correct
   - Check if atom types match those in data file

2. **Incorrect coordination numbers**
   - Verify cutoff distances are appropriate
   - Check for periodic boundary issues
   - Verify units in data/trajectory files

3. **Performance issues**
   - Try adjusting `CELL_UPDATE_FREQ`
   - Enable `SELECTIVE_REBUILD=yes`
   - Check if `USE_TRAJECTORY_BOX` is needed

### Cell List Statistics

The program reports detailed cell list statistics at the end of the run:

```
CELL LIST STATISTICS
-------------------
Complete cell list rebuilds:      42
Cell list maintenance updates:    158
Selective cell rebuilds:          107
Average empty cells:              4071 (48.6%)
Average atoms per non-empty cell: 4.45
Min/Max atoms per cell:           1/12
```

These statistics help diagnose performance issues:
- High rebuild counts may indicate unstable simulation or inappropriate update frequency
- High empty cell percentage may suggest inefficient partitioning
- Very high max atoms per cell may indicate clustering issues

## Technical Details

### Project Structure

```
├── types_mod.f90           # Core data types and structures
├── error_mod.f90           # Error handling utilities
├── config_mod.f90          # Configuration management
├── atom_selection_mod.f90  # Atom selection utilities
├── cell_list_mod.f90       # Spatial partitioning
├── coordination_mod.f90    # Core coordination calculation
├── data_file_io_mod.f90    # LAMMPS data file parsing
├── trajectory_io_mod.f90   # Trajectory file reading
├── output_io_mod.f90       # Results output
├── output_cache_mod.f90    # Output caching for performance 
├── progress_mod.f90        # Progress display
├── io_mod.f90              # Integrated I/O management
├── main.f90                # Main program
└── benchmark_main.f90      # Performance benchmarking
```


## Contact

For questions, bug reports, or feature requests:
- Email: author@example.com
- GitHub Issues: https://github.com/username/coordination-analysis/issues
