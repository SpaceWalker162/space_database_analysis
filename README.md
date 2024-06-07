# Space Database and Analysis
## Acknowledgment

Please use the following to acknowledge use of Space Database and Analysis in your publications:

> Data analysis was performed using the Space Database and Analysis package available at https://github.com/SpaceWalker162/space_database_analysis

## Installation

```sh
git clone https://github.com/SpaceWalker162/space_database_analysis.git
```

Any time you want to update the code to the latest version, run from the command line

```sh
git pull
```

Install dependent packages not included in the anaconda base environment
```sh
conda activate
pip install cdasws cdflib
```

If you do not use command line there are github programs for Windows and Mac, see github web page.

## Usage

Each time starting new Python3 session execute in Matlab:

```python
import dataAnalysisTools as dat
```

To use a routine/function, such as `mvab`, execute in Python3 Session:

```python
eigen_system = dat.mvab(magnetic_field_data)
```

## Issues

If you experience any issues with Space Database and Analysis please [submit an issue](https://github.com/SpaceWalker162/space_database_analysis/issues) via our github page.

When submitting a new issue, please provide information required to replicate the issue as well as information regarding operating system and version of Python used.
