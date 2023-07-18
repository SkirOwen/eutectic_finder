# Eutectic Finder


## Installation
### pip

```shell
pip install -r requirements.txt
```

### Conda

`conda`/`mamba`

```shell
conda install -c conda-forge --file requirements.txt
```


## Usage
Use the `-c`/`--celsius` flag to set to use the temperature in Celsius.
**The default use Kelvin.**

### Using CSV
Using the `-f`/`--file` flag to specify the `.csv` to use for a run.  
The flag `-r`/`--run` for the run number.  
More information in the section `DATA`.  

```shell
python main.py -f data.csv -r 1
```


### Manual/CLI

`-e`/`--enthalpy` to input the enthalpy of each element, separated by a space (in J/mol).  
`-t`/`--temperature` to input the temperature of each element, separated by a space, (**default in Kelvin**).
See the `-c` flag to toggle to Celsius.

```shell
python main.py -m -e 11300 46000 -t 1234 1688
```

### Prompt
To give the input as a prompt use:
`-c` can be used.

```shell
python main.py 
```

## DATA

The `DATA` can be read form a `.csv`.
The file need to be of the following format:

```csv
run,element,enthalpy,temp_c,temp_k,exp_x,cal_x
1,Ag,11300,960.85,1234,85,85
1,Si,46000,1414.85,1688,15,15
2,KCl,25500,772,-1,24,22
2,LiCl,13400,610,-1,43,61
2,NaCl,28500,808,-1,33,17
```

`run` specify the numbers of element in a mixture.  
`element` is the symbol of the element present.  
`enthalpy` is the melting enthalpy in J/mol.  
`temp_c` is the melting temperature in Celsius, (if not know, set -300).   
`temp_k` is the melting temperature in Kelvin, (if not know, set -300).  
`exp_x` is the experimental molar fraction at en eutectic, (if not know, set -300).  
`cal_c` is the calculate molar fraction at en eutectic, (if not know, set -300).


## License

This code is under Apache 2.0