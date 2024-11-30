# HocusPockets: Holistic Computacional Utility for Scoring Pockets

HocusPockets is a tool designed to prioritize biological targets within the proteome of microorganisms based on the druggability of interaction sites and their chemical reactivity, utilizing descriptors such as Fukui indices to evaluate potential binding sites.

## Prerequisites

Before using HocusPockets, ensure the following packages and external programs are installed:

### Python packages

> * os  
> * sys  
> * argparse  
> * subprocess  
> * textwrap  
> * shutil  
> * csv  
> * pandas  
> * numpy  
> * scipy.spatial  

### External software

> * [fpocket](https://github.com/Discngine/fpocket/tree/master)  
> * [ConCavity](https://compbio.cs.princeton.edu/concavity/)  
> * [RADI](https://hal.science/hal-03798049)  

## Installing and Configuring RADI

To properly install and use the **RADI** software, follow these steps:

### 1. Download and Unpack the Source Files

```bash
tar -xzvf file.tar.gz
```

### 2. Install Intel Fortran Compiler

Install the Intel Fortran Compiler (ifx) via the Intel OneAPI repository:
```bash
sudo apt install intel-oneapi-compiler-fortran
```

Configure the environment for Intel compilers:
```bash
source /opt/intel/oneapi/setvars.sh
```

### 3. Prepare and Compile RADI Libraries

#### libmol

- **Step 1:** Create a directory for `libmol` and move the source file:
    ```bash
    mkdir libmol
    mv libmol.t libmol
    cd libmol
    ```
- **Step 2:** Extract the archive:
    ```bash
    ar -x libmol.t
    ```
- **Step 3:** Compile the source code:
    ```bash
    ifx -c *.f
    ```
- **Step 4:** Create the `libmol.a` library:
    ```bash
    ar rcs libmol.a *.o
    ```

#### libuti

- **Step 1:** Create a directory for `libuti` and move the source file:
    ```bash
    mkdir libuti
    mv libuti.t libuti
    cd libuti
    ```
- **Step 2:** Extract the archive:
    ```bash
    ar -x libuti.t
    ```
- **Step 3:** Compile the source code:
    ```bash
    ifx -c *.f
    ```
- **Step 4:** Create the `libuti.a` library:
    ```bash
    ar rcs libuti.a *.o
    ```

#### libfor

- **Step 1:** Create a directory for `libfor` and move the source file:
    ```bash
    mkdir libfor
    mv libfor.t libfor
    cd libfor
    ```
- **Step 2:** Extract the archive:
    ```bash
    ar -x libfor.t
    ```
- **Step 3:** Copy the previously compiled libraries (`libmol.a` and `libuti.a`) to the `libfor` directory:
    ```bash
    cp ../libmol/libmol.a .
    cp ../libuti/libuti.a .
    ```
- **Step 4:** Compile the RADI executable:
    ```bash
    ifx -c RADI.f
    ifx RADI.o libmol.a libuti.a -o radi
    ```
    
## Usage

You can use HocusPockets with either a .pdb or .cif file format:
```bash
python3 ./HocusPockets/ -i file.pdb -o output_path
```

For a full list of options, run:
```bash
python3 ./HocusPockets/ -h
```
