# Setup Guide

Follow these steps to set up a virtual environment and install dependencies for MassQLab. (Or see [Quickstart](quickstart.md))

### 1. Download the Project

**Option A: Clone with Git**  
(Requires [Git](https://git-scm.com/) to be installed)

```bash
git clone https://github.com/JohnsonDylan/MassQLab.git
cd MassQLab
```

**Option B: Download ZIP from GitHub**  
1. Go to [github.com/JohnsonDylan/MassQLab](https://github.com/JohnsonDylan/MassQLab)  
2. Click the green **"Code"** button → **"Download ZIP"**  
3. Extract the ZIP and open a command prompt in the extracted `MassQLab` folder

### 2. Check Python Installation

Run:

```bash
py --list
```

Make sure Python 3.9 is listed. If not:

- Install Python 3.9 from: https://www.python.org/downloads/
- During installation, check the box **"Add Python to PATH"**

If you get an error, Python may not be installed or not added to your system path.

### 3. Use virtualenv OR conda environemnt 
#### 3.A. Create and Activate a Virtual Environment

Ensure virtualenv is installed
```bash
pip install virtualenv
```

Navigate to the project directory:

```bash
cd path\to\MassQLab
```

Create and activate the virtual environment:

```bash
py -3.9 -m venv env
env\Scripts\activate
```

You should now see `(env)` at the start of your terminal prompt.

#### 3.B. Create and Activate a Conda Environment

Ensure Anaconda is installed
https://www.anaconda.com/docs/getting-started/anaconda/install

Navigate to the project directory:

```bash
cd path\to\MassQLab
```

Create and activate the conda environment:

```bash
conda create -n env python=3.9
conda activate env
```

You should now see `(env)` at the start of your terminal prompt.

### 4. Install Dependencies

The following may be sufficient:
```bash
pip install jupyterlab massql reportlab openpyxl matplotlib
```

Otherwise use:
```bash
pip install -r requirements.txt
```

You're now ready to use the notebooks in `notebooks/` or run the workflow using provided scripts.
____
