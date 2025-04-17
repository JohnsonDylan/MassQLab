## Setup Guide (Virtual Environment)

Follow these steps to set up a virtual environment and install dependencies for MassQLab.

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

### 3. Install virtualenv

```bash
pip install virtualenv
```

### 4. Create and Activate the Virtual Environment

Navigate to the project directory if needed:

```bash
cd path\to\MassQLab
```

Create and activate the environment:

```bash
py -3.9 -m venv env
env\Scripts\activate
```

You should now see `(env)` at the start of your terminal prompt.

### 5. Install Dependencies

```bash
pip install -r requirements.txt
```

You're now ready to use the notebooks in `notebooks/` or run the workflow using provided scripts.
