# Script of automated MassQLab processing based on configuration of massqlab_config.json

# pyinstaller command
# pyinstaller --noconfirm --noupx -F --console --collect-all "massql" --collect-all "matchms" --collect-all "pyarrow" --collect-all "pymzml" --exclude-module "kaleido"  "<absolute_path_to_script>"

# Convert jupyter notebok to script
# jupyter nbconvert --to script "<absolute_path_to_notebook>.ipynb"

from MassQLab_workflow import *
run()