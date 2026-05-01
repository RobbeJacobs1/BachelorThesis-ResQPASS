import shutil
# This script is used to copy the "Results" directory to "MatLab" and then remove the original "Results" directory. This is necessary because the "Results" directory contains files that 
# are needed for the MatLab code to run.

src = "Results"
dst = "MatLab"

shutil.copytree(src, dst, dirs_exist_ok=True)
