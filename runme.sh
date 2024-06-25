#!/bin/bash

node=$(hostname -I | cut -f1 -d' ')
user=$(whoami)

for port in {3500..3700}; do
  if ! ss -lptn | grep -q ":$port"; then
    echo "Port $port is available"
    break
  fi
done

echo -e "
Command to create ssh tunnel:
ssh -N -L ${port}:${node}:${port} ${user}@${node}
Use a Browser on your local machine to go to:
localhost:${port}  (prefix w/ https:// if using password)
"

target_dir="/scratch/$USER/qimr-teaching-2024/data"
softlink_name="/data/module2/data"

if [ ! -d "$target_dir" ]; then
  ln -s "$softlink_name" "$target_dir"
fi

# Fix for GLIBC error
# https://stackoverflow.com/questions/58424974/anaconda-importerror-usr-lib64-libstdc-so-6-version-glibcxx-3-4-21-not-fo
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CONDA_PREFIX/lib

cd /scratch/$USER/qimr-teaching-2024
jupyter notebook --no-browser --port=${port} --ip=0.0.0.0
