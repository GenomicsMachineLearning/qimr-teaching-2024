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

cd /scratch/$USER/qimr-teaching-2024
jupyter notebook --no-browser --port=${port} --ip=0.0.0.0
