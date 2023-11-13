#!/bin/bash

# Set username and password as variables
username=""
password=""

# Set the destination folder on the remote server
destination_folder="idfk"

# Clear the contents of the destination folder on the remote server
sshpass -p "$password" ssh "$username"@s7edu.di.uminho.pt "rm -rf ~/$destination_folder/*"

# Copy the contents of the current folder to the remote server
rsync -avz -e "sshpass -p '$password' ssh" ./ "$username"@s7edu.di.uminho.pt:~/$destination_folder/ >> /dev/null

# SSH into the remote server and execute the desired commands
sshpass -p "$password" ssh "$username"@s7edu.di.uminho.pt << EOF
cd ~/$destination_folder/
sbatch --partition=cpar --cpus-per-task=40 -W runner.sh
for file in slurm-*; do
    if [[ -f "\$file" ]]; then
        echo "Processing file: \$file"
        cat "\$file"
    fi
done
EOF