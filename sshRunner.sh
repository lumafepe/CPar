#!/bin/bash

# Set username and password as variables
username=""
password=""

# Set the destination folder on the remote server
destination_folder="idfk"

# Check if the first argument is provided
if [ -z "$1" ]; then
    echo "Usage: $0 [par|seq]"
    exit 1
fi


# Set the number of CPUs based on the first argument
if [ "$1" == "par" ]; then
    cpus=40
    runner="parrunner.sh"
elif [ "$1" == "seq" ]; then
    cpus=1
    runner="seqrunner.sh"
else
    echo "Invalid argument. Use 'par' or 'seq'."
    exit 1
fi

# Clear the contents of the destination folder on the remote server
sshpass -p "$password" ssh "$username"@s7edu.di.uminho.pt "rm -rf ~/$destination_folder/*"

# Copy the contents of the current folder to the remote server
rsync -avz -e "sshpass -p '$password' ssh" ./ "$username"@s7edu.di.uminho.pt:~/$destination_folder/ >> /dev/null

# SSH into the remote server and execute the desired commands
sshpass -p "$password" ssh "$username"@s7edu.di.uminho.pt << EOF
cd ~/$destination_folder/

sbatch --partition=cpar --cpus-per-task=$cpus -W $runner
for file in slurm-*; do
    if [[ -f "\$file" ]]; then
        echo "Processing file: \$file"
        cat "\$file"
    fi
done
EOF
