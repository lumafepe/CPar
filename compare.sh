#!/bin/bash

# ANSI color codes
RED="\e[31m"       # Red text
GREEN="\e[32m"     # Green text
YELLOW="\e[33m"    # Yellow text
RESET="\e[0m"      # Reset text color

printf "\n${YELLOW}Running diff checker...\n${RESET}"

diff_cp_average=$(diff -u cp_average.txt ./output/cp_average.txt 2>&1)
diff_cp_output=$(diff -u cp_output.txt ./output/cp_output.txt 2>&1)
diff_cp_traj=$(diff -u cp_traj.xyz ./output/cp_traj.xyz 2>&1)

# Step 4: Show differences, if any
if [ -z "$diff_cp_average" ] && [ -z "$diff_cp_output" ] && [ -z "$diff_cp_traj" ]; then
    echo -e "\n${GREEN}Good Job! No differences found.${RESET}"
else
    [ -n "$diff_cp_average" ] && echo -e "Difference in cp_average.txt:\n$diff_cp_average" | less -R
    [ -n "$diff_cp_output" ] && echo -e "Difference in cp_output.txt:\n$diff_cp_output" | less -R
    [ -n "$diff_cp_traj" ] && echo -e "Difference in cp_traj.txt:\n$diff_cp_traj" | less -R

    [ -n "$diff_cp_average" ] && printf "${RED}cp_average differs\n${RESET}"
    [ -n "$diff_cp_output" ] && printf "${RED}cp_output differs\n${RESET}"
    [ -n "$diff_cp_traj" ] && printf "${RED}cp_traj differs\n${RESET}"
fi


# Examples:
# sh run_and_compare.sh "make run"
# sh run_and_compare.sh "make run-og"
# sh run_and_compare.sh "make run-vect"
# sh run_and_compare.sh <qualquer comando que gere os ficheiros>


