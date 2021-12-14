#!/bin/bash

session=$1
window=$2
pan1=$3
pan2=$4

#Get width and lenght size of terminal, this is needed if one wants to resize a detached session/window/pane
#with resize-pane command here
set -- $(stty size) #$1=rows, $2=columns

#start a new session in dettached mode with resizable panes
tmux new-session -s $session -n $window -d -x "$2" -y "$(($1 - 1))"
tmux send-keys -t $session 'echo "first command in 1st pane"' C-m

#rename pane 0 with value of $pan1
#tmux select-pane -t $session:$window.0 -T "$pan1"
tmux set -p @mytitle "$pan1"

#split window vertically
tmux split-window -h
tmux send-keys -t $session 'echo "first command in 2nd pane"' C-m
#tmux select-pane -t $session:$window.1 -T "$pan2"
tmux set -p @mytitle "$pan2"

#At the end, attach to the customized session
tmux attach -t $session
