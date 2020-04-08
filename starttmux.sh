#!/bin/sh

session=$1

#start session, -s sessions name
tmux new-session -s $session -d
tmux send-keys -t $session 'cd ~/gitrepo/scripts' C-m
tmux send-keys -t $session 'clear' C-m
tmux split-window -h
tmux send-keys -t $session 'clear' C-m
tmux split-window -v
tmux send-keys -t $session 'clear' C-m
tmux attach -t $session
#tmux select-pane -t 1
#tmux split-window -v
#tmux attach -t python
#tmux new-window -c ~/gitrepo/scripts -n vim
#tmux select-window -t $session:0
#tmux split-window -h
#tmux attach -t $session
#tmux send-keys -t 'cd ~/gitrepo/scripts' C-m
#tmux selectp -t 1
#tmux send-keys "" C-m

