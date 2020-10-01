#!/bin/bash

session=$1

#Get width and lenght size of terminal, this is needed in one wants to resize a dettach session/window/pane
#with resize-pane command here
set -- $(stty size) #$1=rows, $2=columns
host=$(hostname)

tmux new-session -s $session -d -x "$2" -y "$(($1 - 1))" # status line uses a row

if [[ $session == pythonh* ]]; then
    tmux send-keys -t $session 'cd ~/gitrepo/scripts' C-m
    tmux send-keys -t $session 'clear' C-m
    tmux split-window -h
    tmux send-keys -t $session 'clear' C-m
    if [[ $host == gwdu* ]]; then
        tmux send-keys -t $session 'module load conda' C-m
        tmux send-keys -t $session 'source activate env1' C-m
        tmux send-keys -t $session 'python' C-m
    else
        tmux send-keys -t $session 'python' C-m
    fi
    tmux split-window -v
    tmux send-keys -t $session 'clear' C-m
    tmux resize-pane -t $session:0.0 -x 105
    #tmux attach -t python
    #tmux new-window -c ~/gitrepo/scripts -n vim
    #tmux select-window -t $session:0
    #tmux split-window -h
    #tmux attach -t $session
    #tmux send-keys -t 'cd ~/gitrepo/scripts' C-m
    #tmux selectp -t 1
    #tmux send-keys "" C-m
elif [[ $session == pythonv* ]]; then
    tmux send-keys -t $session 'cd ~/gitrepo/scripts' C-m
    tmux send-keys -t $session 'clear' C-m
    tmux split-window -v
    tmux send-keys -t $session 'clear' C-m
    if [[ $host == gwdu* ]]; then
        tmux send-keys -t $session 'module load conda' C-m
        tmux send-keys -t $session 'source activate env1' C-m
        tmux send-keys -t $session 'python' C-m
    else
        tmux send-keys -t $session 'python' C-m
    fi
    tmux split-window -h
    tmux send-keys -t $session 'clear' C-m
    tmux resize-pane -t $session:0.0 -y 86
elif [[ $session == gwdg ]]; then
    tmux send-keys -t $session 'cogwdg' C-m
    tmux send-keys -t $session 'clear' C-m
    tmux split-window -v
    tmux send-keys -t $session 'gwdg' C-m
    tmux send-keys -t $session 'clear' C-m
    tmux split-window -h
    tmux send-keys -t $session 'clear' C-m
elif [[ $session == smaug ]]; then
    tmux send-keys -t $session 'cosmaug' C-m
    tmux send-keys -t $session 'clear' C-m
    tmux split-window -v
    tmux send-keys -t $session 'smaug' C-m
    tmux send-keys -t $session 'clear' C-m
    tmux split-window -h
    tmux send-keys -t $session 'clear' C-m
else
    tmux send-keys -t $session 'clear' C-m
    tmux split-window -v
    tmux send-keys -t $session 'clear' C-m
    tmux split-window -h
    tmux send-keys -t $session 'clear' C-m
fi

tmux attach -t $session
