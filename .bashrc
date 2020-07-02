# ~/.bashrc: executed by bash(1) for non-login shells.
# see /usr/share/doc/bash/examples/startup-files (in the package bash-doc)
# for examples

# If not running interactively, don't do anything
case $- in
    *i*) ;;
      *) return;;
esac

#History
export HISTSIZE=1000000000000000000000
export HISTFILESIZE=100000000000000000
export HISTCONTROL=erasedups:ignoredups
shopt -s histappend
PROMPT_COMMAND="history -n; history -w; history -c; history -r; $PROMPT_COMMAND"

# check the window size after each command and, if necessary,
# update the values of LINES and COLUMNS.
shopt -s checkwinsize

# If set, the pattern "**" used in a pathname expansion context will
# match all files and zero or more directories and subdirectories.
#shopt -s globstar

# make less more friendly for non-text input files, see lesspipe(1)
[ -x /usr/bin/lesspipe ] && eval "$(SHELL=/bin/sh lesspipe)"

# set variable identifying the chroot you work in (used in the prompt below)
if [ -z "${debian_chroot:-}" ] && [ -r /etc/debian_chroot ]; then
    debian_chroot=$(cat /etc/debian_chroot)
fi

# set a fancy prompt (non-color, unless we know we "want" color)
case "$TERM" in
    xterm-color|*-256color) color_prompt=yes;;
esac

# uncomment for a colored prompt, if the terminal has the capability; turned
# off by default to not distract the user: the focus in a terminal window
# should be on the output of commands, not on the prompt
#force_color_prompt=yes

if [ -n "$force_color_prompt" ]; then
    if [ -x /usr/bin/tput ] && tput setaf 1 >&/dev/null; then
	# We have color support; assume it's compliant with Ecma-48
	# (ISO/IEC-6429). (Lack of such support is extremely rare, and such
	# a case would tend to support setf rather than setaf.)
	color_prompt=yes
    else
	color_prompt=
    fi
fi

if [ "$color_prompt" = yes ]; then
    PS1='${debian_chroot:+($debian_chroot)}\[\033[01;32m\]\u@\h\[\033[00m\]:\[\033[01;34m\]\w\[\033[00m\]\$ '
else
    PS1='${debian_chroot:+($debian_chroot)}\u@\h:\w\$ '
fi
unset color_prompt force_color_prompt

# If this is an xterm set the title to user@host:dir
case "$TERM" in
xterm*|rxvt*)
    PS1="\[\e]0;${debian_chroot:+($debian_chroot)}\u@\h: \w\a\]$PS1"
    ;;
*)
    ;;
esac

# enable color support of ls and also add handy aliases
if [ -x /usr/bin/dircolors ]; then
    test -r ~/.dircolors && eval "$(dircolors -b ~/.dircolors)" || eval "$(dircolors -b)"
    alias ls='ls --color=auto'
    #alias dir='dir --color=auto'
    #alias vdir='vdir --color=auto'

    alias grep='grep --color=auto'
    alias fgrep='fgrep --color=auto'
    alias egrep='egrep --color=auto'
    alias tmuxx='bash ~/gitrepo/scripts/starttmux.sh'
fi

# colored GCC warnings and errors
#export GCC_COLORS='error=01;31:warning=01;35:note=01;36:caret=01;32:locus=01:quote=01'

# some more ls aliases
alias ll='ls -alF'
alias la='ls -A'
alias l='ls -CF'

# Add an "alert" alias for long running commands.  Use like so:
#   sleep 10; alert
alias alert='notify-send --urgency=low -i "$([ $? = 0 ] && echo terminal || echo error)" "$(history|tail -n1|sed -e '\''s/^\s*[0-9]\+\s*//;s/[;&|]\s*alert$//'\'')"'

# Alias definitions.
# You may want to put all your additions into a separate file like
# ~/.bash_aliases, instead of adding them here directly.
# See /usr/share/doc/bash-doc/examples in the bash-doc package.

if [ -f ~/.bash_aliases ]; then
    . ~/.bash_aliases
fi

# enable programmable completion features (you don't need to enable
# this, if it's already enabled in /etc/bash.bashrc and /etc/profile
# sources /etc/bash.bashrc).
if ! shopt -oq posix; then
  if [ -f /usr/share/bash-completion/bash_completion ]; then
    . /usr/share/bash-completion/bash_completion
  elif [ -f /etc/bash_completion ]; then
    . /etc/bash_completion
  fi
fi

# added by Anaconda3 installer
export PATH="/home/lapierre/anaconda3/bin:$PATH"

#locale
export LC_ALL="en_US.UTF-8"
export LANG="en_US.UTF-8"
export LANGUAGE="en_US.UTF-8"

#Python path
export PYTHONPATH="/home/lapierre/gitrepo/scripts:$PYTHONPATH"
export PYTHONPATH="/home/lapierre/gitrepo/scripts/wham:$PYTHONPATH"
export PYTHONPATH="/home/lapierre/gitrepo/scripts/miscellaneous:$PYTHONPATH"

export PATH="/home/lapierre/gitrepo/scripts:$PATH"
export PATH="/home/lapierre/gitrepo/scripts/wham:$PATH"

#alias
alias unmount-gwdg='fusermount -uz ~/mnt/gwdg'
alias mount-gwdg='sshfs jlapier@transfer.gwdg.de:/usr/users/jlapier/simulation ~/mnt/gwdg -o auto_cache,reconnect,ServerAliveInterval=500,ServerAliveCountMax=2'
alias unmount-smaug='fusermount -uz ~/mnt/smaug'
alias mount-smaug="sshfs smaug:/data/users/jeremy/simulation /home/lapierre/mnt/smaug/ -o ssh_command='ssh jeremy@alef.lusi.uni-sb.de -t ssh'"
alias ls='ls -ltr --color=auto'
alias sos='cat ~/Documents/command.txt'
alias cosmaug='ssh -t jeremy@alef.lusi.uni-sb.de "ssh -t smaug \"cd /data/users/jeremy || exit 192; exec /bin/bash -login \""'
alias cogwdg='ssh -t -A jlapier@login.gwdg.de ssh -t gwdu103'
alias codema='ssh -t -A jeremy@alef.lusi.uni-sb.de ssh -t dema9'
alias smaug='cd /home/lapierre/mnt/smaug'
alias gwdg='cd ~/mnt/gwdg'
alias cippool='ssh -t -A -XC jela004@gordon.cip.physik.uni-saarland.de "ssh -t -XC 134.96.87.205"'
alias mount-cip='sshfs jela004@134.96.87.194:/home/jela004/tuto ~/mnt/cip -o auto_cache,reconnect'
alias unmount-cip='fusermount -uz ~/mnt/cip'

#plumed
export PATH="~/bin/bin/:$PATH"
export LD_LIBRARY_PATH='~/bin/lib/':$LD_LIBRARY_PATH
export INCLUDE='~/bin/include':$INCLUDE
export plumedir='~/bin/bin/plumed'
