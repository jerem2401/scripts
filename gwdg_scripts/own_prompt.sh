[ "$HOST" = "" ]  && HOST=`hostname | cut -d. -f1`

case "$HOST" in
   gwdu103)
     psname="Broadwell";;
   gwdu101)
     psname="Interlagos";;
   gwdu102)
     psname="Sandy-bridge";;
   gwdu05)
     psname=login ;;
   gwdu60)
     psname=gwdu60 ;;
   jj28l??)
      psname=Juropa ;;
   *)
     psname=""
esac 


export PS1="\[\033[01;38;5;91m\]DarkMagician@\h[\[\033[0m\]\[\033[01;38;5;9m\]${psname}\[\033[0m\]\[\033[01;38;5;91m\]]\[\033[0m\]\[\033[01;38;5;15m\]:\w\[\033[0m\]\n\[\033[01;38;5;91m\]$ \[\033[0m\]"
#"\[\033[01;38;5;91m\]DarkMagician\[\033[0m\]\[\033[01;38;5;91m\]@\[\033[0m\]\[\e[01;38;5;112m\]\h\[\033[0m\]\[\033[01;38;5;91m\]:\[\033[0m\]\[\033[01;38;5;15m\]\w\[\033[0m\]\[\033[01;38;5;91m\]\\$\[\033[0m\] "
