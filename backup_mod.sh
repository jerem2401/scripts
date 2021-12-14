#!/bin/bash
#This script is inspired from Robert Becker and Michael Ramsey scripts (https://whattheserver.com/keep-your-sshfs-mounts-mounted)


#function for calculating the difference between dates in unit days
day_difference (){
    difference=$((($(date -d $1 +%s) - $(date -d $2 +%s))/86400))
}

#save current date as YYYY-MM-DD
new_backup=$(date '+%Y-%m-%d')

#echoing start date for log
echo -e "\nBeginning of backup on $new_backup" >> /media/jeremy/jerem_backup_2/backup_history.log


########################################################Sanity check#####################################################################
#Check if backup hardrive and smaug is mounted, this will also try to mount smaug if not mounted
#Careful: you need to setup passwordless ssh and add a check file named "is_mounted" in both directories

#USER SPECIFIC VARIABLES
#path to backup hardrive
myemail='jeremy.lapierre@uni-saarland.de'
harddrive='/media/jeremy/jerem_backup_2/'
#Set file name to test for in the backup hard drive. Just create an empty file "is_mounted" in hard drive so we can check if it exists.
harddrivetestfile='/media/jeremy/jerem_backup_2/is_mounted'

#Set file name to test for in the smaug mouted directory. Just create an empty file "is_mounted" IN THE DIRECTORY YOU MOUNT  so we can check if it exists.
remotemounttestfile="/home/jeremy/mnt/smaug/is_mounted"
#path to local dir where you mount smaug
smaugmntlocal='/home/jeremy/mnt/smaug'
#path of smaug directory you want to mount
smaugmntremote='/data/users/jeremy/simulation'
#Smaug user id
useridsmaug='jeremy'

if mountpoint -q $harddrive && [ -f $harddrivetestfile ]; then
    echo "Your hard drive is mounted correctly" >> "${harddrive}/backup_history.log"
else
    echo "Your hard drive is not mounted correctly, exiting backup script" >> "${harddrive}/backup_history.log"
    mail -s "Daily backup" $myemail <<< 'Backup not done, hard drive is not mounted correctly'
    exit 1
fi

if mountpoint -q $smaugmntlocal && [ -f $remotemounttestfile ]; then
    echo "Your smaug directory is mounted correctly" >> "${harddrive}/backup_history.log"
else
   echo "Smaug not mounted properly" >> "${harddrive}/backup_history.log"
   #umount gracefully if possible
   umount $smaugmntlocal  > /dev/null 2>&1

   #kill any frozen process on the mount
   fuser -k $smaugmntlocal  > /dev/null 2>&1
   fusermount -u $smaugmntlocal  > /dev/null 2>&1
   umount -l $smaugmntlocal  > /dev/null 2>&1
   umount $smaugmntlocal  > /dev/null 2>&1

   #remount smaug
   sshfs smaug:$smaugmntremote $smaugmntlocal -o ssh_command="ssh ${useridsmaug}@alef.lusi.uni-sb.de -t ssh"

   if mountpoint -q $smaugmntlocal && [ -f $remotemounttestfile ]; then
       echo "Your smaug directory has been mounted" >> "${harddrive}/backup_history.log"
   else
       echo "Smaug not mounted properly needs fixed manually"
       mail -s "Daily backup" $myemail <<< 'Backup not done, smaug not mouted'
       exit 1
   fi
fi
########################################################Sanity check#####################################################################
#scp bash history to back this up
echo "scp-ing smaug bash history" >> "${harddrive}/backup_history.log"
scp jeremy@alef.lusi.uni-sb.de:/home/users/jeremy/.bash_history /home/jeremy/mnt/smaug/

#go to your external hard drive
cd /media/jeremy/jerem_backup_2/

#get the latest backup folder by look for pattern XXXX-XX-XX
latest_backup=$(ls | grep -e "[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]" | tail -n 1)


#write a log entry for your new backup
echo -e "Used ${latest_backup} as template for ${new_backup}" >> /media/jeremy/jerem_backup_2/backup_history.log

#create a new backup with hard links of your latest backup
rsync -avzHP --exclude ".*" --exclude anaconda3 --exclude Desktop --exclude Downloads --exclude opt --exclude Pictures --exclude mnt/gwdg --exclude mnt/cip --exclude src --exclude Videos --exclude Music --exclude Public --exclude snap --exclude Templates --delete --link-dest=$PWD/${latest_backup} /home/jeremy/ /media/jeremy/jerem_backup_2/${new_backup} >> /media/jeremy/jerem_backup_2/backup_history.log 2>&1

#save exit code of rsync
EC=$?

#send email to notify non 0 exit code
if (($EC > 0)); then
    mail -s "Daily backup" jeremy.lapierre@uni-saarland.de <<< 'rsync had non 0 exit code, you should check the log file'
fi



#Now comes the cleanup part after the rsync finished
wait

#echoing for log file
echo "Beginning of clean part of backup script" >> /media/jeremy/jerem_backup_2/backup_history.log

#change variable name for comprehension
today=$new_backup

unset -v 'days'
declare -a days

for i in "[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]"; do
    days+=( "$i" )
done

#setup arrays
declare -a one_week_days
declare -a two_week_days
declare -a month_days
declare -a half_year_days
declare -a year_days

#store the directories (dates) in different arrays which resemble a certain "age-range" (week, 2 weeks, month, half year and older than half year)

for day in ${days[@]}
do
    day_difference $today $day
    if [[ $difference -ge 7 && $difference -lt 14 ]]
    then
	one_week_days+=( "$day" )
    elif [[ $difference -ge 14 && $difference -lt 21 ]]
    then
	two_week_days+=( "$day" )
    elif [[ $difference -ge 21 && $difference -lt 30 ]]
    then
	month_days+=( "$day" )
    elif [[ $difference -ge 30 && $difference -lt 180 ]]
    then
	half_year_days+=( "$day" )
    elif [[ $difference -ge 180 ]]
    then
	year_days+=( "$day" )
    fi
done

#store the length of the arrays
len1week=${#one_week_days[@]}
len2week=${#two_week_days[@]}
len1month=${#month_days[@]}
len6month=${#half_year_days[@]}
len1year=${#year_days[@]}


#remove every folder except the oldest backup in each specific "age-range"
#This means folders will "travel' (or age) through all those ranges, so you always have one folder in each of those ranges
#Thus: rm folder from index 0 to lenght of age-range array minus one
#first line is sanity check

echo "one_week_days are: ${one_week_days[@]}, two_week_days are: ${two_week_days[@]}, month_days are ${month_days[@]}, half_year_days are: ${half_year_days[@]}, year_days are: ${year_days[@]}" >> /media/jeremy/jerem_backup_2/backup_history.log

for i in "${one_week_days[@]:1:$len1week-1}"
do
    echo "removing: $i" >> /media/jeremy/jerem_backup_2/backup_history.log
    rm -r $i
done

for j in "${two_week_days[@]:1:$len2week-1}"
do
    echo "removing: $j" >> /media/jeremy/jerem_backup_2/backup_history.log
    rm -r $j
done

for k in "${month_days[@]:1:$len1month-1}"
do
    echo "removing: $k" >> /media/jeremy/jerem_backup_2/backup_history.log
    rm -r $k
done

for l in "${half_year_days[@]:1:$len6month-1}"
do
    echo "removing: $l" >> /media/jeremy/jerem_backup_2/backup_history.log
    rm -r $l
done

for m in "${year_days[@]:1:$len1year-1}"
do
    echo "removing: $m" >> /media/jeremy/jerem_backup_2/backup_history.log
    rm -r $m
done

echo "End of cleaning part and end of backup of date: $today" >> /media/jeremy/jerem_backup_2/backup_history.log
