#!/bin/bash
running=true
touch dingdong.txt

while [ $running ]
do 
    (qme | grep pAu | wc -l)>dingdong.txt
    while read line
    do
	count=$line
    done < dingdong.txt
    if [ $count -eq 0 ];
    then
	running=false
	printf '\7'
	echo "  "
	echo "    _____    ____   _   _  ______  _ "
	echo "   |  __ \  / __ \ | \ | ||  ____|| |"
	echo "   | |  | || |  | ||  \| || |__   | |"
	echo "   | |  | || |  | || . ' ||  __|  | |"
	echo "   | |__| || |__| || |\  || |____ |_|"
	echo "   |_____/  \____/ |_| \_||______|(_)"
	break
    else
	sleep 1
    fi
done

rm dingdong.txt
