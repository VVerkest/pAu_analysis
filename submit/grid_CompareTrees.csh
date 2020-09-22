#!/bin/csh

# used to submit sequential jobs on the grid

# first make sure program is updated and exists
make bin/CompareTrees || exit

set noglob
    
set ExecPath = `pwd`
set execute = './bin/CompareTrees'
set maxevents = 920795

set OutBase = CompareTrees
set uscore = "_"
    
set fileNo = 0
set nEventsPerFile = 50000
set startEvent = ${fileNo}
set firstArg = 'start'

set jobName = ${OutBase}${uscore}${fileNo}
    
# Make the output file                                                                                                                                               
set outName = ${OutBase}/${OutBase}${uscore}${startEvent}.root

# Logfiles. Thanks cshell for this "elegant" syntax to split err and out                                                                                                       
set LogFile     = ${OutBase}/${OutBase}${uscore}${startEvent}.log
set ErrFile     = ${OutBase}/${OutBase}${uscore}${startEvent}.err
    
echo "Logging output to " $LogFile
echo "Logging errors to " $ErrFile
echo "now submitting this script: "
    
echo sbatch --mem-per-cpu=32GB -q express -p erhip --nodes=4 -o $LogFile -e $ErrFile --time=01:00:00 --job-name=${jobName} -- ${ExecPath}/submit/qwrap.sh ${ExecPath} $execute $startEvent

sbatch --mem-per-cpu=32GB -q express -p erhip --nodes=4 -o $LogFile -e $ErrFile --time=01:00:00 --job-name=${jobName} -- ${ExecPath}/submit/qwrap.sh ${ExecPath} $execute $startEvent

@ fileNo = $fileNo + 1

echo " "
    
set startEvent = $nEventsPerFile
    
while ( $startEvent <= 920795 )

set jobName = ${OutBase}${uscore}${fileNo}
    
# Make the output file                                                                                                                                               
set outName = ${OutBase}/${OutBase}${uscore}${startEvent}.root

# Logfiles. Thanks cshell for this "elegant" syntax to split err and out                                                                                                       
set LogFile     = ${OutBase}/${OutBase}${uscore}${startEvent}.log
set ErrFile     = ${OutBase}/${OutBase}${uscore}${startEvent}.err
    
echo "Logging output to " $LogFile
echo "Logging errors to " $ErrFile
echo "now submitting this script: "

echo sbatch --mem-per-cpu=32GB -q express -p erhip --nodes=4 -o $LogFile -e $ErrFile --time=01:00:00 --job-name=${jobName} -- ${ExecPath}/submit/qwrap.sh ${ExecPath} $execute $fileNo

sbatch --mem-per-cpu=32GB -q express -p erhip --nodes=4 -o $LogFile -e $ErrFile --time=01:00:00 --job-name=${jobName} -- ${ExecPath}/submit/qwrap.sh ${ExecPath} $execute $fileNo

@ startEvent = $startEvent + $nEventsPerFile
@ fileNo = $fileNo + 1

echo " "

end
