#!/bin/csh

# used to submit sequential jobs on the grid

# first make sure program is updated and exists
 make bin/pAu_analysis || exit

set ExecPath = `pwd`
set execute = './bin/pAu_analysis'
set numevents = 1000000
set base = /wsu/home/el/el98/el9852/physics/analysis/pAu_analysis/production_pAu200_2015/MB/pAu_2015_200_MB
set outFile = MB

# Create the folder name for output
#set outFile = stock
# Make the directories since they may not exist...                                                                                                                             
if ( ! -d out/${outFile} ) then
mkdir -p out/${outFile}
endif

if ( ! -d log/${outFile} ) then
mkdir -p log/${outFile}
endif

#echo ${base}
# Now Submit jobs for each data file                                                                                                                                           
foreach input ( ${base}* )

# Create the output file base name                                                                                                                                             
set OutBase = `basename $input | sed 's/.root//g'`
set uscore = "_"
set OutBase = "$OutBase$uscore$outFile"
    
# Make the output names and path                                                                                                                                               
set outLocation = out/${outFile}/
set outName = ${OutBase}.root

# Input files                                                                                                                                                                  
set Files = ${input}

# Logfiles. Thanks cshell for this "elegant" syntax to split err and out                                                                                                       
set LogFile     = log/${outFile}/${OutBase}.log
set ErrFile     = log/${outFile}/${OutBase}.err

echo "Logging output to " $LogFile
echo "Logging errors to " $ErrFile
    
set arg = "$Files $outLocation$outName $numevents"

echo "now submitting this script: "
echo qsub -V -l mem=4GB -o $LogFile -e $ErrFile -N $1 -- ${ExecPath}/submit/qwrap.sh ${ExecPath} $execute $arg
    

qsub -V -q erhiq -l mem=4GB -o $LogFile -e $ErrFile -N pAu_analysis -- ${ExecPath}/submit/qwrap.sh ${ExecPath} $execute $arg

@ i++

end
