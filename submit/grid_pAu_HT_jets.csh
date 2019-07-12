#!/bin/csh

# used to submit sequential jobs on the grid

# first make sure program is updated and exists
 make bin/pAu_HT_jets || exit

set ExecPath = `pwd`
set execute = './bin/pAu_HT_jets'
set numevents = 1000
set base = /wsu/home/el/el98/el9852/physics/analysis/pAu_analysis/production_pAu200_2015/HT/pAu_2015_200_HT
set BackgroundChargeBias = allBG
set JetChargeBias = allJets
set outFile = HT_JP2

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
set OutBase = "$OutBase$uscore$outFile$uscore$uscore$BackgroundChargeBias$uscore$JetChargeBias"
    
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
    
set arg = "$Files $outLocation$outName $numevents $BackgroundChargeBias $JetChargeBias"

echo "now submitting this script: "
echo qsub -V -l mem=4GB -o $LogFile -e $ErrFile -N $1 -- ${ExecPath}/submit/qwrap.sh ${ExecPath} $execute $arg
    

qsub -V -q wsuq -l mem=4GB -o $LogFile -e $ErrFile -N pAu_analysis -- ${ExecPath}/submit/qwrap.sh ${ExecPath} $execute $arg

@ i++

end
