#!/bin/csh

# used to submit sequential jobs on the grid

# first make sure program is updated and exists
 make bin/sim || exit

set ExecPath = `pwd`
set execute = './bin/sim'
set numevents = -1
#set base = /tier2/home/groups/rhi/STAR/Data/EmbedPythiaRun15pAu200_picos/P18ih/pAu_200_REREproduction_2015/out/HT2/pt-hat
set base = /tier2/home/groups/rhi/STAR/Data/EmbedPythiaRun15pAu200_picos/P18ih/pAu_200_production_2015/out/pt-hat
set outDir = sim
set outFile = pAu2015embedding

# Create the folder name for output
#set outFile = stock
# Make the directories since they may not exist...                                                                                                                             
if ( ! -d out/${outDir} ) then
mkdir -p out/${outDir}
endif

if ( ! -d log/${outDir} ) then
mkdir -p log/${outDir}
endif

# Now Submit jobs for each data file                                                                                                                                           
foreach input ( ${base}* )

# Create the output file base name                                                                                                                                             
set OutBase = `basename $input | sed 's/.root//g'`
set uscore = "_"
set OutBase = "$outFile$uscore$OutBase"

# Make the output names and path                                                                                                                                               
set outLocation = out/${outDir}/
set outName = ${OutBase}.root

# Input files                                                                                                                                                                  
set Files = ${input}

# Logfiles. Thanks cshell for this "elegant" syntax to split err and out                                                                                                       
set LogFile     = log/${outDir}/${OutBase}.log
set ErrFile     = log/${outDir}/${OutBase}.err
set matchFlag    = "match"

echo "Logging output to " $LogFile
echo "Logging errors to " $ErrFile
    
    # //! [0]: output directory
    # //! [1]: output filename
    # //! [2]: flag: require match? Options: "nomatch" (match = 0), or "match" (match = 1).
    # //! [3]: input data

set arg = "$outLocation$outName  $matchFlag $Files"

echo "now submitting this script: "
echo sbatch --mem-per-cpu=8GB -q express -p erhip -o $LogFile -e $ErrFile -t 45 --job-name=sim $1 -- ${ExecPath}/submit/qwrap.sh ${ExecPath} $execute $arg
    

#qsub -V -q erhiq -l mem=4GB -o $LogFile -e $ErrFile -N $1 -- ${ExecPath}/submit/qwrap.sh ${ExecPath} $execute $arg
sbatch --mem-per-cpu=8GB -q express -p erhip -o $LogFile -e $ErrFile -t 45 --job-name=sim $1 -- ${ExecPath}/submit/qwrap.sh ${ExecPath} $execute $arg
    
end
