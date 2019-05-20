#!/bin/csh

# used to submit sequential jobs on the grid

# first make sure program is updated and exists
 make bin/pAu_analysis || exit

set ExecPath = `pwd`
set analysis = $1
set execute = './bin/pAu_analysis'

# Now Submit jobs for each data file
set i = 0
while ( $i < 27 )

# Create the output file base name
set OutBase = out/outfile_${i}

# Make the output names and path

set outName = ${OutBase}.root


# Logfiles. Thanks cshell for this "elegant" syntax to split err and out
set LogFile     = log/pAu_analysis_${i}.log
set ErrFile     = log/pAu_analysis_${i}.err

echo "Logging output to " $LogFile
echo "Logging errors to " $ErrFile

#set arg = "$outName"

qsub -V -q erhiq -l mem=2GB -o $LogFile -e $ErrFile -N jetfinderAnalysis -- ${ExecPath}/submit/qwrap.sh ${ExecPath} $execute

@ i++

end
