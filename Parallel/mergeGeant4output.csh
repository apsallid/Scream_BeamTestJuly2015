#!/bin/tcsh

setenv suddirineoslistguns "pi- mu- e-"
#energies in GeV
setenv suddirineoslistenergies "50 100 150 200"

foreach subdirineosguns  ($suddirineoslistguns)

echo "===================================================================================="

foreach subdirineosenergies  ($suddirineoslistenergies)

echo "------------------------"
echo "${subdirineosguns} ${subdirineosenergies} GeV"

setenv workpath "/afs/cern.ch/work/a/apsallid/CMS/Geant4/Scream/Parallel/$subdirineosguns/$subdirineosenergies/results"

setenv output testcalo_${subdirineosguns}_${subdirineosenergies}.root
setenv outputtree testcalo_tree_${subdirineosguns}_${subdirineosenergies}.root

hadd $output ${workpath}/geantoutput_${subdirineosguns}_${subdirineosenergies}GeV_*.root 

hadd $outputtree ${workpath}/geantoutput_tree_${subdirineosguns}_${subdirineosenergies}GeV_*.root

end 

end

foreach subdirineosguns  ($suddirineoslistguns)

echo "===================================================================================="

foreach subdirineosenergies  ($suddirineoslistenergies)

echo "------------------------"
echo "${subdirineosguns} ${subdirineosenergies} GeV"

foreach num (`seq 2 100`)

setenv workpathadc "/afs/cern.ch/work/a/apsallid/CMS/Geant4/Scream/Parallel/$subdirineosguns/$subdirineosenergies/results"

setenv filewithadc ${workpathadc}/adc_${subdirineosguns}_${subdirineosenergies}GeV_${num}.data

#cat $filewithadc >> ${workpathadc}/adc_${subdirineosguns}_${subdirineosenergies}GeV_1.data

end

mv ${workpathadc}/adc_${subdirineosguns}_${subdirineosenergies}GeV_1.data adc_${subdirineosguns}_${subdirineosenergies}GeV.data

end

end


