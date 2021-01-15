#!/bin/bash

subjectList="102954 118095 123821 134116 144389 164872 184537 193775 111581 118393 125506 138945 145064 167212 191985 113214 122359 128751 139050 145235 171880 197354 116716 123779 128991 141973 149233 149235 180653 193773 195590 193728"

# These participants have no MB9: 195590, 193728. Will create errors for theMB 9 run for these participants (but it's fine).

for subjectID in $subjectList
do

cd ~/risk_share/OptimalSMS_rsfMRI/rsfmri_multiband_nii_dicom/${subjectID}
sed -e s:subjectID:"${subjectID}":g <~/risk_share/OptimalSMS_rsfMRI/Programs/Voxel2nodesystem_singlesubj.m >Voxel2nodesystem_${subjectID}.m

sed -e s:subjectID:"${subjectID}":g <~/risk_share/OptimalSMS_rsfMRI/Programs/Voxel2nodesystem_singlesubj.sh >Voxel2nodesystem_${subjectID}.sh


qsub  -cwd -N vox2node_${subjectID} Voxel2nodesystem_${subjectID}.sh

done
