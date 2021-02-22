# RestingStateMultiband

Scripts supporting "Which multiband factor should you choose for your resting-state fMRI study?"

A) Shell scripts. These must be run first:

1. `preprocessing_loopSubjects.sh` loops through the subjects and creates a copy of 
	`preprocessing_singlesubj_3dTproject.sh` for each subject in which the "subjectID" is replaced with the numeric value. It then submits (via qsub) the preprocessing scripts for all subjects. The script performs preprocessing for four pipelines: without / with temporal filtering x without / with spatial smoothing. The main preprocessing script utlizes elements of the Functional Connectomes 1000 project scripts (https://www.nitrc.org/projects/fcon_1000/) but with a number of updates including topup distortion correction, fnirt, and simultaneous filtering and nuisance regression with AFNI 3dTproject. This takes approximately four hours (assuming the number of nodes=number of subjects).

5. `Voxel2nodesystem_loopSubjects.sh` calculates the average time course for each power ROI. Loops across subjects. 

	Creates files
	`Voxel2nodesystem_singlesubj.sh` and 
	`Voxel2nodesystem_singlesubj.m` where the former is used to qsub the latter. 
	Also calculates the subject correlation matrices. Takes about 15 minutes.


B) R and Matlab Scripts: these were run interactively

`1_MotionExamine_cluster.R`: flags scans that failed motion control

`2_CreateAveSDMaps.R`: creates standard deviation maps

`3_GFactorTables.R` create g-factor tables
 
`4_run_SeedMaps.m` calls the function SeedMaps_create.m, creating seed maps for the specified ROIs and each pipeline

`5_ConsolidateCorrMat_cluster.m`  gathers the subject correlation matrices into a single array. Also calcualtes the g*-factors. 

`6a_CorrelationsVGfactorsAndMB_GAMsGEEs.R`
 Analyzes impacts of g-factor and multiband factor for a subset of intra-commmunity edges (6 regions/communities)
 1. GAMMs analyzing impact of gfactor on correlations
 2. GEEs analyzing subset of edges for different communities and MB factor
 Creates figures "Z correlation versus g*-factor" from gamm


`6b_CorrelationsVGfactorsAndMB_GAMsGEEs_ss.R`
	same as 6a but run on the datasets with 6-mm smoothing

`7a_ExamineCorrelationMatrices_plots.R`
-creates plots of MB factors versus SB 2 mm colored by g* factor
-creates plots of correlation matrices
-creates the dataset somemb.z.cor.srt.RData for input to the permutation tests

`7b_ExamineCorrelationMatrices_plots_sm.R`
	same as 5a but run on the datasets with 6-mm smoothing

C) Permutation tests were run on the cluster using qsub, e.g., 
```
cd ~/risk_share/OptimalSMS_rsfMRI/Programs
qsub -cwd -N permtest_notf a6a_CompareMBFactor_PermutationTets_notf.sh
```

`a8a_CompareMBFactors_PermutationTests_notf.sh` simple shell script to execute using the scheduler and qsub

`8a_CompareMBFactors_PermutationTests_notf.R` permutation tests comparing proportion activated in MB 8 to other acquisitions
Takes approximately 2.25 hours to run

`a8a_CompareMBFactors_PermutationTests_notf.sh`

`8b_CompareMBFactors_PermutationTests_tf.R`

`a8c_CompareMBFactors_PermutationTests_notf_vs_tf`

`8c_CompareMBFactors_PermutationTests_notf_versus_tf`
Permutation tests to support the hypothesis that temporal filtering decreases the proportion active

`a8d_CompareMBFactors_PermutationTests_notf_sm.sh` simple shell script to execute using the scheduler and qsub

`8d_CompareMBFactors_PermutationTests_notf_sm.R` identical to 6a but using smoothed data


D) This script is run interactively:

`9_PermutationTests_tables.R`  creates tables summarizing the permutation tests from 8a-8d. 



