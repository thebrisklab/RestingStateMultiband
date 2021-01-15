#!/bin/bash 

## Controls which modules are being run:
# MODULE 1:
# motion correction, topup, and slice timing:
run_motion=0

# MODULE 2:
# anatomical registration: FNIRT to MNI
run_anatomical=0

# MODULE 3:
# Functional registration to MNI
run_functional_registration=0

# MODULE 4:
# Nuisance regression
run_nuisance_regression=1

## Modify for each subject:
subject=subjectID
FWHM=6
hp=0.009
lp=0.08

## Runs four piplines: temporal filtering (yes / no) by spatial smoothing (yes / no)

# Note: CSIC directory structure currently in use
# cwd=~/OptimalSMS_rsfMRI/Data/${subject}/nii
cwd=~/risk_share/OptimalSMS_rsfMRI/rsfmri_multiband_nii_dicom/${subject}/nii

anat_dir=$cwd
standard=/usr/local/fsl/data/standard/MNI152_T1_2mm_brain.nii.gz
standard_nomask=/usr/local/fsl/data/standard/MNI152_T1_2mm.nii.gz
# NOTE: Need to edit the paths in the T1_2_MNI152_2mm.cnf as well

prior_dir=/home/benjamin.risk/risk_share/OptimalSMS_rsfMRI/tissuepriors
# fsl priors thresholded >0.5

anat=t1_mprage_SAG

# afnidir=/Users/ben/abin
afnidir=/usr/local/AFNI/linux_xorg7_64

# fsldir=/Users/brisk/Applications2/fsl/bin/
fsldir=/usr/local/fsl/bin

filelistAll="rsfMRI_3_3mm_iso_MB1_AP rsfMRI_2mm_iso_MB1_AP rsfMRI_2mm_iso_MB2_AP rsfMRI_2mm_iso_MB3_AP rsfMRI_2mm_iso_MB4_AP rsfMRI_2mm_iso_MB6_AP rsfMRI_2mm_iso_MB8_AP rsfMRI_2mm_iso_MB9_AP rsfMRI_2mm_iso_MB12_AP"
filelistMB="rsfMRI_2mm_iso_MB2_AP rsfMRI_2mm_iso_MB3_AP rsfMRI_2mm_iso_MB4_AP rsfMRI_2mm_iso_MB6_AP rsfMRI_2mm_iso_MB8_AP rsfMRI_2mm_iso_MB9_AP rsfMRI_2mm_iso_MB12_AP"
filelist72="rsfMRI_2mm_iso_MB1_AP rsfMRI_2mm_iso_MB2_AP rsfMRI_2mm_iso_MB3_AP rsfMRI_2mm_iso_MB4_AP rsfMRI_2mm_iso_MB6_AP rsfMRI_2mm_iso_MB8_AP rsfMRI_2mm_iso_MB9_AP rsfMRI_2mm_iso_MB12_AP"

# slice timing file in same order as filelistAll:
slicetimefiles=(rsfMRI3.3mmisoMB1AP_slicetiming.txt rsfMRI2mmisoMB1AP_slicetiming.txt rsfMRI2mmisoMB2AP_slicetiming.txt rsfMRI2mmisoMB3AP_slicetiming.txt rsfMRI2mmisoMB4AP_slicetiming.txt rsfMRI2mmisoMB6AP_slicetiming.txt rsfMRI2mmisoMB8AP_slicetiming.txt rsfMRI2mmisoMB9AP_slicetiming.txt rsfMRI2mmisoMB12AP_slicetiming.txt)
n_volslist=(118 64 125 188 248 364 478 520 688)
TRlist=(3.0 5.67 2.85 1.91 1.44 0.962 0.736 0.675 0.512)

cd $cwd
## ######################
## 1. MOTION CORRECTION, TOP-UP, AND SLICE TIMING:
## ######################

if [ $run_motion -eq 1 ]
then
	echo ------------------------------
	echo ---- RUNNING FUNCTIONAL PREPROCESSING ----
	echo ------------------------------

	## NOTE: run <DicomConversion.m> to create nifti files before running this script
	index=0
	for rest in $filelistAll
	do
		n_vols=`echo ${n_volslist[$index]}-1 | bc`
		mkdir ./${rest}
		## 1. Dropping first 10 seconds of data

		echo "Dropping first TRs"
		if [ $index -eq 0 ] 
		then 
		$afnidir/3dcalc -a ${rest}.nii.gz[4..${n_vols}] -expr 'a' -prefix ./${rest}/${rest}_dr.nii.gz -overwrite
		elif [ $index -eq 1 ]
		then 
		$afnidir/3dcalc -a ${rest}.nii.gz[2..${n_vols}] -expr 'a' -prefix ./${rest}/${rest}_dr.nii.gz -overwrite
		elif [ $index -eq 2 ]
		then 
		$afnidir/3dcalc -a ${rest}.nii.gz[4..${n_vols}] -expr 'a' -prefix ./${rest}/${rest}_dr.nii.gz -overwrite
		elif [ $index -eq 3 ]
		then 
		$afnidir/3dcalc -a ${rest}.nii.gz[6..${n_vols}] -expr 'a' -prefix ./${rest}/${rest}_dr.nii.gz -overwrite
		elif [ $index -eq 4 ]
		then 
		$afnidir/3dcalc -a ${rest}.nii.gz[7..${n_vols}] -expr 'a' -prefix ./${rest}/${rest}_dr.nii.gz -overwrite
		elif [ $index -eq 5 ]
		then 
		$afnidir/3dcalc -a ${rest}.nii.gz[11..${n_vols}] -expr 'a' -prefix ./${rest}/${rest}_dr.nii.gz -overwrite
		elif [ $index -eq 6 ]
		then 
		$afnidir/3dcalc -a ${rest}.nii.gz[14..${n_vols}] -expr 'a' -prefix ./${rest}/${rest}_dr.nii.gz -overwrite
		elif [ $index -eq 7 ]
		then 
		$afnidir/3dcalc -a ${rest}.nii.gz[15..${n_vols}] -expr 'a' -prefix ./${rest}/${rest}_dr.nii.gz -overwrite
		elif [ $index -eq 8 ]
		then 
		$afnidir/3dcalc -a ${rest}.nii.gz[20..${n_vols}] -expr 'a' -prefix ./${rest}/${rest}_dr.nii.gz -overwrite
		fi

		index=$(($index+1))
	done


	## Extract twelfth volume from single-band images, corresponds to 8 in dropped version; this will serve as the "SBRef" for registration,
	## and facilitates the use of the same code for single and multiband
	$afnidir/3dcalc -a rsfMRI_3_3mm_iso_MB1_AP.nii.gz[11] -expr 'a' -prefix rsfMRI_3_3mm_iso_MB1_AP_SBRef.nii.gz -overwrite
	$afnidir/3dcalc -a rsfMRI_2mm_iso_MB1_AP.nii.gz[11] -expr 'a' -prefix rsfMRI_2mm_iso_MB1_AP_SBRef.nii.gz -overwrite



	## Motion correction
	## mcflirt: use -plots options to retain 6 parameters
	for rest in $filelistAll
	do

	## NOTE: Takes 15 hours to do all nine datasets with spline_final option
	##  $fsldir/mcflirt -in ${rest}.nii.gz -out ./${rest}/${rest}_mcf.nii.gz -refvol ${rest}_SBRef.nii.gz -plots -spline_final
	# Defaults to 6 dof, i.e., rigid body
	# use -rmsrel to obtain framewise displacement.
	$fsldir/mcflirt -in ./${rest}/${rest}_dr.nii.gz -out ./${rest}/${rest}_mcf.nii.gz -refvol ${rest}_SBRef.nii.gz -plots -rmsrel

	done

	## Estimate distortion field
	## Topup correction
	## Note: Alternatively, we could do a single motion correction to the SBREF3 image.
	## However, this would result in inflated estimates of subject motion using conventional measures
	echo -----------------------------------------------------
	echo -------- Applying TOPUP
	echo -----------------------------------------------------


	## motion correct PA to the AP reference:
	$fsldir/flirt -in rsfMRI_2mm_iso_MB3_PA_SBRef.nii.gz  -ref rsfMRI_2mm_iso_MB3_AP_SBRef.nii.gz -out rsfMRI_2mm_iso_MB3_PA_SBRef_mc.nii.gz -dof 6 -interp spline

	## fslmerge: concatenate the AP and PA single band reference images
	$fsldir/fslmerge -t rsfMRI_2mm_iso_MB3_AP_PA_SBRef.nii.gz  rsfMRI_2mm_iso_MB3_AP_SBRef.nii.gz rsfMRI_2mm_iso_MB3_PA_SBRef_mc.nii.gz 


	$fsldir/topup --imain=rsfMRI_2mm_iso_MB3_AP_PA_SBRef.nii.gz --datain=../../acquisition_parameters.txt --config=../../b02b0.cnf --out=my_output

	## register 6 dof sbref mb to sbref mb=3 AP, apply to volumes
	## This is also done to mb=3 for consistency
	## SBREF --> SBREF 3
	for rest in $filelist72
	do
	## $fsldir/flirt -ref rsfMRI_2mm_iso_MB3_AP_SBRef.nii.gz -in ${rest}_SBRef.nii.gz -out ./${rest}/${rest}_SBRef_toSBRef3.nii.gz -omat ./${rest}/${rest}_SBRef_toSBRef3.mat -dof 6 -interp spline
	$fsldir/flirt -ref rsfMRI_2mm_iso_MB3_AP_SBRef.nii.gz -in ${rest}_SBRef.nii.gz -out ./${rest}/${rest}_SBRef_toSBRef3.nii.gz -omat ./${rest}/${rest}_SBRef_toSBRef3.mat -dof 6


	## Register 4D to SBREF3 using previous transformation
	## $fsldir/flirt -in ${rest}/${rest}_mcf.nii.gz -ref rsfMRI_2mm_iso_MB3_AP_SBRef.nii.gz -out ${rest}/${rest}_mcf_tomb3.nii.gz -applyxfm -init ./${rest}/${rest}_SBRef_toSBRef3.mat -interp spline
	$fsldir/flirt -in ${rest}/${rest}_mcf.nii.gz -ref rsfMRI_2mm_iso_MB3_AP_SBRef.nii.gz -out ${rest}/${rest}_mcf_tomb3.nii.gz -applyxfm -init ./${rest}/${rest}_SBRef_toSBRef3.mat

	## topup defaults to spline
	$fsldir/applytopup --imain=./${rest}/${rest}_mcf_tomb3.nii.gz --inindex=1 --datain=../../acquisition_parameters.txt --topup=my_output --method=jac --out=./${rest}/${rest}_tu.nii.gz

	done


	## ###################################
	## Special treatment for 3.3
	## Apply topup AFTER slice timing correction
	## This is slice timing is applied in the 3.3 mm space,
	## followed by applytopup on the resample 2mm data
	## -------------------------->
	index=0
	rest=rsfMRI_3_3mm_iso_MB1_AP
	cd ${cwd}

	echo -----------------------------------------------------
	echo -------- Slice time correction for ${rest} ----------
	echo -------- Using ${slicetimefiles[$index]} ------------
	echo -----------------------------------------------------

	$fsldir/slicetimer -i ./${rest}/${rest}_dr.nii.gz -o ${rest}/${rest}_prest.nii.gz --tcustom=../${slicetimefiles[$index]}

	## Note: here we use spline interpolation because the resolution is changing
	$fsldir/flirt -ref rsfMRI_2mm_iso_MB3_AP_SBRef.nii.gz -in ${rest}_SBRef.nii.gz -out ./${rest}/${rest}_SBRef_toSBRef3.nii.gz -omat ./${rest}/${rest}_SBRef_toSBRef3.mat -dof 6 -interp spline

	## Register 4D to SBREF3 using previous transformation
	## For 3.3 mm and MB=1, use slice time corrected data:
	$fsldir/flirt -in ${rest}/${rest}_prest.nii.gz -ref rsfMRI_2mm_iso_MB3_AP_SBRef.nii.gz -out ${rest}/${rest}_mcf_tomb3.nii.gz -applyxfm -init ./${rest}/${rest}_SBRef_toSBRef3.mat -interp spline

	## topup defaults to spline
	$fsldir/applytopup --imain=./${rest}/${rest}_mcf_tomb3.nii.gz --inindex=1 --datain=../../acquisition_parameters.txt --topup=my_output --method=jac --out=./${rest}/${rest}_tu.nii.gz
	## $fsldir/applytopup --imain=./${rest}/${rest}_mcf_tomb3.nii.gz --inindex=1 --datain=../../acquisition_parameters.txt --topup=my_output --method=jac --out=./${rest}/${rest}_tu.nii.gz --interp=trilinear


	cp ./${rest}/${rest}_tu.nii.gz ./${rest}/${rest}_st.nii.gz

	## END SPECIAL TREATMENT FOR 3.3  #
	##<--------------------------


	index=1
	for rest in $filelist72
	do
	cd ${cwd}/${rest}

	## Slice timing correction
	## NOTE: --repetition flag does not do anything when specifying custom slice timing acquisition
	echo -----------------------------------------------------
	echo -------- Slice time correction for ${rest} ----------
	echo -------- Using ${slicetimefiles[$index]} ------------
	echo -----------------------------------------------------

	$fsldir/slicetimer -i ${rest}_tu.nii.gz -o ${rest}_st.nii.gz --tcustom=../../${slicetimefiles[$index]}
	index=$(($index+1))
	done

fi









## ######################
## 2. Structural registration
## ######################

if [ $run_anatomical -eq 1 ]
then
	echo ------------------------------
	echo ---- RUNNING ANATOMICAL ----
	echo ------------------------------

	## Skull-strip anatomical:
	cd $cwd

	## rm ${anat}_surf.nii.gz* ${anat}_brain.nii.gz
	$afnidir/3dSkullStrip -input ${anat}.nii.gz -o_ply ${anat}_surf.nii.gz -overwrite
	$afnidir/3dcalc -a ${anat}.nii.gz -b ${anat}_surf.nii.gz -expr 'a*step(b)' -prefix ${anat}_brain.nii.gz -overwrite
	## NOTE:  bet grabbed part of the throat, but didn't grab the marrow, while afni skullstrip grabbed 
	## some marrow, but not the throat. Prefer afni.
	## $fsldir/bet ${anat}.nii.gz ${anat}_brain_bet.nii.gz

	## 4. T1->STANDARD
	## NOTE: flirt works better with BET while fnirt works better without bet

	## $fsldir/flirt -ref ${standard} -in ${anat_dir}/${anat}_brain.nii.gz -out highres2standard -omat highres2standard.mat -cost corratio -searchcost corratio -dof 12 -interp spline
	$fsldir/flirt -ref ${standard} -in ${anat_dir}/${anat}_brain.nii.gz -out highres2standard -omat highres2standard.mat -cost corratio -searchcost corratio -dof 12

	## $fsldir/fnirt --in=${anat_dir}/${anat} --iout=highres2standardnl --aff=highres2standard.mat --interp=spline --config=../../T1_2_MNI152_2mm.cnf 
	$fsldir/fnirt --in=${anat_dir}/${anat} --iout=highres2standardnl --aff=highres2standard.mat --config=../../T1_2_MNI152_2mm.cnf 



	## Create mat file for conversion from standard to high res
#	$fsldir/convert_xfm -inverse -omat standard2highres.mat highres2standard.mat
fi




## ######################
## 3. Functional Registration
## ######################
if [ $run_functional_registration -eq 1 ]
then
	echo ------------------------------
	echo !!!! RUNNING REGISTRATION !!!!
	echo ------------------------------


	## Register FUNCTION MB3 to T1. Note all MB data have already been registered
	## to MB3.
	reg_dir=${cwd}/rsfMRI_2mm_iso_MB3_AP/reg
	mkdir ${reg_dir}
	cd $cwd

	## create example_func for compatability with pipeline:
	## NOTE: data have been aligned to MB3_SBRef
	## Thus, use MB3_SBRef for all 4D datasets:
	## cp ${rest}_SBRef.nii.gz ./${rest}/example_func.nii.gz
	cp rsfMRI_2mm_iso_MB3_AP_SBRef.nii.gz ${reg_dir}/example_func.nii.gz


	## 1. Copy required images into reg directory
	### copy anatomical
	cp ${anat_dir}/${anat}_brain.nii.gz ${reg_dir}/highres.nii.gz
	### copy standard
	cp ${standard} ${reg_dir}/standard.nii.gz

	cp ${cwd}/highres2standard.mat ${reg_dir}/highres2standard.mat



	## 2. cd into reg directory
	cd ${reg_dir}

	## 3. FUNC->T1

	## Need this for the initial transformation for applywarp:
	$fsldir/flirt -ref highres -in example_func -out example_func2highres -omat example_func2highres.mat -cost corratio -dof 6 

	## 5. FUNC->STANDARD
	## Create mat file for registration of functional to standard
	$fsldir/convert_xfm -omat example_func2standard.mat -concat highres2standard.mat example_func2highres.mat

	for rest in $filelistAll
	do

	func_dir=${cwd}/${rest}

	## reg_dir=${cwd}/${rest}/reg
	## mkdir ${reg_dir}

	## Applywarp defaults to trilinear
	##$fsldir/applywarp --ref=${standard_nomask} --in=${func_dir}/${rest}_rpp --out=${func_dir}/${rest}_res2standard_fnirt --warp=${anat_dir}/${anat}_warpcoef.nii --premat=${reg_dir}/example_func2highres.mat --interp=spline
	$fsldir/applywarp --ref=${standard_nomask} --in=${func_dir}/${rest}_st --out=${func_dir}/${rest}_res2standard_fnirt --warp=${anat_dir}/${anat}_warpcoef.nii --premat=${reg_dir}/example_func2highres.mat
	$afnidir/3dTstat -mean -prefix ${func_dir}/${rest}_res2standard_fnirt_mean.nii.gz ${func_dir}/${rest}_res2standard_fnirt.nii.gz  -overwrite

	done

	############################
	## Basic preprocessing: 
	## Skull strip
	############################


	for rest in $filelistAll
	do
	func_dir=${cwd}/${rest}
	cd ${func_dir}
	## Remove skull/edge detect
	echo "Skull stripping ${subject} ${rest}"
	$afnidir/3dAutomask -prefix ${rest}_mask.nii.gz -dilate 1 ${rest}_res2standard_fnirt.nii.gz -overwrite
	$afnidir/3dcalc -a ${rest}_res2standard_fnirt.nii.gz -b ${rest}_mask.nii.gz -expr 'a*b' -prefix ${rest}_ss.nii.gz -overwrite

	## Grandmean scaling
	echo "Grand-mean scaling ${subject} ${rest}"
	$fsldir/fslmaths ${rest}_ss.nii.gz -ing 10000 ${rest}_pp.nii.gz -odt float

	# Generate binary mask, this is slightly bigger than 3dAutomask for an MB 12 that I checked
	echo "Generating mask of preprocessed data for ${rest}"
	$fsldir/fslmaths ${rest}_pp.nii.gz -Tmin -bin ${rest}_pp_mask.nii.gz -odt char

	## BRisk: Clean-up
	## rm ${rest}_st.nii.gz ${rest}_ss.nii.gz ${rest}_sm.nii.gz ${rest}_gms.nii.gz ${rest}_filt.nii.gz ${rest}_filt_mean.nii.gz ${rest}_dt.nii.gz ${rest}_mcf_tomb3.nii.gz ${rest}_mcf.nii.gz ${rest}_mask.nii.gz
	done
fi

##########################################################################################################################
##---NUISANCE SIGNAL REGRESSION ----------------------------------------------------------------------------------------------------##
##########################################################################################################################

## BRisk: The F1000 pipelines uses segmentation and then intersects with FSL templates.
## The fsl templates have narrowly defined wm and csf using the mask (_bin) >=0.51. The Segmentation in fast captures a 
## lot of csf between the gray matter and skull, and is a little tricky to use due to partial volume effects. 
## Since the templates are conservative, and since we are using FNIRT, I do the averaging in 
## MNI space aligned volumes using the FSL (harvard) tissue priors

if [ $run_nuisance_regression -eq 1 ]
then


	echo --------------------------------------------
	echo !!!! RUNNING NUISANCE SIGNAL REGRESSION !!!!
	echo --------------------------------------------

	segment_dir=${cwd}/segment
	reg_dir=${cwd}/rsfMRI_2mm_iso_MB3_AP/reg
	func_dir=${cwd}/rsfMRI_2mm_iso_MB3_AP

	# First, intersect tissue priors with global mask from the MB3 acquisition
	rest=rsfMRI_2mm_iso_MB3_AP
	echo ${prior_dir}

	PRIOR_WHITE=${prior_dir}/avg152T1_white_bin.nii.gz
	PRIOR_CSF=${prior_dir}/avg152T1_csf_bin.nii.gz

	mkdir ${segment_dir}

	## 4. Copy functional mask from FSLpreproc step 5 - this is the global signal mask
	echo "Creating global mask"
	cp ${func_dir}/${rest}_pp_mask.nii.gz ${segment_dir}/global_mask.nii.gz

	## Mask CSF template by subject's global mask and copy to segment directory:
	$fsldir/fslmaths ${PRIOR_CSF} -mas ${segment_dir}/global_mask ${segment_dir}/csf_mask

	## White matter mask
	$fsldir/fslmaths ${PRIOR_WHITE} -mas ${segment_dir}/global_mask ${segment_dir}/wm_mask

	index=0
	for rest in $filelistAll
	do
	func_dir=${cwd}/${rest}
	nuisance_dir=${func_dir}/nuisance

	## 1. make nuisance directory
	mkdir -p ${nuisance_dir}

	# Extract signal for global, csf, and wm
	## 2. Global
	echo "Extracting global signal for ${subject}"
	$afnidir/3dmaskave -mask ${segment_dir}/global_mask.nii.gz -quiet ${func_dir}/${rest}_pp.nii.gz > ${nuisance_dir}/global.1D -overwrite

	## 3. csf
	echo "Extracting signal from csf for ${subject}"
	$afnidir/3dmaskave -mask ${segment_dir}/csf_mask.nii.gz -quiet ${func_dir}/${rest}_pp.nii.gz > ${nuisance_dir}/csf.1D -overwrite

	## 4. wm
	echo "Extracting signal from white matter for ${subject}"
	$afnidir/3dmaskave -mask ${segment_dir}/wm_mask.nii.gz -quiet ${func_dir}/${rest}_pp.nii.gz > ${nuisance_dir}/wm.1D -overwrite

	## 5. Perform linear and quadratic detrending, nuisance regression, temporal filtering
	echo "Performing nuisance regression, spatial smoothing, and temporal filtering"

## No spatial smoothing:
	$afnidir/3dTproject -input ${func_dir}/${rest}_pp.nii.gz -prefix ${func_dir}/${rest}_pp_9p.nii.gz -polort 2 -ort ${nuisance_dir}/global.1D ${nuisance_dir}/wm.1D ${nuisance_dir}/csf.1D ${func_dir}/${rest}_mcf.nii.gz.par -automask -overwrite -verb
	$afnidir/3dTproject -input ${func_dir}/${rest}_pp.nii.gz -prefix ${func_dir}/${rest}_pp_9p_tf.nii.gz -polort 2 -ort ${nuisance_dir}/global.1D ${nuisance_dir}/wm.1D ${nuisance_dir}/csf.1D ${func_dir}/${rest}_mcf.nii.gz.par -passband $hp $lp -automask -overwrite -verb

## sdev maps
	$afnidir/3dTstat -stdev -prefix ${func_dir}/${rest}_pp_stdev.nii.gz ${func_dir}/${rest}_pp.nii.gz  -overwrite

	$afnidir/3dcalc -a ${func_dir}/${rest}_pp_9p_stdev.nii.gz -b  $base_data -expr 'b/a' -prefix ${func_dir}/${rest}_pp_9p_gfactor.nii.gz -overwrite
	$afnidir/3dcalc -a ${func_dir}/${rest}_pp_9p_tf_stdev.nii.gz -b  $base_data_tf -expr 'b/a' -prefix ${func_dir}/${rest}_pp_9p_tf_gfactor.nii.gz -overwrite


## Spatial smoothing:
	$afnidir/3dTproject -input ${func_dir}/${rest}_pp.nii.gz -prefix ${func_dir}/${rest}_pp_9p_sm.nii.gz -polort 2 -ort ${nuisance_dir}/global.1D ${nuisance_dir}/wm.1D ${nuisance_dir}/csf.1D ${func_dir}/${rest}_mcf.nii.gz.par -blur $FWHM -automask -overwrite -verb
	$afnidir/3dTproject -input ${func_dir}/${rest}_pp.nii.gz -prefix ${func_dir}/${rest}_pp_9p_tf_sm.nii.gz -polort 2 -ort ${nuisance_dir}/global.1D ${nuisance_dir}/wm.1D ${nuisance_dir}/csf.1D ${func_dir}/${rest}_mcf.nii.gz.par -passband $hp $lp -blur $FWHM -automask -overwrite -verb

## sdev maps
	$afnidir/3dTstat -stdev -prefix ${func_dir}/${rest}_pp_9p_sm_stdev.nii.gz ${func_dir}/${rest}_pp_9p_sm.nii.gz  -overwrite
	$afnidir/3dTstat -stdev -prefix ${func_dir}/${rest}_pp_9p_tf_sm_stdev.nii.gz ${func_dir}/${rest}_pp_9p_tf_sm.nii.gz -overwrite

	done

fi
