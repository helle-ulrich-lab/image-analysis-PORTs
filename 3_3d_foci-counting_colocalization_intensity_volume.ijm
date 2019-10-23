//--------------------------------------
// This script performs an object-based colocalization analysis for two
// target proteins imaged in two different channels. It determines the
// number of foci per nucleus for both targets, their intensities and
// volumes as well as the degree of overlap with the other target. The
// background fluorescence of one of the targets is used to create
// nuclear masks; if this is not possible, a fluorescent nuclear marker
// should be introduced by other means.
// author: Ronald Wong
// Acknowledgement: Fabrice P Cordelières. Part of the code is taken 
// from EMBO Bioimage Data Analysis course 2017
//-------------------------------------

//Initialize
setBatchMode(true);
print("\\Clear");
run("Colors...", "foreground=white background=black selection=magenta");
roiManager("reset");

//Create tables to write data
run("New... ", "name=Volume_Int_Colocalization_Channel1_combined type=Table");
print("[Volume_Int_Colocalization_Channel1_combined]","\\Headings: \tLabel\tFull object\tCommon part\tRatio\tObjectInt\tObjectMeanInt\tstrain\ttreatment\tpoint\treplicate\tnucleus\ttime frame");
run("New... ", "name=Volume_Int_Colocalization_Channel2_combined type=Table");
print("[Volume_Int_Colocalization_Channel2_combined]","\\Headings: \tLabel\tFull object\tCommon part\tRatio\tObjectInt\tObjectMeanInt\tstrain\ttreatment\tpoint\treplicate\tnucleus\ttime frame");
run("New... ", "name=Number_of_foci_Channel1_combined type=Table");
print("[Number_of_foci_Channel1_combined]","\\Headings: \tFoci count\tstrain\ttreatment\tpoint\treplicate\tnucleus\ttime frame");
run("New... ", "name=Number_of_foci_Channel2_combined type=Table");
print("[Number_of_foci_Channel2_combined]","\\Headings: \tFoci count\tstrain\ttreatment\tpoint\treplicate\tnucleus\ttime frame");

//create a dialog for user to input paramters
Dialog.create("Colocalization");
Dialog.addNumber("Minimum size of objects on channel1 (in voxels)", 5);
Dialog.addNumber("Minimum size of objects on channel2 (in voxels)", 5);
Dialog.addNumber("Lower threshold value for channel1", 2000);
Dialog.addNumber("Lower threshold value for channel2", 2600);
Dialog.show();

minsize_Channel1=Dialog.getNumber();
minsize_Channel2=Dialog.getNumber();
threshold_Channel1=Dialog.getNumber();
threshold_Channel2=Dialog.getNumber();

//get directory with image files
dir=getDirectory("Choose directory with the Image files");
print("processing directory ",dir);	
fileList=getFileList(dir);
nFiles=lengthOf(fileList);
print(nFiles+" files found");

for(k=0;k<nFiles;k++) {  //loop over all files in the directory
	if (matches(fileList[k],".*R3D_D3D.dv")) { //analyze deconvolved DelataVision files
		print("Processing "+fileList[k]);
		run("Bio-Formats Importer", "open="+dir+fileList[k]+" color_mode=Colorized view=Hyperstack stack_order=XYCZT");
		//get name of file to write to the table
		m=indexOf(fileList[k],"_");
		o=indexOf(fileList[k],"_", m+1);
		p=indexOf(fileList[k],"_", o+1);
		q=indexOf(fileList[k],"_", p+1);
		r=indexOf(fileList[k],"_", q+1);
		s=indexOf(fileList[k],"_", r+1);
		t=indexOf(fileList[k],"_", s+1);

		date=substring(fileList[k],m+1,o);
		strain=substring(fileList[k],o+1,p);
		treatment=substring(fileList[k],p+1,q);
		timepoint=substring(fileList[k],q+1,r);
		point=substring(fileList[k],s+1,t);
		point_num=parseFloat(point);
		replicate=point_num%3+1;
		
		//crop image
		makeRectangle(32, 32, 960, 960);
		run("Crop");
		ori=getTitle();
		
		//nuclear segmentation based on the nuclear signals from a nuclear protein
		selectWindow(ori);	
		run("Duplicate...", "title=[nucleus count] duplicate channels=2");
		run("Z Project...", "projection=[Average Intensity]");
		rename("SUM_nucleus count");
		run("Gaussian Blur...", "sigma=3");
		setAutoThreshold("Moments dark");
		run("Convert to Mask");
		run("Dilate");
		run("Analyze Particles...", "size=3-17 clear display exclude");				
		nucleus_count=nResults;
		run("Clear Results");
		print("nucleus count",nucleus_count);
		print("roimanager count",roiManager("count"));
		selectWindow("nucleus count");
		close();
		print(nucleus_count+"nuclei found in time point "+timepoint);
		//Each nucleus identified will be analyzed individually
		for (n=0; n<nucleus_count; n++){
			selectWindow("SUM_nucleus count");
			run("Select None");			
			run("Analyze Particles...", "size=3-17 clear display exclude add");
			print(n);
			print("nucleus count",nResults);
			print("roimanager count",roiManager("count"));
			selectWindow(ori);
			roiManager("select", n);
			run("Duplicate...", "title="+ori+"-nucleus-"+n+" duplicate");
					
			//Normalize names		
			run("Split Channels");
			selectWindow("C1-"+ori+"-nucleus-"+n);
			rename("Channel1");
			selectWindow("C2-"+ori+"-nucleus-"+n);
			rename("Channel2");	
				
			//Isolate Objects in 3D
			selectWindow("Channel1");
			run("Duplicate...", "title=Channel1_threshold duplicate");
			setThreshold(threshold_Channel1, 70000);
			run("Convert to Mask", "method=Default background=Dark black");
			run("3D OC Options", "  dots_size=5 font_size=10 redirect_to=none");
			run("3D Objects Counter", "threshold=0 slice=1 min.="+minsize_Channel1+" max.=100000 objects");
			rename("Tagged_map_Channel1");
				
			selectWindow("Channel2");
			run("Duplicate...", "title=Channel2_threshold duplicate");
			setThreshold(threshold_Channel2, 70000);
			run("Convert to Mask", "method=Default background=Dark black");
			run("3D OC Options", "  dots_size=5 font_size=10 redirect_to=none");
			run("3D Objects Counter", "threshold=0 slice=1 min.="+minsize_Channel2+" max.=100000 objects");
			rename("Tagged_map_Channel2");
		
			selectWindow("Tagged_map_Channel1");
			run("Duplicate...", "title=Mask_Channel1 duplicate");
			setThreshold(1, 255);
			run("Convert to Mask", "method=Default background=Dark");
			run("Divide...", "value=255 stack");
			resetMinAndMax();
				
			selectWindow("Tagged_map_Channel2");
			run("Duplicate...", "title=Mask_Channel2 duplicate");
			setThreshold(1, 255);
			run("Convert to Mask", "method=Default background=Dark");
			run("Divide...", "value=255 stack");
			resetMinAndMax();
				
			imageCalculator("AND create stack", "Mask_Channel1","Mask_Channel2");
			rename("Common_volumes");
				
			//Get colocalization values,intensities and volumes of foci
			objectsVolume1=getValues("Tagged_map_Channel1", "Mask_Channel1");
			commonVolume1=getValues("Tagged_map_Channel1", "Common_volumes");
			objectsInt1=getIntensity("Tagged_map_Channel1", "Channel1");
					
			generateOutputs(objectsVolume1, commonVolume1, objectsInt1, "Tagged_map_Channel1");
			if (nResults>0) {				
				String.copyResults();
				print("[Volume_Int_Colocalization_Channel1_combined]",String.paste);
				run("Clear Results");
			} else	print("[Volume_Int_Colocalization_Channel1_combined]","NaN\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\t"+strain+"\t"+treatment+"\t"+point+"\t"+replicate+"\t"+n+"\t"+timepoint);
			selectWindow("Coloc_Map");
			rename("Volume_Colocalization_Channel1");
									
			objectsVolume2=getValues("Tagged_map_Channel2", "Mask_Channel2");
			commonVolume2=getValues("Tagged_map_Channel2", "Common_volumes");
			objectsInt2=getIntensity("Tagged_map_Channel2", "Channel2");
									
			generateOutputs(objectsVolume2, commonVolume2, objectsInt2, "Tagged_map_Channel2");
			if (nResults>0) {					
				String.copyResults();
				print("[Volume_Int_Colocalization_Channel2_combined]",String.paste);
				run("Clear Results");
			} else	print("[Volume_Int_Colocalization_Channel2_combined]","NaN\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\t"+strain+"\t"+treatment+"\t"+point+"\t"+replicate+"\t"+n+"\t"+timepoint);
			selectWindow("Coloc_Map");
			rename("Volume_Colocalization_Channel2");
				
			//Get number of foci
			getFociCount("Tagged_map_Channel1");
			String.copyResults();
			print("[Number_of_foci_Channel1_combined]",String.paste);
			run("Clear Results");

			getFociCount("Tagged_map_Channel2");
			String.copyResults();
			print("[Number_of_foci_Channel2_combined]",String.paste);
			run("Clear Results");

			//close windows to prepare for next analysis
			selectWindow("Channel1");
			close();
			selectWindow("Channel2");
			close();
			selectWindow("Channel1_threshold");
			close();
			selectWindow("Channel2_threshold");
			close();
			selectWindow("Tagged_map_Channel1");
			close();
			selectWindow("Tagged_map_Channel2");
			close();
			selectWindow("Mask_Channel1");
			close();
			selectWindow("Mask_Channel2");
			close();
			selectWindow("Volume_Colocalization_Channel1");
			close();
			selectWindow("Volume_Colocalization_Channel2");
			close();
			selectWindow("Common_volumes");
			close();
			
		}//loop for each nucleus
		run("Close All");
	} //loop for each image
}//loop for each file
print("all finished");


//--------------------------------------
//Functions
//--------------------------------------

function getFociCount(objectsMap){
	//Empties any pre-existing results table
	run("Clear Results");
	
	//Activate objects’ map
	selectWindow(objectsMap);
	run("Z Project...", "projection=[Max Intensity]");

	//Get and store the number of objects
	getStatistics(area, mean, min, nObjects, std, histogram);
	selectWindow("MAX_"+objectsMap);
	close();
	
	setResult("Foci Count", 0, nObjects);
	setResult("strain", 0, strain);
	setResult("treatment", 0, treatment);
	setResult("point", 0, point);
	setResult("replicate", 0, replicate);
	setResult("nucleus", 0, n);
	setResult("time frame", 0, timepoint);
	updateResults();	
}

//Retrieve volumes object per object
function getValues(objectsMap, imageToQuantify){
	//Activate objects’ map
	selectWindow(objectsMap);
	run("Z Project...", "projection=[Max Intensity]");

	//Get and store the number of objects
	getStatistics(area, mean, min, nObjects, std, histogram);
	selectWindow("MAX_"+objectsMap);
	close();
	
	//Create an output array, properly dimensioned
	measures=newArray(nObjects);

	//For each object
	for(i=1; i<=nObjects; i++){
		//Activate the objects’ map 
		selectWindow(objectsMap);

		//Set the threshold to select the current object
		setThreshold(i, i);

		//Empty the ROI Manager
		roiManager("Reset");

		//Run analyze particles, adding the outlines to the ROI Manager
		run("Analyze Particles...", "add stack");

		//Create a variable to store the volume and initialise it to zero
		singleMeasure=0;

		//For each outline
		for(j=0; j<roiManager("Count"); j++){
			//Activate the image on which to measure
			selectWindow(imageToQuantify);

			//Select the ROI
			roiManager("Select", j);

			//Measure the volume
			getStatistics(area, mean, min, max, std, histogram);

			//Add the volume to the variable
			singleMeasure+=area*mean;
		}
		//End for each outline
	
		//Push the measure to the output array
		measures[i-1]=singleMeasure;
	
	//End for each object
	}
	
	//Return the output array
	return measures;
}

//Retrieve intensity per object
function getIntensity(objectsMap, imageToQuantify){
	//Activate objects’ map
	selectWindow(objectsMap);
	run("Z Project...", "projection=[Max Intensity]");

	//Get and store the number of objects
	getStatistics(area, mean, min, nObjects, std, histogram);
	selectWindow("MAX_"+objectsMap);
	close();
	
	//Create an output array, properly dimensioned
	measures=newArray(nObjects);

	//For each object
	for(i=1; i<=nObjects; i++){
		//Activate the objects’ map 
		selectWindow(objectsMap);

		//Set the threshold to select the current object
		setThreshold(i, i);

		//Empty the ROI Manager
		roiManager("Reset");

		//Run analyze particles, adding the outlines to the ROI Manager
		run("Analyze Particles...", "add stack");

		//Create a variable to store the volume and initialise it to zero
		singleMeasure=0;

		//For each outline
		for(j=0; j<roiManager("Count"); j++){
			//Activate the image on which to measure
			selectWindow(imageToQuantify);

			//Select the ROI
			roiManager("Select", j);

			//Measure the volume
			getStatistics(area, mean, min, max, std, histogram);

			//Add the volume to the variable
			singleMeasure+=area*mean;
		}
		//End for each outline
	
		//Push the measure to the output array
		measures[i-1]=singleMeasure;
	
	//End for each object
	}
	
	//Return the output array
	return measures;
}

//Generates two types of outputs: a results table and 2 co-localisation maps
function generateOutputs(objectsMeasures, commonMeasures, ObjectInt, objectsMap){
	//Empties any pre-existing results table
	run("Clear Results");

	//Duplicate the objects map
	selectWindow(objectsMap);
	run("Select None"); //Needed to remove any ROI from the image
	run("Duplicate...", "title=Coloc_Map duplicate");
	run("32-bit"); //Needed to accomodate decimal intensities
	
	for(i=0; i<objectsMeasures.length; i++){
		//Calculate the ratio
		ratio=commonMeasures[i]/objectsMeasures[i];
		ObjectMeanInt=ObjectInt[i]/objectsMeasures[i];
		
		//Fill the results table with data
		setResult("Label", nResults, "Object_"+(i+1));
		setResult("Full object", nResults-1, objectsMeasures[i]);
		setResult("Common part", nResults-1, commonMeasures[i]);
		setResult("Ratio", nResults-1, ratio);
		setResult("ObjectInt", nResults-1, ObjectInt[i]);
		setResult("ObjectMeanInt", nResults-1, ObjectMeanInt);
		setResult("strain", nResults-1, strain);
		setResult("treatment", nResults-1, treatment);
		setResult("point", nResults-1, point);
		setResult("replicate", nResults-1, replicate);
		setResult("nucleus", nResults-1, n);
		setResult("time frame", nResults-1, timepoint);
		updateResults();	
		//Replace each object's tag by the corresponding colocalisation ratio
		selectWindow("Coloc_Map");
		run("Macro...", "code=[if(v=="+(i+1)+") v="+ratio+"] stack");
	}
	resetMinAndMax();
}