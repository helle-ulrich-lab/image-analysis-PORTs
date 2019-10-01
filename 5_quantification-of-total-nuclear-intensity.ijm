// This script performs quantification of total and mean nuclear
// fluorescent intensities. It requires a fluorescent channel for
// generating a nuclear mask (DAPI) and a target fluorescent protein.
// author: Ronald Wong

//Initialize
setBatchMode(true);
print("\\Clear");
run("Colors...", "foreground=white background=black selection=magenta");
roiManager("reset");

//Create tables to write data into
run("New... ", "name=Total_nuclear_Intensity_of_Channel1_combined type=Table");
print("[Total_nuclear_Intensity_of_Channel1_combined]","\\Headings: \tArea \tIntDen \tRawIntDen \tstrain\ttreatment\treplicate\tnucleus\ttime frame");

//get directory with tif files
dir=getDirectory("Choose directory with the DeltaVision files");
print("processing directory ",dir);	
fileList=getFileList(dir);
nFiles=lengthOf(fileList);
print(nFiles+" DeltaVision files found");
for(k=0;k<nFiles;k++) {  //loop over all files in the directory
	if (matches(fileList[k],".*R3D_D3D.dv")) {
		print("Processing "+fileList[k]);
		run("Bio-Formats Importer", "open="+dir+fileList[k]+" color_mode=Colorized view=Hyperstack stack_order=XYCZT");
		//Get file name to write into table
		m=indexOf(fileList[k],"_");
		o=indexOf(fileList[k],"_", m+1);
		p=indexOf(fileList[k],"_", o+1);
		q=indexOf(fileList[k],"_", p+1);
		r=indexOf(fileList[k],"_", q+1);
		s=indexOf(fileList[k],"_", r+1);
		strain=substring(fileList[k],m+1,o);
		treatment=substring(fileList[k],o+1,p);
		timepoint=substring(fileList[k],p+1,q);
		replicate=substring(fileList[k],q+1,r);
				
		//crop image
		makeRectangle(32, 32, 960, 960);
		run("Crop");
		ori=getTitle();
		selectWindow(ori);
		//Normalize names		
		run("Split Channels");
		selectWindow("C1-"+ori);
		rename("Channel1");
		selectWindow("C2-"+ori);
		rename("Channel2");
		selectWindow("C3-"+ori);
		rename("Channel3");
	
		//nuclear segmentation
		roiManager("reset");	
		selectWindow("Channel3");
		run("Z Project...", "projection=[Max Intensity]");
		run("Median...", "radius=5");
		setAutoThreshold("Triangle dark");
		setOption("BlackBackground", true);
		run("Convert to Mask");
		run("Dilate");
		run("Dilate");
		run("Analyze Particles...", "size=3-30 clear display exclude add");
		if (roiManager("count")>0) {
			//nuclear signal quantification
			run("Clear Results");
			selectWindow("Channel1");
			run("Z Project...", "projection=[Average Intensity]");
			run("Set Measurements...", "area integrated redirect=None decimal=3");
			for(i=0; i<roiManager("count"); i++){
				selectWindow("AVG_Channel1");
				roiManager("select", i);
				run("Measure");
				setResult("strain", nResults-1, strain);
				setResult("treatment", nResults-1, treatment);
				setResult("replicate", nResults-1, replicate);
				setResult("nucleus", nResults-1, i);
				setResult("time frame", nResults-1, timepoint);
				updateResults();	
			}
		String.copyResults();
		print("[Total_nuclear_Intensity_of_Channel1_combined]",String.paste);
		
		run("Clear Results");
		run("Close All");
		} else //loop through all nuclei if any found
		run("Close All");
	}//loop through all R3D_D3D.dv files
} //loop through all files
print("all finished");