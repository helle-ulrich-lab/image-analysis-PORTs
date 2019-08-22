//--------------------------------------
//This macro performs foci counting per cell. It uses DAPI to create a nuclear mask and performs segmentation of foci in 3D
//author: Ronald Wong
//--------------------------------------

//Initialize
print("\\Clear");
setBatchMode(true);

//Create table to write data
run("New... ", "name=results-nuclei type=Table");
print("[results-nuclei]","\\Headings:\t \tIntDen\tRawIntDen\tSeriesNumber\tFileNumber\tFileName");

//get directory with image files
dir=getDirectory("Choose directory with the image files");
print("processing directory ",dir);	
fileList=getFileList(dir);
nFiles=lengthOf(fileList);
print(nFiles+" files found");
for(k=0;k<nFiles;k++) {  //loop over all files in the directory


	if (matches(fileList[k],".*lif")) { //look for Lif files for analysis
		//get name of file to write to the table
		m=indexOf(fileList[k],"_");
		o=indexOf(fileList[k],"_", m+1);
		p=indexOf(fileList[k],"_", o+1);
		q=indexOf(fileList[k],".");
		filename=substring(fileList[k],0,q);
		print("Processing "+fileList[k]);
		//open image file
		run("Bio-Formats Macro Extensions"); 
		Ext.setId(dir+fileList[k]);
		Ext.getSeriesCount(seriesCount);
		print("there are "+seriesCount+" series in your lif file");

		for (seriesNumber=1;seriesNumber<seriesCount+1;seriesNumber++){ //analyze all series in a Lif file
			roiManager ("reset");
			run("Bio-Formats Importer", "open="+dir+fileList[k]+" color_mode=Colorized view=Hyperstack stack_order=XYCZT series_"+seriesNumber);
 			print("opened series "+seriesNumber);
			rename("stack");
			//normalize name
			run("Split Channels");
			close("C3-stack");
			selectWindow("C1-stack");
			rename("foci");
			selectWindow("C2-stack");
			rename("nuclei");
			//segment nuclei
			selectWindow("nuclei");
			run("Z Project...", "projection=[Max Intensity]");
			run("FeatureJ Laplacian", "compute smoothing=8");
			setAutoThreshold("MaxEntropy");
			run("Convert to Mask");
			run("Set Measurements...", "area redirect=None decimal=3");
			run("Analyze Particles...", "size=2-Infinity add display exclude clear");
			close("nuclei");
			print("nuclei segmentation finished");
			numberNuclei=roiManager("count");
			print(""+numberNuclei+" nuclei found");
			
			if (numberNuclei>0) {
				//foci segmentation with TopHat filter and auto-thresolding
				selectWindow("foci");
				run("3D Fast Filters","filter=TopHat radius_x_pix=2.0 radius_y_pix=2.0 radius_z_pix=2.0 Nb_cpus=4");
				setAutoThreshold("MaxEntropy dark stack");
				run("Convert to Mask", "method=MaxEntropy background=Dark black");
				
				run("3D OC Options", "volume nb_of_obj._voxels centre_of_mass dots_size=5 font_size=10 show_numbers white_numbers redirect_to=[foci]");
				run("3D Objects Counter", "threshold=128 slice=13 min.=5 max.=37158912 exclude_objects_on_edges objects statistics summary");
				numberFoci=nResults();
				//count foci- a blank image stack will be created to put each focus as a single pixel for counting in 3D
				if (numberFoci>0){
					selectWindow("foci");
					getDimensions(width, height, channels, slices, frames);
					newImage("foci-center", "8-bit grayscale-mode", width , height, slices);
					run("Select All");
					run("Clear","stack");
					for (i=0;i<numberFoci;i++) { 
						x=getResult("XM", i);
						y=getResult("YM", i);
						z=getResult("ZM",i);
						toUnscaled(x,y);
						setSlice(z);
						setPixel(x,y,1);
						}
					//count foci in each nucleus
					run("Z Project...", "projection=[Sum Slices]");
					run("Clear Results");
					run("Set Measurements...", "integrated redirect=None decimal=3");
					for (i=0;i<numberNuclei;i++) {
						selectWindow("SUM_foci-center");
						roiManager("select",i);
						run("Measure");
						setResult("SeriesNumber", i, seriesNumber);
						setResult("FileNumber", i, k);
						setResult("FileName", i, filename);
						}
					//copy data to table
					updateResults();
	        		String.copyResults();
	        		print("[results-nuclei]",String.paste);
	       			run("Clear Results");
	     			print("foci counting finished");
   
				} // if numberFoci>0
			}run("Close All"); //if numberNuclei>0 
		} //loop for all series
	} //if file is lif
} //loop for all files
print("all finished");