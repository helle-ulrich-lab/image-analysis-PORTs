//--------------------------------------
//This macro performs foci counting per cell. It uses DAPI to create a nuclear mask and performs segmentation of foci in 2D
//author: Ronald Wong
//--------------------------------------

//Initialize
print("\\Clear");
setBatchMode(true);

//Create table to write data
run("New... ", "name=results-nuclei type=Table");
print("[results-nuclei]","\\Headings:\t \tIntDen\tRawIntDen\tseriesNumber\tFile");


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
		filename=substring(fileList[k],m+1,q);
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
				run("Z Project...", "projection=[Max Intensity]");
				run("Split Channels");
				close("C3-MAX_stack");
				selectWindow("C1-MAX_stack");
				rename("foci");
				selectWindow("C2-MAX_stack");
				rename("nuclei");
				//segment nuclei
				selectWindow("nuclei");
				run("FeatureJ Laplacian", "compute smoothing=5");
				setAutoThreshold("MaxEntropy");
				run("Convert to Mask");
				run("Set Measurements...", "area redirect=None decimal=3");
				run("Analyze Particles...", "size=2-10 add display exclude clear");
				close("nuclei");
				print("nuclei segmentation finished");
				numberNuclei=roiManager("count");
				print(""+numberNuclei+" nuclei found");
						
				if (numberNuclei>0) {
				
					//segment foci
					run("Clear Results");
					selectWindow("foci");
					run("Duplicate...","title=gaussian");
					run("Gaussian Blur...","sigma=2");
					imageCalculator("Subtract create", "foci","gaussian");
					rename("foci-segmented");
					setAutoThreshold("MaxEntropy dark");
					run("Convert to Mask");
					close("gaussian");
					print("foci segmentation finished");
					run("Set Measurements...", "centroid redirect=None decimal=3");
					run("Analyze Particles...", "size=0.2-0.8 display");
					numberFoci=nResults();
					print("foci quantification finished");  
					print(""+numberFoci+ " foci found");
					if (numberFoci>0) {
  
			   		//count foci- a blank image will be created to put each focus as a single pixel for counting
			   		selectWindow("foci");
					run("Duplicate...","title=foci-center");
					setBackgroundColor(0, 0, 0);
					run("Select All");
					run("Clear");
					for (i=0;i<numberFoci;i++) { 
						x=getResult("X", i);
						y=getResult("Y", i);
						toUnscaled(x,y);
						setPixel(x,y,1);
					}
						//count foci in each nucleus
						run("Properties...", "channels=1 slices=1 frames=1 unit=pixel pixel_width=1 pixel_height=1 voxel_depth=1");
						run("Clear Results");
						run("Set Measurements...", "integrated redirect=None decimal=3");
						for (i=0;i<numberNuclei;i++) {
							selectWindow("foci-center");
							roiManager("select",i);
							run("Measure");
							setResult("seriesNumber", i, seriesNumber);
							setResult("File", i, filename);
						}
					//copy data to table
					updateResults();
		        	String.copyResults();
		        	print("[results-nuclei]",String.paste);
		       		run("Clear Results");
		     	
		    		print("foci counting finished");
		
			
					} //if numberFoci>0 
	
				run("Close All");
				} //if number of nuclei >0 
			} //loop for series
	} //if file is lif
} //loop for all files
print("all finished");