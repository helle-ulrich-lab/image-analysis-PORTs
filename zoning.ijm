//--------------------------------------
// This macro counts foci and determines the location of each focus
// with respect to three nuclear zones of equal volume. It requires a
// fluorescent channel that marks the nuclear volume and a target
// fluorescent protein.
// author: Ronald Wong and Mária Hanulová
//--------------------------------------

//Initialize
setBatchMode(true);
print("\\Clear");

//basic parameters to be set- size of nuclei and foci before scaling in pixels
//size limit (in pixel) for nuclei
min_size=500;
max_size=7500;
//size limit (in pixel) for foci
min_focus=5;
max_focus=100;
//pixels of bounding box
pad=20; 
//z extension for padding
ext=5;
//exapnsion factor for creating an expanded ellipsoid before dividing into three zones
scaleXY=4;
scaleZ=2*scaleXY;
scaleV=scaleXY*scaleXY*scaleZ;
//threshold values for foci segmentation
threshold_foci=1200;
minsize_foci=5;
//circularity limit (from 0-1, 1 being a perfect sphere)
circ_3d=0.2;

//create table to write data
run("New... ", "name=zones_table type=Table"); //zone dimensions desired and real
print("[zones_table]","\\Headings:object\twhole\t2thirds_desired\t2thirds_real\t2thirds_difference\t2hirds_difference_fraction\t1third_desired\t1third_real\t1third_difference\t1third_difference_fraction\tstrain\ttreatment\tpoint\treplicate\ttime frame");
run("New... ", "name=foci_table type=Table"); //foci in zones
print("[foci_table]","\\Headings:nucleusNo\tfocusNo\tfocusZone\tx\ty\tz\tstrain\ttreatment\tpoint\treplicate\ttime frame");
run("New... ", "name=Number_of_foci_combined type=Table"); //foci counting
print("[Number_of_foci_combined]","\\Headings: \tFoci count\tstrain\ttreatment\tpoint\treplicate\tnucleus\ttime frame");


//get directory with tif files
dir=getDirectory("Choose directory with the DeltaVision files");
print("processing directory ",dir);	
fileList=getFileList(dir);
nFiles=lengthOf(fileList);
print(nFiles+" files found");
for(k=0;k<nFiles;k++) {  //loop over all files in the directory
	if (matches(fileList[k],".*R3D_D3D.dv")) { //opne DeltaVision file
		print("Processing "+fileList[k]);
		run("Bio-Formats Importer", "open="+dir+fileList[k]+" color_mode=Colorized view=Hyperstack stack_order=XYCZT");
		//get name of file to write to the table
		m=indexOf(fileList[k],"_");
		o=indexOf(fileList[k],"_", m+1);
		p=indexOf(fileList[k],"_", o+1);
		q=indexOf(fileList[k],"_", p+1);
		r=indexOf(fileList[k],"_", q+1);
		strain=substring(fileList[k],m+1,o);
		treatment=substring(fileList[k],o+1,p);
		timepoint=substring(fileList[k],p+1,q);
		point=substring(fileList[k],q+1,r);
		point_num=parseFloat(point);
		replicate=point_num%3+1;

		//crop image
		makeRectangle(32, 32, 960, 960);
		run("Crop");
		//Get image name		
		name=getTitle();
		run("Split Channels");
		
		//nucleus segmentation
		selectWindow("C1-"+name);
		run("3D Fast Filters","filter=Median radius_x_pix=2.0 radius_y_pix=2.0 radius_z_pix=2.0 Nb_cpus=4");
		run("FeatureJ Laplacian", "compute smoothing=6");
		setAutoThreshold("Otsu no-reset stack");
		setOption("BlackBackground", true);
		run("Convert to Mask", "method=Otsu background=Light black");
		run("Fill Holes", "stack");
		
		//test circularity in 3D and take out nuclei that are not suitable for the analysis
		run("3D OC Options", "volume surface nb_of_obj._voxels nb_of_surf._voxels dots_size=5 font_size=10 redirect_to=none");
		run("3D Objects Counter", "threshold=128 slice=1 min.="+min_size+" max.="+max_size+" exclude_objects_on_edges objects statistics");
		rename("mask_3d_initial");
		getDimensions(width,height,channels,slices,frames);
		for (i=0;i<nResults;i++) {
			surface=getResult("Surface (micron^2)",i);
			volume=getResult("Volume (micron^3)",i);
			print(i,36*PI*pow(volume,2)/pow(surface,3));
			if (36*PI*pow(volume,2)/pow(surface,3)<circ_3d) {
				print("nucleus",i,"removed");
				for (j=1;j<=slices;j++) {
					Stack.setSlice(j);
					changeValues(i+1, i+1, 0); 
				}	
			}
		}
		//If nuclei are found, an array of coordinates with the bounding box for the nuclei will be created
		selectWindow("mask_3d_initial");
		Stack.getStatistics(voxelCount, mean, min, max);
		if (max>0) {
			setThreshold(1,nResults);
			run("Convert to Mask", "method=RenyiEntropy background=Dark");
			run("3D OC Options", "volume nb_of_obj._voxels bounding_box dots_size=5 font_size=10 redirect_to=none");
			run("3D Objects Counter", "threshold=128 slice=1 min.="+min_size+" max.="+max_size+" objects statistics");
			rename("mask_3d");
			close("C1-"+name);
			close("3D_Maximum");
			print(nResults+" nuclei found");
			//arrays for bounding box - x,y,width,heigth
			bx=newArray(nResults);
			by=newArray(nResults);
			bw=newArray(nResults);
			bh=newArray(nResults);
			for (i=0;i<nResults;i++) {
				bx[i]=getResult("BX",i)-pad/2;
				by[i]=getResult("BY",i)-pad/2;
				bw[i]=getResult("B-width",i)+pad;
				bh[i]=getResult("B-height",i)+pad; 
			}
			Array.print(bx);
			Array.print(by);
			Array.print(bw);
			Array.print(bh);
						
			//segment foci
			selectWindow("C2-"+name);
			setThreshold(threshold_foci, 70000);
			run("Convert to Mask", "method=Default background=Dark black");
			run("3D OC Options", "  dots_size=5 font_size=10 redirect_to=none");
			run("3D Objects Counter", "threshold=0 slice=1 min.="+minsize_foci+" max.=100000 exclude_objects_on_edges objects");
			setThreshold(1, 70000);
			run("Convert to Mask", "method=Default background=Dark black");
			rename("foci");
			close("C2-"+name);
			
			//Each nucleus is analyzed individually
			nObjects=bx.length;
			print(nObjects);
			for (i=1;i<=nObjects;i++) {
				print("object",i);
				
				//count foci
				//crop out nucleus from the segmented foci image
				selectWindow("foci");
				run("Duplicate...","title=focus duplicate");
				run("Specify...", "width="+bw[i-1]+" height="+bh[i-1]+" x="+bx[i-1]+" y="+by[i-1]);
				run("Crop");
				//create object map to measure how many objects in this box
				selectWindow("focus");
				run("3D Objects Counter", "threshold=0 slice=1 min.="+minsize_foci+" max.=100000 exclude_objects_on_edges objects");
				rename("ObjectMap");
				getFociCount("ObjectMap");
				String.copyResults();
				print("[Number_of_foci_combined]",String.paste);
				close("ObjectMap");
				run("Clear Results");

				//zone analysis
				//check if nucleus contains foci
				selectWindow("focus");
				Stack.getStatistics(area,mean,min,max);
				if (max==255) { //foci present, continue with analysis
					//The foci image is expanded according to the expansion factor to create a more accurate zoning analysis
					scaleAndPad("focus",bw[i-1],bh[i-1],scaleXY,scaleZ,pad);
					//crop out nucleus
					selectWindow("mask_3d");
					run("Duplicate...","title=nucleus duplicate");
					setThreshold(i,i);
					run("Convert to Mask", "method=RenyiEntropy background=Dark");
					run("Specify...", "width="+bw[i-1]+" height="+bh[i-1]+" x="+bx[i-1]+" y="+by[i-1]);
					run("Crop");									
					//The nucleus image is expanded according to the expansion factor to create a more accurate zoning analysis
					scaleAndPad("nucleus",bw[i-1],bh[i-1],scaleXY,scaleZ,pad);
			
					//fit ellipse to the nucleus
					selectWindow("nucleus");
					run("3D Ellipsoid Fitting", " ");
					selectWindow("Ellipsoids");
					resetMinAndMax;
					run("8-bit");
					close("nucleus");	
					run("3D OC Options", "nb_of_obj._voxels dots_size=5 font_size=10 redirect_to=none");
					//generate zones
					zones("Ellipsoids",i); 
					//Measure in which zones foci are situated
					zones_measure("all_zones","focus",i); 
					run("Merge Channels...", "c4=all_zones c6=focus create");
					close();
				} //if foci present
				else {
					close("focus");
					close("nucleus");
					print("[foci_table]",""+i+"\tNA\tNA\tNA\tNA\tNA\t"+strain+"\t"+treatment+"\t"+point+"\t"+replicate+"\t"+timepoint); //if no foci
				}
			} 	run("Close All");
		} run("Close All");
	} //loop for all Deltavision files
}//loop for all files in directory
print("all finished");

///////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////functions////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
function zones(all,n) { //parameter: nucleus to zone and nucleus number
	selectWindow(all);
	run("3D Objects Counter", "threshold=128 slice=1 min.=500 max.=2000000 statistics");
	vol_current=getResult("Nb of obj. voxels",0);
	vol_next=vol_current;
	vol_whole=vol_current;
	vol_2thirds=2*vol_current/3;
	vol_1third=vol_current/3;
	run("Duplicate...","duplicate title=current");
	run("Duplicate...","duplicate title=next");
	while (vol_next>vol_2thirds) {
		close("current");
		vol_current=vol_next;
		selectWindow("next");
		rename("current");
		run("3D Fast Filters","filter=Minimum radius_x_pix=1.0 radius_y_pix=1.0 radius_z_pix=1.0 Nb_cpus=4");//minimum = erosion in 3D; creates new window, changes lut
		rename("next");
		run("3D Objects Counter", "threshold=128 slice=1 min.="+scaleV*min_size/3+" max.="+scaleV*max_size+" statistics");
		vol_next=getResult("Nb of obj. voxels",0); //volume in voxels
		print(vol_current,vol_next);
	}
	print("2/3");
	if ((vol_2thirds-vol_next)<(vol_current-vol_2thirds)) {
		selectWindow("next");
		run("Duplicate...","duplicate title=2thirds");
		vol_2thirds_real=vol_next;
		}
	else {
		close("next");
		selectWindow("current");
		run("Duplicate...","duplicate title=2thirds");
		run("Duplicate...","duplicate title=next");
		vol_2thirds_real=vol_current;
		vol_next=vol_current;
	}
	while (vol_next>vol_1third) {
		close("current");
		vol_current=vol_next;
		selectWindow("next");
		rename("current");
		run("3D Fast Filters","filter=Minimum radius_x_pix=1.0 radius_y_pix=1.0 radius_z_pix=1.0 Nb_cpus=4");//creates new window, changes lut
		rename("next");
		run("3D Objects Counter", "threshold=128 slice=1 min.="+scaleV*min_size/3+" max.="+scaleV*max_size+" statistics");
		vol_next=getResult("Nb of obj. voxels",0); //volume in voxels
		print(vol_current,vol_next);
	}
	print("1/3");
	if ((vol_1third-vol_next)<(vol_current-vol_1third)) {
		selectWindow("next");
		rename("zone3");
		close("current");
		vol_1third_real=vol_next;
		}
	else {
		close("next");
		selectWindow("current");
		rename("zone3");
		vol_1third_real=vol_current;
	}
	print("desired volumes",vol_whole,vol_2thirds,vol_1third);
	print("actual volumes",vol_whole,vol_2thirds_real,vol_1third_real);

	print("[zones_table]",""+n+"\t"+vol_whole+"\t"+vol_2thirds+"\t"+vol_2thirds_real+"\t"+
	vol_2thirds-vol_2thirds_real+"\t"+(vol_2thirds-vol_2thirds_real)/vol_2thirds+"\t"+vol_1third+"\t"+vol_1third_real+"\t"+vol_1third-vol_1third_real+"\t"+(vol_1third-vol_1third_real)/vol_1third+"\t"+strain+"\t"+treatment+"\t"+point+"\t"+replicate+"\t"+timepoint);
	
	imageCalculator("Subtract create stack", all,"2thirds");
	rename("zone1");
	imageCalculator("Subtract create stack", "2thirds","zone3");
	rename("zone2");
	run("Invert LUT");
	
	//visualize
	selectWindow("zone1");
	run("Duplicate...","title=z1 duplicate");
	run("Subtract...", "value=254 stack");
	selectWindow("zone2");
	run("Duplicate...","title=z2 duplicate");
	run("Subtract...", "value=253 stack");
	selectWindow("zone3");
	run("Invert LUT");
	run("Duplicate...","title=z3 duplicate");
	run("Subtract...", "value=252 stack");
	imageCalculator("Add stack", "z1","z2");
	imageCalculator("Add stack", "z1","z3");
	rename("all_zones");
	close("z3");
	close("z2");
	close("2thirds");
	close("zone1");
	close("zone2");
	close("zone3");
	close("Ellipsoids");
} //function zones

function zones_measure(z,f,n) { //zones and foci image
	//after scaling the image, scale is gone, everything is in pixels
	//measure intensity of the central focus pixel in the zones image
	selectWindow(f);
	run("3D OC Options", "centroid dots_size=5 font_size=10 redirect_to=none");
	run("3D Objects Counter", "threshold=128 slice=1 min.="+scaleV*min_focus+" max.="+scaleV*max_focus+" statistics");
	for (i=0;i<nResults;i++) {
		fx=getResult("X",0);
		fy=getResult("Y",0);
		fz=getResult("Z",0);
		selectWindow(z);
		Stack.setSlice(fz+0.5); //whole part is taken
		run("Specify...", "width=1 height=1 x="+fx+" y="+fy); //rounding here
		getStatistics(area, mean);
		print("[foci_table]",""+n+"\t"+i+1+"\t"+mean+"\t"+fx+"\t"+fy+"\t"+fz+"\t"+strain+"\t"+treatment+"\t"+point+"\t"+replicate+"\t"+timepoint);
	}
}

function scaleAndPad(image,w,h,sxy,sz,p) {// scaling and padding of image
	selectWindow(image);
	//extend in z
	b=bitDepth(); 
	newImage("ext1", ""+b+"-bit black", w, h, ext);
	newImage("ext2", ""+b+"-bit black", w, h, ext);
	run("Concatenate...", "open image1=ext1 image2="+image+" image3=ext2 image4=[-- None --]");
	rename(image);
	//scale
	selectWindow(image);
	run("Scale...", "x="+sxy+" y="+sxy+" z="+sz+" interpolation=None process");
	close(image);	
	selectWindow(image+"-1");
	rename(image);
}

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
	setResult("nucleus", 0, i);
	setResult("time frame", 0, timepoint);
	updateResults();	
}