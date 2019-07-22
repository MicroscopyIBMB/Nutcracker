// Name: Dm_Nutcracker.ijm Nov.2018
// Authors: Elena Rebollo & Jaume Boix, Molecular Imaging Platform IBMB, Barcelona
// Fiji version: Fiji lifeline 22 Dec 2015
/*This macro removes the viteline membrane autofluorescence signal from 4D Drosophila embryo movies 
  acquired in confocal or multiphoton microscopes*/

/*Short Instructions:
Open the 4D Hyperstack and run the macro
Dialog box: 
o	Substract Background Radius: Use 20-50 for medium and high autofluorescence intensity in the vitelline membrane. In case of low intensity set a high value (1000). 
o	CLAHE Block Size: Use 8-10 for weak autofluorescence halo. Use 30-40 for medium and high autofluorescence halo.
o	Gaussian: Use 3 for pixel size of 1 micrometer. 
o	LoG Radius: Use 3 as standard. This value could be increased when there are no bright objects close to vitelline membrane 
o   Fitting Line Length: Length of the perpendicular line of the fitting algorithm. 1/3 of its length is placed outside the ellipse.
o	Adjust ROI: number of pixels to shrink the ROI just underneath the vitelline membrane
o	Embryo Close to edge: Activate this option if the embryo is very close or touches the image borders. It avoids some artifacts during segmentation.
o	Remove floating particles: Activate when floating particles appear outside the embryo

Next dialog: select the image plane from which the sections should be removed (macro will not work for planes made parallel to the embryo viteline membrane.

The macro will run until the end. It will take a long time (30-60’) for big files like 1 GB.
 */



//CREATE DIALOG TO ADJUST PARAMETERS
Dialog.create("Choose parameters");
Dialog.addNumber("Substract Background radius:", 20);
Dialog.addNumber("CLAHE block size:", 30);
Dialog.addNumber("Gaussian Blur Radius:",3);
Dialog.addNumber("LoG Radius:",3);
Dialog.addNumber("Fitting Line Length (micrometer):",50);
Dialog.addNumber("Adjust ROI (px):", -6);
Dialog.addCheckbox("Embryo close to edge?", false);
Dialog.addCheckbox("Remove floating particles?", false);
Dialog.show;
Radius=Dialog.getNumber();
Block=Dialog.getNumber();
Gaussian=Dialog.getNumber();
LoGRadius=Dialog.getNumber();
Line=Dialog.getNumber();
FitROI=Dialog.getNumber();
Canvas=Dialog.getCheckbox();
particles=Dialog.getCheckbox();

if (particles==true) {
	Dialog.create("Remove particles parameters:");
	Dialog.addSlider("Minimum circularity", 0.5, 1.0, 0.8);
	Dialog.addNumber("Maximum particle size (px):",100);
	Dialog.show;
	circ=Dialog.getNumber();
	psize=Dialog.getNumber();
}


//GET FILE NAME & DIMENSIONS
RawName = getTitle;
Name = File.nameWithoutExtension;
getDimensions(width, height, channels, slices, frames);
getPixelSize(unit,pixelW,pixelH);
pixelsLine=round(Line/pixelW);

//SELECT FIRST PLANE TO DELETE EXTRA SLICES
waitForUser("Select plane to delete extra slices");
Stack.getPosition(channel, slice, frame);

//CONVERT HYPERSTACK TO STACK
run("Make Substack...", "slices="+slice+"-"+slices+" frames=1-"+frames);
getDimensions(width2, height2, channels2, slices2, frames2);
run("Hyperstack to Stack");
rename("Stack");
selectWindow(RawName);
run("Close");

//CREATE IMAGE TEMPLATE
if (Canvas==true) {
selectWindow("Stack");
width3 = width2 + 24;
height3 = height2 + 24;
run("Canvas Size...", "width="+d2s(width3, 0)+" height="+d2s(height3, 0)+" position=Center zero");
newImage("temp", "8-bit black", width3, height3, 1);
}
else {
newImage("temp", "8-bit black", width2, height2, 1);
}


//REMOVE OUTLINE AUTOFLUORESCENCE SIGNAL IN ALL FRAMES
selectWindow("Stack");
//Create For loop and insert workflow to be applied to ech image plane
setBatchMode(true);
for(i=1; i<=nSlices; i++) {
	setSlice(i);
	// Copy frame i to mask image
	run("Select All");
	run("Copy", "slice");
	selectWindow("temp");
	run("Paste");
	roiManager("reset");
	run("Subtract Background...", "rolling="+d2s(Radius,0));
	run("Enhance Local Contrast (CLAHE)", "blocksize="+d2s(Block,0)+" histogram=8 maximum=8 mask=*None*");
	run("Gaussian Blur...", "sigma="+d2s(Gaussian,0));
	run("FeatureJ Laplacian", "compute smoothing="+d2s(LoGRadius,0));
	run("8-bit");
	run("Invert");
	setAutoThreshold("Otsu dark");
	setOption("BlackBackground", false);
	run("Convert to Mask");

	//Remove floating particles outside the embryo (option)
	if (particles==true) {
	run("Analyze Particles...", "size=0-Infinity circularity="+d2s(circ,2)+"-1.00 clear add");
	run("Analyze Particles...", "size=0-"+d2s(psize,0)+" circularity=0.00-1.00 add");
	setForegroundColor(255, 255, 255);
	roiManager("deselect");
	roiManager("fill");
	roiManager("reset");
	}
	run("Skeletonize");
	run("Options...", "iterations=1 count=1");
	run("Dilate");
	rename("skeleton");
	run("Points from Mask");
	run("Convex Hull");
	run("Fit Ellipse");

// SHAPE FITTING ALGORITHM
	// Get reference ellipse coordinates into arrays and create more arrays to store future coordinates
	roiManager("add");
	Roi.getCoordinates(X,Y);
	Xref=Array.concat(X,X[0]);
	Yref=Array.concat(Y,Y[0]);
	Xfit=newArray();
	Yfit=newArray();
	Xline=newArray(pixelsLine);
	Yline=newArray(pixelsLine);

	// Finding coordinates along perpendicular imaginary line
	for (k=0;k<(Xref.length-1);k++) {
    	dx=Xref[k+1]-Xref[k];
    	dy=Yref[k+1]-Yref[k];
    	a=atan2(-dy,dx);
    	z=a-(PI/2);
		b=searchMax(Xline,Yline,Xref[k],Yref[k],z,pixelsLine);
		if (b<pixelsLine) {
		Xfit=Array.concat(Xfit,Xline[b-1]);
    	Yfit=Array.concat(Yfit,Yline[b-1]);
		} else {
		Xfit=Array.concat(Xfit,Xref[k]);
		Yfit=Array.concat(Yfit,Yref[k]);	
		}
	}

	// Extract fitted coordinates into polygonal ROI
	makeSelection("polygon",Xfit,Yfit);
	roiManager("reset");
	roiManager("add");

	// AUTOFLUORESCENT HALO REMOVAL ON ORIGINAL IMAGE
	roiManager("Select", 0);
	run("Enlarge...", "enlarge="+d2s(FitROI, 0)+" pixel");
	roiManager("update");
	selectWindow("Stack");
	roiManager("Select", 0);
	setBackgroundColor(0, 0, 0);
	run("Clear Outside", "slice");
	roiManager("reset");
	selectWindow("skeleton");
	run("Close");
	selectWindow("Stack");
}

setBatchMode(false);

//CONVERT PROCESSED STACK INTO HYPERSTACK
selectWindow("temp");
run("Close");
run("Stack to Hyperstack...", "order=xyczt(default) channels=1 slices="+slices2+" frames=frames display=Color");
if (Canvas==true) {
run("Canvas Size...", "width="+d2s(width2, 0)+" height="+d2s(height2, 0)+" position=Center zero");
}
rename(Name+"_processed.tif");

// USER-DEFINED FUNCTION TO SEARCH FOR A VALUE ALONG COORDINATES ON AN IMAGINARY LINE
function searchMax(Xline,Yline, Xref,Yref,z,pixelsLine) {
   b=0;
   Int=0;
   do {
   	Xline[b]=round(Xref+(b-pixelsLine/3)*cos(z));       
    Yline[b]=round(Yref+(b-pixelsLine/3)*sin(-z));
    selectWindow("skeleton");
    Int=getPixel(Xline[b],Yline[b]);
    b++;
   } while (Int<255 && b<pixelsLine)
   return b,Xline,Yline
}
    
          
        