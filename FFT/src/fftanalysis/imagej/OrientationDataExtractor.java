package fftanalysis.imagej;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;

public class OrientationDataExtractor implements PlugInFilter {
	
	protected ImagePlus image;
	
	private int width;
	private int height;
	
	public double value;
	public String name;

	@Override
	public void run(ImageProcessor ip) {
		width = ip.getWidth();
		height = ip.getHeight();
		
		if (!showDialog()) return;
		// Placeholder behavior: log pixel info. (Real processing lives in FIBA_Orientation.)
		IJ.log("OrientationDataExtractor: width=" + width + ", height=" + height + ", value=" + value + ", name=" + name);
		if (image != null) image.updateAndDraw();
	}
	
	private boolean showDialog() {
		GenericDialog gd = new GenericDialog("Process pixels");
		
		gd.addNumericField("value", 0.00, 2);
		gd.addStringField("name", "Mike");
		gd.showDialog();
		if (gd.wasCanceled()) return false;
		value = gd.getNextNumber();
		name = gd.getNextString();
		return true;
	}

	@Override
	public int setup(String arg0, ImagePlus arg1) {
		this.image = arg1;
		return DOES_ALL;
	}

}
