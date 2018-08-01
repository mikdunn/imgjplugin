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
		// TODO Auto-generated method stub
		width = ip.getWidth();
		height = ip.getHeight();
		
		if (showDialog()) {
			process(ip);
			image.updateAndDraw();
		}
		
		
	}
	
	private boolean showDialog() {
		GenericDialog gd = new GenericDialog("Process pixels");
		
		gd.addNumericField("value", 0.00, 2);
		gd.addStringField("name", "Mike");
		
		
	}

	@Override
	public int setup(String arg0, ImagePlus arg1) {
		// TODO Auto-generated method stub
		return 0;
	}

}
