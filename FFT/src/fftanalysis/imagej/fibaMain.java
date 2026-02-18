package fftanalysis.imagej;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;

import java.lang.Math.*;
import java.lang.Object.*;
import java.awt.*;
import java.awt.image.*;
import java.awt.image.BufferedImage;
import java.io.*;

public class fibaMain {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// Minimal launcher for ImageJ when running from an IDE.
		// (The actual plugin entry point is fftanalysis.imagej.FIBA_Orientation.)
		new ImageJ();
		IJ.log("ImageJ started. Open an image and run Plugins > FIBA_Orientation.");

	}

}
