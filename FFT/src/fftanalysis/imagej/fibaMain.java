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
		int switchCase = 2; //default
		int displayType = 1; //default
		if(args.length == 1){
			switchCase = Integer.parseInt(args[0]);
		} else if(args.length == 2){
			switchCase = Integer.parseInt(args[0]);
			displayType = Integer.parseInt(args[1]);
		} else{
			System.out.println("Usage: java fiba " + "CaseNumber DisplayType");
			System.out.println("CaseNumber from 0 to 13 (Inclusive)");
			System.out.println("DisplayType from 0 to 2 (Inclusive)");
			System.out.println("Running Case " + switchCase + " by default.");
			System.out.println("Running DisplayType " + displayType + " by default");
		}//end else statement
		//CreateArray of test data
		int rows = 41;
		int cols = 41;
		//get a test surface in the spatial domain
		double[][] spatialData = getSpatialData(switchCase,rows,cols);

	}

}
