package fftanalysis.imagej;

import java.lang.Math.*;
import java.lang.Object.*;
import java.awt.*;
import java.awt.image.*;
import java.awt.image.BufferedImage;
import java.io.*;

public class fiba3D {
	
	
	public static void main(String[] args) {
		int numberRows = 59;
		int numberCols = 59;
		double[][] data = new double[numberRows][numberCols];
		int blockSize = 2;
		for(int row = 0; row < numberRows; row++){
			for(int col = 0; col < numberCols;col++){
				int xSquare = col*col;
				int ySquare = row*row;
				data[row][col] = xSquare + ySquare;
			} // end column loop
		} // end row loop
		fiba3D x = new fiba3D();
		

	}

}
