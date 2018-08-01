package fftanalysis.imagej;

import java.lang.Math.*;
import java.lang.Object.*;
import java.awt.*;
import java.awt.image.*;
import java.awt.image.BufferedImage;
import java.io.*;


// This program contains the FFT and inverse FFT code (xform2D, inverseXform2D)
// Also includes shift origin, making the the coordinates centered at the origin
// Tukey window (window functions) and Histogram normalization will be on another program

public class fiba {
	static void xform2D(double[][] inputData, double[][] realOut, double[][] imagOut, double[][] amplitudeOut){
		int height = inputData.length;
		int width = inputData[0].length;
		System.out.println("height = "+height);
		System.out.println("width = "+width);
		if (width != height){
			System.exit(0);
		}
		else {
			for(int yWave = 0; yWave < height; yWave++){
				for(int xWave = 0; xWave < width; xWave++){
					for(int ySpace = 0; ySpace < height; ySpace++){
						for (int xSpace = 0; xSpace < width; xSpace++){
							realOut[yWave][xWave] += (inputData[ySpace][xSpace]*Math.cos(2*Math.PI*((1.0*xWave*xSpace
									/width) +(1.0*yWave*ySpace/height))))/Math.sqrt(width*height);
							imagOut[yWave][xWave] -= (inputData[ySpace][xSpace]*Math.sin(2*Math.PI*((1.0*xWave*xSpace
									/width) +(1.0*yWave*ySpace/height))))/Math.sqrt(width*height);
						} //end xSpace loop
					} //end ySpace loop
				} // end xWave loop
			}  // end yWave loop
		} // end xform2D
	}
	
	static void inverseXform2D(double[][] real, double[][] imag, double[][]dataOut){
		int height = real.length;
		int width = real[0].length;
		System.out.println("height = "+ height);
		System.out.println("width = "+ width);
		if (width != height){
			System.exit(0);
		}
		else {
			for(int ySpace = 0; ySpace < height; ySpace++){
				for(int xSpace = 0; xSpace < width; xSpace++){
					for(int yWave = 0; yWave < height; yWave++){
						for (int xWave = 0; xWave < width; xWave++){
							dataOut[ySpace][xSpace] += (real[yWave][xWave]*Math.cos(2*Math.PI*((1.0*xSpace*xWave
									/width) +(1.0*yWave*ySpace/height))))/Math.sqrt(width*height)-
							(imag[yWave][xWave]*Math.sin(2*Math.PI*((1.0*xSpace*xWave
									/width) + (1.0*ySpace*yWave/height))))/Math.sqrt(width*height);
						} //end xWave loop
					} //end yWave loop
				} // end xSpace loop
			}  // end ySpace loop
		} // end inverseXform2D
	}
	
	static double[][] shiftOrigin(double[][] data){
		int numberOfRows = data.length;
	    int numberOfCols = data[0].length;
	    int newRows;
	    int newCols;
	    double[][] output = new double[numberOfRows][numberOfCols];
	    if(numberOfRows%2 != 0){
	    	newRows = numberOfRows + (numberOfRows + 1)/2;
	    } else {
	    	newRows = numberOfRows + numberOfRows/2;
	    }
	    
	    if(numberOfCols%2 != 0){
	    	newCols = numberOfCols + (numberOfCols +1)/2;
	    } else {
	    	newCols = numberOfCols + numberOfCols/2;
	    }
	    
	    double[][] temp = new double[newRows][newCols];
	    
	    for (int row = 0;row < numberOfRows;row++){
	    	for (int col = 0;col < numberOfCols;col++){
	    		temp[row][col] = data[row][col];
	    	}//col loop (evens
	    }//row loop
	    
	    if(numberOfCols%2 != 0){//odd shift
	    	for(int row = 0; row < numberOfRows;row++){
	    		for (int col = 0; col < (numberOfCols+1)/2;col++){
	    			temp[row][col + numberOfCols] = temp[row][col];
	    		}//col loop (odds)
	    	}// row loop (odds)
	    	for(int row = 0;row < numberOfRows;row++){
	            for(int col = 0;col < numberOfCols;col++){
	              temp[row][col] = temp[row][col+(numberOfCols + 1)/2];
	            }//col loop (odd shift back)
	          }//row loop (odd shift back)
	    }
	    else{//even shift
	    	for(int row = 0;row < numberOfRows;row++){
	            for(int col = 0; col < numberOfCols/2;col++){
	              temp[row][col + numberOfCols] = temp[row][col];
	            }//col loop
	          }//row loop
	          for(int row = 0;row < numberOfRows;row++){
	            for(int col = 0; col < numberOfCols;col++){
	              temp[row][col] = temp[row][col + numberOfCols/2];
	            }//col loop
	          }//row loop
	        }//end else
	        //Now do the vertical shift
	        if(numberOfRows%2 != 0){//shift for odd
	         for(int col = 0;col < numberOfCols;col++){
	            for(int row = 0; row < (numberOfRows+1)/2;row++){
	              temp[row + numberOfRows][col] = temp[row][col];
	            }//row loop
	          }//col loop
	         for(int col = 0;col < numberOfCols;col++){
	            for(int row = 0; row < numberOfRows;row++){
	              temp[row][col] = temp[row+(numberOfRows + 1)/2][col];
	            }//row loop
	          }//col loop
	         
	        }else{//shift for even
	          //Slide topmost (numberOfRows/2) rows down
	          // by numberOfRows rows
	          for(int col = 0;col < numberOfCols;col++){
	            for(int row = 0; row < numberOfRows/2;row++){
	              temp[row + numberOfRows][col] = temp[row][col];
	            }//row loop
	          }//col loop
	          for(int col = 0;col < numberOfCols;col++){
	            for(int row = 0;row < numberOfRows;row++){
	              temp[row][col] = temp[row + numberOfRows/2][col];
	            }//row loop
	          }//col loop
	        }//end else
	        for(int row = 0;row < numberOfRows;row++){
	          for(int col = 0;col < numberOfCols;col++){
	            output[row][col] = temp[row][col];
	          }//col loop
	        }//row loop
	        return output;
	}

	static void histogram(double[][] inputData){
// Make this method run by the ImageJ histogram equalization		
		
	}
	static void selectROI(double[][] inputData){
		// This will approximate the minumum number of boxes to use to make
		// Square Regions of interest
				
			}
	
	public static void main(String[] args) {
		fiba x = new fiba();
		
		x.xform2D();
	}

}
