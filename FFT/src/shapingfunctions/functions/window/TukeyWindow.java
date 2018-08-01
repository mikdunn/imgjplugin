package shapingfunctions.functions.window;


import java.lang.Math.*;
import java.lang.Object.*;
import java.awt.*;
import java.awt.image.*;
import java.awt.image.BufferedImage;
import java.io.*;
import java.util.*;
import processing.core.PApplet;
import processing.core.PConstants;
import shapingfunctions.library.Function;


public class TukeyWindow extends Function {
	
	private float alpha;
	
	public TukeyWindow(){
		this(0);
	}
	
	public TukeyWindow(float alpha){
		super();
		
		this.alpha = alpha;
	}
	
	@Override
	public float applyFunction(float x, boolean clamp) {
		float ah = alpha / 2.0f;
		float omah = 1.0f - ah;
		
		float y = 1.0f;
		
		if (Float.compare(x, ah) <= 0){
			y = 0.5f * (1.0f + PApplet.cos(PConstants.PI * ((2*x / alpha) - 1.0f)));
		}
		else if (x > omah){
			y = 0.5f * (1.0f + PApplet.cos(PConstants.PI * ((2*x / alpha) - (2.0f/alpha) + 1.0f)));

		}
		
		return clamp(y, clamp);
	}
	public float getAlpha(){
		return alpha;
	}
	public void setAlpha(float alpha){
		this.alpha = alpha;
	}
}
