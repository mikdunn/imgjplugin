package fftanalysis.imagej;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.DirectoryChooser;
import ij.gui.GenericDialog;
import ij.gui.Plot;
import ij.io.FileInfo;
import ij.io.FileSaver;
import ij.measure.ResultsTable;
import ij.plugin.filter.PlugInFilter;
import ij.process.ByteProcessor;
import ij.process.ColorProcessor;
import ij.process.ImageProcessor;

import java.awt.Color;
import java.io.File;

/**
 * ImageJ plugin entry point.
 *
 * Runs a MATLAB-modeled FFT orientation analysis and reconstruction based on the
 * provided fiba.m reference.
 */
public class FIBA_Orientation implements PlugInFilter {

    private ImagePlus imp;

    @Override
    public int setup(String arg, ImagePlus imp) {
        this.imp = imp;
        return DOES_ALL;
    }

    @Override
    public void run(ImageProcessor ip) {
        if (imp == null || ip == null) {
            IJ.error("No image");
            return;
        }

        // Ensure grayscale
        ImageProcessor gray = ip;
        if (!(ip instanceof ByteProcessor)) {
            gray = ip.convertToByte(true);
        }

        // Center-crop to square if needed
        final int w0 = gray.getWidth();
        final int h0 = gray.getHeight();
        final int n = Math.min(w0, h0);
        if (n < 8) {
            IJ.error("Image too small");
            return;
        }
        final int x0 = (w0 - n) / 2;
        final int y0 = (h0 - n) / 2;
        gray.setRoi(x0, y0, n, n);
        gray = gray.crop();

        final int w = n / 2;

        final FibaMatlabProcessor.Params params = new FibaMatlabProcessor.Params();
        params.rmin = 4;
        params.rmax = w - 1;
        params.alpha = 0.4;
        params.beta = 0.3;
        params.gamma = 0.3;

        final OutputOptions outOpts = new OutputOptions();
        outOpts.saveOutputs = true;
        outOpts.showComposite = true;
        outOpts.showPlot = true;

        if (!showDialog(params, outOpts, w)) {
            return;
        }

        final double[][] j = toDouble(gray);

        final FibaMatlabProcessor.Result res;
        try {
            res = FibaMatlabProcessor.process(j, params);
        } catch (Exception e) {
            IJ.handleException(e);
            return;
        }

        // Output logging like MATLAB disp()
        IJ.log("Weighted average fiber: " + res.pAng + " degree");
        IJ.log("Width of the statistical significant peak: " + res.spWid + " degree");
        IJ.log("Strength of the significant peak: " + Math.round(res.bandStrength * 100.0) + "%");
        IJ.log("The 30% peak bandwidth is " + res.pWidth + " degree");
        if (res.warnPk == 1) {
            IJ.log("Warning, there might be more than one peak!");
        }
        IJ.log("Peak boundary: " + res.ang1 + "-" + res.ang2);

        // Results table row
        final ResultsTable rt = ResultsTable.getResultsTable();
        rt.incrementCounter();
        rt.addValue("pAng_deg", res.pAng);
        rt.addValue("spWid_deg", res.spWid);
        rt.addValue("bandStrength", res.bandStrength);
        rt.addValue("pWidth_deg", res.pWidth);
        rt.addValue("warnPk", res.warnPk);
        rt.addValue("ang1_deg", res.ang1);
        rt.addValue("ang2_deg", res.ang2);
        rt.show("Results");

        final String baseName = stripExtension(imp.getTitle());
        final String outputDir = determineOutputDir(imp, outOpts);

        ImagePlus composite = null;
        if (outOpts.showComposite || outOpts.saveOutputs) {
            composite = new ImagePlus(baseName + "_rec", buildComposite(res));
            if (outOpts.showComposite) composite.show();
        }

        ImagePlus plotImp = null;
        if (outOpts.showPlot || outOpts.saveOutputs) {
            plotImp = buildPlotImage(baseName + "_dat", res);
            if (outOpts.showPlot) plotImp.show();
        }

        if (outOpts.saveOutputs) {
            if (outputDir == null) {
                IJ.error("No output directory selected");
                return;
            }
            if (composite != null) {
                final File outFile = new File(outputDir, baseName + "_rec.jpg");
                new FileSaver(composite).saveAsJpeg(outFile.getAbsolutePath());
            }
            if (plotImp != null) {
                final File outFile = new File(outputDir, baseName + "_dat.jpg");
                new FileSaver(plotImp).saveAsJpeg(outFile.getAbsolutePath());
            }
        }
    }

    private static final class OutputOptions {
        boolean saveOutputs;
        boolean showComposite;
        boolean showPlot;
        String outputDirOverride;
    }

    private static boolean showDialog(FibaMatlabProcessor.Params params, OutputOptions out, int w) {
        final GenericDialog gd = new GenericDialog("FIBA (MATLAB-modeled)");
        gd.addMessage("Parameters are modeled after fiba.m");
        gd.addNumericField("rmin", params.rmin, 0);
        gd.addNumericField("rmax (<=0 for default w-1)", params.rmax, 0);
        gd.addNumericField("alpha (Tukey on image)", params.alpha, 3);
        gd.addNumericField("beta (radial mask Tukey)", params.beta, 3);
        gd.addNumericField("gamma (angular mask Tukey)", params.gamma, 3);

        gd.addMessage("Defaults: w=" + w + ", rmax=w-1");

        gd.addCheckbox("Show composite output", out.showComposite);
        gd.addCheckbox("Show orientation plot", out.showPlot);
        gd.addCheckbox("Save *_rec.jpg and *_dat.jpg", out.saveOutputs);

        gd.showDialog();
        if (gd.wasCanceled()) return false;

        params.rmin = (int) Math.round(gd.getNextNumber());
        params.rmax = (int) Math.round(gd.getNextNumber());
        params.alpha = gd.getNextNumber();
        params.beta = gd.getNextNumber();
        params.gamma = gd.getNextNumber();

        out.showComposite = gd.getNextBoolean();
        out.showPlot = gd.getNextBoolean();
        out.saveOutputs = gd.getNextBoolean();

        return true;
    }

    private static String determineOutputDir(ImagePlus imp, OutputOptions opts) {
        if (opts.outputDirOverride != null && opts.outputDirOverride.trim().length() > 0) {
            return opts.outputDirOverride;
        }

        final FileInfo fi = imp.getOriginalFileInfo();
        if (fi != null && fi.directory != null) {
            return fi.directory;
        }

        final DirectoryChooser dc = new DirectoryChooser("Choose output folder");
        final String dir = dc.getDirectory();
        return dir;
    }

    private static String stripExtension(String name) {
        if (name == null) return "image";
        int dot = name.lastIndexOf('.');
        if (dot > 0) return name.substring(0, dot);
        return name;
    }

    private static double[][] toDouble(ImageProcessor ip) {
        final int w = ip.getWidth();
        final int h = ip.getHeight();
        final double[][] out = new double[h][w];
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                out[y][x] = ip.getPixelValue(x, y);
            }
        }
        return out;
    }

    private static ColorProcessor buildComposite(FibaMatlabProcessor.Result res) {
        final int n = res.n;
        final int gap = 5;
        final int width = 5 * n + 4 * gap;
        final int height = n;
        final ColorProcessor cp = new ColorProcessor(width, height);
        cp.setColor(Color.white);
        cp.fill();

        int x = 0;
        blitGray(cp, res.origNorm, x, 0);
        x += n + gap;
        blitGray(cp, res.imgS, x, 0);
        x += n + gap;
        blitGray(cp, res.imgFDisp, x, 0);
        x += n + gap;
        blitGray(cp, res.imgR2, x, 0);
        x += n + gap;
        blitRGB(cp, res.overlay, x, 0);

        return cp;
    }

    private static void blitGray(ColorProcessor cp, double[][] img01, int x0, int y0) {
        final int h = img01.length;
        final int w = img01[0].length;
        final int[] rgb = new int[3];
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                int v = (int) Math.round(clamp01(img01[y][x]) * 255.0);
                rgb[0] = v;
                rgb[1] = v;
                rgb[2] = v;
                cp.putPixel(x0 + x, y0 + y, rgb);
            }
        }
    }

    private static void blitRGB(ColorProcessor cp, double[][][] img01, int x0, int y0) {
        final int h = img01.length;
        final int w = img01[0].length;
        final int[] rgb = new int[3];
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                rgb[0] = (int) Math.round(clamp01(img01[y][x][0]) * 255.0);
                rgb[1] = (int) Math.round(clamp01(img01[y][x][1]) * 255.0);
                rgb[2] = (int) Math.round(clamp01(img01[y][x][2]) * 255.0);
                cp.putPixel(x0 + x, y0 + y, rgb);
            }
        }
    }

    private static double clamp01(double v) {
        if (v < 0) return 0;
        if (v > 1) return 1;
        return v;
    }

    private static ImagePlus buildPlotImage(String title, FibaMatlabProcessor.Result res) {
        final double[] x = new double[180];
        final double[] y = new double[180];
        for (int i = 0; i < 180; i++) {
            x[i] = i;
            y[i] = res.sol[i];
        }

        final Plot plot = new Plot(title, "Fibril Orientation (deg)", "Signal Strength", x, y);
        plot.setLineWidth(2);

        // mean line
        plot.setColor(Color.red);
        plot.add("line", x, constantArray(180, res.meanSol));

        // mean + std line
        plot.setColor(Color.green.darker());
        plot.add("line", x, constantArray(180, res.meanSol + res.stdSol));

        // highlight reconstruction bandwidth (similar to MATLAB aind1/aind2)
        plot.setColor(Color.red);
        plot.setLineWidth(4);
        if (res.aind1 != null) {
            addHighlight(plot, res.aind1, res);
        }
        if (res.aind2 != null && res.aind2.length > 1) {
            addHighlight(plot, res.aind2, res);
        }

        return plot.getImagePlus();
    }

    private static void addHighlight(Plot plot, int[] aind, FibaMatlabProcessor.Result res) {
        final int len = aind.length;
        final double[] hx = new double[len];
        final double[] hy = new double[len];
        for (int i = 0; i < len; i++) {
            int ang = aind[i] - 1; // MATLAB uses aind-1 on x-axis
            if (ang < 0) ang = 0;
            if (ang > 179) ang = 179;
            hx[i] = ang;
            hy[i] = res.sol[ang];
        }
        plot.add("line", hx, hy);
    }

    private static double[] constantArray(int n, double v) {
        final double[] a = new double[n];
        for (int i = 0; i < n; i++) a[i] = v;
        return a;
    }
}
