package fftanalysis.imagej;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.gui.Plot;
import ij.Macro;
import ij.io.FileInfo;
import ij.io.DirectoryChooser;
import ij.io.FileSaver;
import ij.measure.ResultsTable;
import ij.plugin.filter.PlugInFilter;
import ij.process.ByteProcessor;
import ij.process.ColorProcessor;
import ij.process.ImageProcessor;

import java.awt.Color;
import java.awt.GraphicsEnvironment;
import java.io.File;
import java.io.FileWriter;
import java.util.HashMap;
import java.util.Map;

/**
 * ImageJ plugin entry point.
 *
 * Runs a MATLAB-modeled FFT orientation analysis and reconstruction based on the
 * provided fiba.m reference.
 */
public class FIBA_Orientation implements PlugInFilter {

    private ImagePlus imp;
    private String argOptions;

    @Override
    public int setup(String arg, ImagePlus imp) {
        this.imp = imp;
        this.argOptions = arg;
        return DOES_ALL;
    }

    @Override
    public void run(ImageProcessor ip) {
        try {
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
        // Suppress an unnaturally sharp spike at exactly 90deg (common axial artifact) before peak+mask.
        params.suppressAngleSpike = true;
        params.suppressAngleDeg = 90;
        params.suppressHalfWidthDeg = 0;
        params.suppressIfOverMedianRatio = 6.0;

        // Mask-stage artifact suppression: remove any perfectly-vertical column that spans the full image height.
        params.removeFullHeightVerticalLine = true;
        params.removeFullHeightVerticalLineMinCoverage = 1.0;

        final OutputOptions outOpts = new OutputOptions();
        outOpts.saveOutputs = true;
        outOpts.showComposite = true;
        outOpts.showPlot = true;

        // Allow non-interactive execution from macros / batch mode.
        // Priority: explicit argOptions (from setup) then Macro options.
        String macroOpts = (argOptions != null && argOptions.trim().length() > 0)
                ? argOptions
                : Macro.getOptions();
        if (macroOpts == null || macroOpts.trim().isEmpty()) {
            final String sysOpts = System.getProperty("fiba.options");
            if (sysOpts != null && sysOpts.trim().length() > 0) {
                macroOpts = sysOpts;
            }
        }

        // In ImageJ "-batch" runs, GraphicsEnvironment is typically *not* headless,
        // but IJ.getInstance() can still be null. Use a separate flag for UI safety.
        final boolean headless = GraphicsEnvironment.isHeadless();
        final boolean interactiveUI = !headless && IJ.getInstance() != null;
        if (macroOpts != null && macroOpts.trim().length() > 0) {
            applyOptions(macroOpts, params, outOpts, w);
        }

        final String baseName = stripExtension(imp.getTitle());
        final String outputDir = determineOutputDir(imp, outOpts, interactiveUI);
        if (outOpts.debug) {
            debugToFile(outOpts, outputDir, "[FIBA] start interactiveUI=" + interactiveUI + " headless=" + headless);
            debugToFile(outOpts, outputDir, "[FIBA] opts='" + (macroOpts == null ? "" : macroOpts) + "'");
        }

        if (outOpts.debug) {
            // Console output is unreliable in some ImageJ batch launches; keep it, but also log to file.
            System.out.println("[FIBA] interactiveUI=" + interactiveUI + " headless=" + headless + " opts='" + (macroOpts == null ? "" : macroOpts) + "'");
        }

        if (interactiveUI && (macroOpts == null || macroOpts.trim().isEmpty())) {
            if (!showDialog(params, outOpts, w)) {
                return;
            }
        } else {
            // In batch/headless mode always save outputs unless explicitly disabled.
            // (Macro options can override.)
            // Also suppress windows.
            if (macroOpts == null || !macroOpts.toLowerCase().contains("save=")) {
                outOpts.saveOutputs = true;
            }
            if (macroOpts == null || !macroOpts.toLowerCase().contains("showcomposite=")) {
                outOpts.showComposite = false;
            }
            if (macroOpts == null || !macroOpts.toLowerCase().contains("showplot=")) {
                outOpts.showPlot = false;
            }
        }

        final double[][] j = toDouble(gray);

        if (outOpts.saveOutputs && (outputDir == null || outputDir.trim().isEmpty())) {
            // Never prompt for a directory in batch/non-interactive mode.
            IJ.error("No output directory available (set outputDir=... in options)");
            debugToFile(outOpts, outputDir, "[FIBA] ERROR: outputDir is null/empty");
            return;
        }

        final FibaMatlabProcessor.Result res;
        try {
            debugToFile(outOpts, outputDir, "[FIBA] processing start n=" + n);
            res = FibaMatlabProcessor.process(j, params);
            debugToFile(outOpts, outputDir, "[FIBA] processing done pAng=" + res.pAng);
        } catch (Exception e) {
            if (outOpts.debug || headless) {
                e.printStackTrace();
            }
            debugToFile(outOpts, outputDir, "[FIBA] EXCEPTION: " + e);
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
        if (interactiveUI) {
            rt.show("Results");
        }

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
            if (composite != null) {
                final File outFile = new File(outputDir, baseName + "_rec.jpg");
                new FileSaver(composite).saveAsJpeg(outFile.getAbsolutePath());
                debugToFile(outOpts, outputDir, "[FIBA] saved: " + outFile.getAbsolutePath());
            }
            if (plotImp != null) {
                final File outFile = new File(outputDir, baseName + "_dat.jpg");
                new FileSaver(plotImp).saveAsJpeg(outFile.getAbsolutePath());
                debugToFile(outOpts, outputDir, "[FIBA] saved: " + outFile.getAbsolutePath());
            }
        }

        debugToFile(outOpts, outputDir, "[FIBA] done");

        } catch (Throwable t) {
            t.printStackTrace();
            try {
                // Best-effort: log to a file if we can infer an output directory.
                String outDir = null;
                try {
                    final String sysOpts = System.getProperty("fiba.options");
                    if (sysOpts != null) {
                        final Map<String, String> kv = parseKeyValueOptions(sysOpts);
                        final String d = kv.get("outputdir");
                        if (d != null && d.trim().length() > 0) outDir = d.trim();
                    }
                } catch (Throwable ignore) {
                    // ignore
                }
                if (outDir == null && imp != null) {
                    final FileInfo fi = imp.getOriginalFileInfo();
                    if (fi != null) outDir = fi.directory;
                }
                debugToFileAlways(outDir, "[FIBA] FATAL: " + t);
            } catch (Throwable ignore) {
                // ignore
            }
            IJ.handleException(t);
        }
    }

    private static final class OutputOptions {
        boolean saveOutputs;
        boolean showComposite;
        boolean showPlot;
        String outputDirOverride;
        boolean debug;
    }

    private static void applyOptions(String opts, FibaMatlabProcessor.Params params, OutputOptions out, int w) {
        final Map<String, String> kv = parseKeyValueOptions(opts);

        params.rmin = parseInt(kv.get("rmin"), params.rmin);
        params.rmax = parseInt(kv.get("rmax"), params.rmax);
        if (params.rmax <= 0) params.rmax = w - 1;

        params.alpha = parseDouble(kv.get("alpha"), params.alpha);
        params.beta = parseDouble(kv.get("beta"), params.beta);
        params.gamma = parseDouble(kv.get("gamma"), params.gamma);

        // Artifact suppression at an exact angle before peak+mask.
        params.suppressAngleSpike = parseBoolean(firstNonNull(kv.get("suppressanglespike"), kv.get("suppress90"), kv.get("removeninety")), params.suppressAngleSpike);
        params.suppressAngleDeg = parseInt(firstNonNull(kv.get("suppressangledeg"), kv.get("suppressangle"), kv.get("suppresstheta")), params.suppressAngleDeg);
        params.suppressHalfWidthDeg = parseInt(firstNonNull(kv.get("suppresshalfwidthdeg"), kv.get("suppresswidth"), kv.get("suppressw")), params.suppressHalfWidthDeg);
        params.suppressIfOverMedianRatio = parseDouble(firstNonNull(kv.get("suppressifovermedianratio"), kv.get("suppressratio"), kv.get("suppressr")), params.suppressIfOverMedianRatio);

        params.removeFullHeightVerticalLine = parseBoolean(
            firstNonNull(kv.get("removefullheightverticalline"), kv.get("removeverticalfullline"), kv.get("removeverticalcolumn")),
            params.removeFullHeightVerticalLine);
        params.removeFullHeightVerticalLineMinCoverage = parseDouble(
            firstNonNull(kv.get("removefullheightverticallinemincoverage"), kv.get("removeverticalmincoverage"), kv.get("verticalmincoverage")),
            params.removeFullHeightVerticalLineMinCoverage);

        out.saveOutputs = parseBoolean(kv.get("save"), out.saveOutputs);
        out.showComposite = parseBoolean(kv.get("showcomposite"), out.showComposite);
        out.showPlot = parseBoolean(kv.get("showplot"), out.showPlot);
        out.debug = parseBoolean(kv.get("debug"), out.debug);

        final String outDir = kv.get("outputdir");
        if (outDir != null && outDir.trim().length() > 0) {
            out.outputDirOverride = outDir.trim();
        }
    }

    private static String firstNonNull(String a, String b, String c) {
        if (a != null) return a;
        if (b != null) return b;
        return c;
    }

    private static Map<String, String> parseKeyValueOptions(String opts) {
        // ImageJ macro options are typically "key=value key2=value2".
        // Keys are case-insensitive.
        final Map<String, String> map = new HashMap<String, String>();
        if (opts == null) return map;
        final String[] tokens = opts.trim().split("\\s+");
        for (int i = 0; i < tokens.length; i++) {
            final String t = tokens[i];
            final int eq = t.indexOf('=');
            if (eq <= 0) continue;
            final String k = t.substring(0, eq).trim().toLowerCase();
            final String v = t.substring(eq + 1).trim();
            if (k.length() == 0) continue;
            map.put(k, v);
        }
        return map;
    }

    private static int parseInt(String s, int def) {
        if (s == null) return def;
        try {
            return (int) Math.round(Double.parseDouble(s));
        } catch (Exception e) {
            return def;
        }
    }

    private static double parseDouble(String s, double def) {
        if (s == null) return def;
        try {
            return Double.parseDouble(s);
        } catch (Exception e) {
            return def;
        }
    }

    private static boolean parseBoolean(String s, boolean def) {
        if (s == null) return def;
        final String v = s.trim().toLowerCase();
        if (v.equals("true") || v.equals("1") || v.equals("yes") || v.equals("y")) return true;
        if (v.equals("false") || v.equals("0") || v.equals("no") || v.equals("n")) return false;
        return def;
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

    private static String determineOutputDir(ImagePlus imp, OutputOptions opts, boolean interactiveUI) {
        if (opts.outputDirOverride != null && opts.outputDirOverride.trim().length() > 0) {
            return opts.outputDirOverride;
        }

        final FileInfo fi = imp.getOriginalFileInfo();
        if (fi != null && fi.directory != null) {
            return fi.directory;
        }

        if (!interactiveUI) {
            return null;
        }

        final DirectoryChooser dc = new DirectoryChooser("Choose output folder");
        return dc.getDirectory();
    }

    private static void debugToFile(OutputOptions out, String outputDir, String message) {
        if (out == null || !out.debug) return;
        debugToFileAlways(outputDir, message);
    }

    private static void debugToFileAlways(String outputDir, String message) {
        if (outputDir == null || outputDir.trim().isEmpty()) return;
        FileWriter fw = null;
        try {
            final File f = new File(outputDir, "fiba_plugin_debug.txt");
            fw = new FileWriter(f, true);
            fw.write(message);
            fw.write(System.lineSeparator());
        } catch (Exception ignore) {
            // ignore
        } finally {
            if (fw != null) {
                try { fw.close(); } catch (Exception ignore) { /* ignore */ }
            }
        }
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
