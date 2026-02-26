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

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.GraphicsEnvironment;
import java.awt.RenderingHints;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileWriter;
import java.util.HashMap;
import java.util.Locale;
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
        // Artifact suppression is OFF by default. We will deal with artifacts ONLY with alpha values.
        params.suppressAngleSpike = false;
        params.removeFullHeightVerticalLine = false;
        params.removeFullWidthHorizontalLine = false;

        final OutputOptions outOpts = new OutputOptions();
        outOpts.saveOutputs = true;
        outOpts.saveSolCsv = true;
        outOpts.showComposite = true;
        outOpts.showPlot = true;
        outOpts.plotApplyTukey = false;
        // Fixed polar scaling so amplitude differences remain visible across alpha sweeps.
        // Typical SOL peaks are on the order of 0.01–0.03.
        outOpts.polarScaleMax = 0.03;

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
            final String lower = (macroOpts == null) ? "" : macroOpts.toLowerCase();
            if (!(lower.contains("savesolcsv=") || lower.contains("savesol=") || lower.contains("writesolcsv="))) {
                outOpts.saveSolCsv = outOpts.saveOutputs;
            }
            if (macroOpts == null || !macroOpts.toLowerCase().contains("showcomposite=")) {
                outOpts.showComposite = false;
            }
            if (macroOpts == null || !macroOpts.toLowerCase().contains("showplot=")) {
                outOpts.showPlot = false;
            }
        }

        final double[][] j = toDouble(gray);

        if ((outOpts.saveOutputs || outOpts.saveSolCsv) && (outputDir == null || outputDir.trim().isEmpty())) {
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
            debugToFile(outOpts, outputDir, String.format(Locale.US,
                    "[FIBA] metrics pAng=%d spWid=%d pWidth=%d bandStrength=%.6f warnPk=%d ang1=%d ang2=%d meanSol=%.10g stdSol=%.10g",
                    res.pAng, res.spWid, res.pWidth, res.bandStrength, res.warnPk, res.ang1, res.ang2, res.meanSol, res.stdSol));
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
            final double[] solPlot = solForPlot(res, params, outOpts);
            plotImp = buildPlotImage(baseName + "_dat", res, solPlot);
            if (outOpts.showPlot) plotImp.show();
        }

        ImagePlus polarImp = null;
        if (outOpts.showPlot || outOpts.saveOutputs) {
            final double[] solPlot = solForPlot(res, params, outOpts);
            polarImp = buildPolarPlotImage(baseName + "_polar", res, solPlot, outOpts.polarScaleMax);
            // Keep behavior simple: when user requests plots, show both.
            if (outOpts.showPlot) polarImp.show();
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
            if (polarImp != null) {
                final File outFile = new File(outputDir, baseName + "_polar.jpg");
                new FileSaver(polarImp).saveAsJpeg(outFile.getAbsolutePath());
                debugToFile(outOpts, outputDir, "[FIBA] saved: " + outFile.getAbsolutePath());
            }
        }

        if (outOpts.saveSolCsv) {
            try {
                final File outFile = new File(outputDir, baseName + "_sol.csv");
                writeSolCsv(outFile, res, params);
                debugToFile(outOpts, outputDir, "[FIBA] saved: " + outFile.getAbsolutePath());
            } catch (Exception e) {
                debugToFile(outOpts, outputDir, "[FIBA] WARN: failed to save SOL CSV: " + e);
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
        boolean saveSolCsv;
        boolean showComposite;
        boolean showPlot;
        boolean plotApplyTukey;
        double polarScaleMax;
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
        params.suppressAnglesDeg = parseIntList(firstNonNull(kv.get("suppressangles"), kv.get("suppressanglelist"), kv.get("suppresslist")), params.suppressAnglesDeg);
        params.suppressHalfWidthDeg = parseInt(firstNonNull(kv.get("suppresshalfwidthdeg"), kv.get("suppresswidth"), kv.get("suppressw")), params.suppressHalfWidthDeg);
        params.suppressIfOverMedianRatio = parseDouble(firstNonNull(kv.get("suppressifovermedianratio"), kv.get("suppressratio"), kv.get("suppressr")), params.suppressIfOverMedianRatio);

        params.removeFullHeightVerticalLine = parseBoolean(
            firstNonNull(kv.get("removefullheightverticalline"), kv.get("removeverticalfullline"), kv.get("removeverticalcolumn")),
            params.removeFullHeightVerticalLine);
        params.removeFullHeightVerticalLineMinCoverage = parseDouble(
            firstNonNull(kv.get("removefullheightverticallinemincoverage"), kv.get("removeverticalmincoverage"), kv.get("verticalmincoverage")),
            params.removeFullHeightVerticalLineMinCoverage);

        params.removeFullWidthHorizontalLine = parseBoolean(
            firstNonNull(kv.get("removefullwidthhorizontalline"), kv.get("removehorizontalfullline"), kv.get("removehorizontalrow")),
            params.removeFullWidthHorizontalLine);
        params.removeFullWidthHorizontalLineMinCoverage = parseDouble(
            firstNonNull(kv.get("removefullwidthhorizontallinemincoverage"), kv.get("removehorizontalmincoverage"), kv.get("horizontalmincoverage")),
            params.removeFullWidthHorizontalLineMinCoverage);

        out.saveOutputs = parseBoolean(kv.get("save"), out.saveOutputs);
        out.saveSolCsv = parseBoolean(firstNonNull(kv.get("savesolcsv"), kv.get("savesol"), kv.get("writesolcsv")), out.saveSolCsv);
        out.showComposite = parseBoolean(kv.get("showcomposite"), out.showComposite);
        out.showPlot = parseBoolean(kv.get("showplot"), out.showPlot);
        out.plotApplyTukey = parseBoolean(firstNonNull(kv.get("plotapplytukey"), kv.get("applytukeytoplot"), kv.get("tukeyplot")), out.plotApplyTukey);
        out.polarScaleMax = parseDouble(firstNonNull(kv.get("polarscalemax"), kv.get("polarrmax"), kv.get("polarscale")), out.polarScaleMax);
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

    private static int[] parseIntList(String s, int[] def) {
        if (s == null) return def;
        final String t = s.trim();
        if (t.isEmpty()) return def;
        final String[] parts = t.split(",");
        int[] out = new int[parts.length];
        int n = 0;
        for (int i = 0; i < parts.length; i++) {
            final String p = parts[i].trim();
            if (p.isEmpty()) continue;
            try {
                out[n++] = (int) Math.round(Double.parseDouble(p));
            } catch (Exception ignore) {
                // ignore token
            }
        }
        if (n == 0) return def;
        if (n == out.length) return out;
        int[] trimmed = new int[n];
        System.arraycopy(out, 0, trimmed, 0, n);
        return trimmed;
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
        gd.addCheckbox("Save *_sol.csv", out.saveSolCsv);

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
        out.saveSolCsv = gd.getNextBoolean();

        return true;
    }

    private static void writeSolCsv(File outFile, FibaMatlabProcessor.Result res, FibaMatlabProcessor.Params params) throws Exception {
        if (outFile == null || res == null || res.sol == null || res.sol.length < 180) return;

        final int rminUsed = (params == null) ? -1 : Math.max(0, params.rmin);
        final int rmaxUsed = (params == null)
                ? -1
                : ((params.rmax > 0) ? params.rmax : (res.w - 1));

        try (FileWriter fw = new FileWriter(outFile, false)) {
            fw.write("angle_deg,sol,meanSol,stdSol,pAng_deg,spWid_deg,bandStrength,pWidth_deg,warnPk,ang1_deg,ang2_deg,n,rmin,rmax,alpha,beta,gamma");
            fw.write(System.lineSeparator());
            for (int ang = 0; ang < 180; ang++) {
                fw.write(Integer.toString(ang));
                fw.write(',');
                fw.write(Double.toString(res.sol[ang]));
                fw.write(',');
                fw.write(Double.toString(res.meanSol));
                fw.write(',');
                fw.write(Double.toString(res.stdSol));
                fw.write(',');
                fw.write(Integer.toString(res.pAng));
                fw.write(',');
                fw.write(Integer.toString(res.spWid));
                fw.write(',');
                fw.write(Double.toString(res.bandStrength));
                fw.write(',');
                fw.write(Integer.toString(res.pWidth));
                fw.write(',');
                fw.write(Integer.toString(res.warnPk));
                fw.write(',');
                fw.write(Integer.toString(res.ang1));
                fw.write(',');
                fw.write(Integer.toString(res.ang2));
                fw.write(',');
                fw.write(Integer.toString(res.n));
                fw.write(',');
                fw.write(Integer.toString(rminUsed));
                fw.write(',');
                fw.write(Integer.toString(rmaxUsed));
                fw.write(',');
                fw.write(Double.toString(params == null ? Double.NaN : params.alpha));
                fw.write(',');
                fw.write(Double.toString(params == null ? Double.NaN : params.beta));
                fw.write(',');
                fw.write(Double.toString(params == null ? Double.NaN : params.gamma));
                fw.write(System.lineSeparator());
            }
        }
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

    private static ImagePlus buildPlotImage(String title, FibaMatlabProcessor.Result res, double[] solForPlot) {
        final double[] x = new double[180];
        final double[] y = new double[180];
        for (int i = 0; i < 180; i++) {
            x[i] = i;
            y[i] = (solForPlot != null && solForPlot.length >= 180) ? solForPlot[i] : res.sol[i];
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

    private static ImagePlus buildPolarPlotImage(String title, FibaMatlabProcessor.Result res, double[] solForPlot, double polarScaleMax) {
        final int size = 600;
        final int margin = 40;
        final int cx = size / 2;
        final int cy = size / 2;
        final int maxR = (size / 2) - margin;

        // Replicate 180-degree SOL to 360 for a full-circle visualization.
        // Use a fixed scale so alpha sweeps are visually comparable.
        final double[] s = (solForPlot != null && solForPlot.length >= 180) ? solForPlot : res.sol;
        double scaleMax = polarScaleMax;
        if (!(scaleMax > 0)) scaleMax = 0.03;

        final BufferedImage bi = new BufferedImage(size, size, BufferedImage.TYPE_INT_RGB);
        final Graphics2D g = bi.createGraphics();
        try {
            g.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
            g.setColor(Color.white);
            g.fillRect(0, 0, size, size);

            // Grid
            g.setColor(new Color(220, 220, 220));
            g.setStroke(new BasicStroke(1f));
            for (int k = 1; k <= 4; k++) {
                final int r = (int) Math.round(maxR * (k / 4.0));
                g.drawOval(cx - r, cy - r, 2 * r, 2 * r);
            }
            g.drawLine(cx - maxR, cy, cx + maxR, cy);
            g.drawLine(cx, cy - maxR, cx, cy + maxR);

            // Label orientation (matches our theta-from-vertical convention)
            g.setColor(Color.darkGray);
            // Standard polar convention: 0° at top, 90° at right, increasing clockwise.
            g.drawString("0°", cx - 8, cy - maxR - 8);
            g.drawString("90°", cx + maxR + 6, cy + 4);
            g.drawString("180°", cx - 18, cy + maxR + 16);
            g.drawString("270°", cx - maxR - 28, cy + 4);

            // SOL polyline
            g.setColor(new Color(0, 102, 204));
            g.setStroke(new BasicStroke(2.5f));
            int prevX = Integer.MIN_VALUE;
            int prevY = Integer.MIN_VALUE;
            for (int deg = 0; deg < 360; deg++) {
                final double v = s[deg % 180];
                double rr = maxR * (v / scaleMax);
                if (rr < 0) rr = 0;
                if (rr > maxR) rr = maxR;

                final double theta = Math.toRadians(deg);
                // x is column, y is row; theta=0 points up (negative row) for standard polar plots.
                final int x = (int) Math.round(cx + rr * Math.sin(theta));
                final int y = (int) Math.round(cy - rr * Math.cos(theta));

                if (prevX != Integer.MIN_VALUE) {
                    g.drawLine(prevX, prevY, x, y);
                }
                prevX = x;
                prevY = y;
            }

            // Center marker
            g.setColor(Color.black);
            g.fillOval(cx - 2, cy - 2, 5, 5);
        } finally {
            g.dispose();
        }

        return new ImagePlus(title, new ColorProcessor(bi));
    }

    private static double[] solForPlot(FibaMatlabProcessor.Result res, FibaMatlabProcessor.Params params, OutputOptions out) {
        if (res == null || res.sol == null || res.sol.length < 180) return null;
        if (out == null || !out.plotApplyTukey) return res.sol;

        final double alpha = (params == null) ? 0.0 : params.alpha;
        final double[] w = tukeyEdgeVectorPlot(180, 90, alpha);
        final double[] y = new double[180];
        for (int i = 0; i < 180; i++) {
            y[i] = res.sol[i] * w[i];
        }
        return y;
    }

    private static double[] tukeyEdgeVectorPlot(int n, int w, double alpha) {
        // Mirror the Tukey edge vector logic used in the core processor (MATLAB-faithful).
        if (alpha <= 0) {
            final double[] s = new double[n];
            for (int i = 0; i < n; i++) s[i] = 1;
            return s;
        }

        final double[] s = new double[n];
        for (int i1 = 1; i1 <= n; i1++) {
            final double i = i1;
            final double v;
            if (i <= alpha * w) {
                v = 0.5 * (1.0 + Math.cos(Math.PI * (((i - 1.0) / alpha / w) - 1.0)));
            } else if (i <= 2.0 * w * (1.0 - alpha / 2.0)) {
                v = 1.0;
            } else {
                v = 0.5 * (1.0 + Math.cos(Math.PI * ((i / alpha / w) - (2.0 / alpha) - 1.0)));
            }
            s[i1 - 1] = v;
        }
        return s;
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
