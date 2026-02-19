package fftanalysis.imagej;

import ij.IJ;
import ij.ImagePlus;
import ij.Macro;
import ij.gui.GenericDialog;
import ij.io.DirectoryChooser;
import ij.plugin.PlugIn;
import ij.process.ByteProcessor;
import ij.process.ImageProcessor;

import java.awt.GraphicsEnvironment;
import java.io.File;
import java.io.FileWriter;
import java.util.HashMap;
import java.util.Map;

/**
 * MATLAB fibaall.m equivalent.
 *
 * Processes a folder containing numbered cropbox images (1.jpg, 2.jpg, ...),
 * runs the MATLAB-modeled fiba pipeline per image, and writes:
 *  - SOL CSV: first row 0..179, then one row per image with sol[0..179]
 *  - pAng CSV: one value per image (wrapped to [-90,90] like fibaall.m)
 */
public class FIBA_All_FromFolder implements PlugIn {

    @Override
    public void run(String arg) {
        final boolean headless = GraphicsEnvironment.isHeadless();
        final boolean interactiveUI = !headless && IJ.getInstance() != null;

        Options opts = new Options();
        opts.count = 10;
        opts.saveOutputs = true;
        opts.wrap90 = true;
        opts.debug = false;

        FibaMatlabProcessor.Params params = new FibaMatlabProcessor.Params();
        params.rmin = 4;
        params.rmax = 0;
        params.alpha = 0.4;
        params.beta = 0.3;
        params.gamma = 0.3;
        // Suppress an unnaturally sharp spike at exactly 90deg (common axial artifact) before peak+mask.
        params.suppressAngleSpike = true;
        params.suppressAngleDeg = 90;
        params.suppressHalfWidthDeg = 0;
        params.suppressIfOverMedianRatio = 6.0;

        String macroOpts = (arg != null && arg.trim().length() > 0) ? arg : Macro.getOptions();
        if (macroOpts == null || macroOpts.trim().isEmpty()) {
            final String sys = System.getProperty("fiba.options");
            if (sys != null && sys.trim().length() > 0) macroOpts = sys;
        }
        if (macroOpts != null && macroOpts.trim().length() > 0) {
            applyOptions(macroOpts, opts, params);
        }

        if (interactiveUI && (opts.dir == null || opts.dir.trim().isEmpty())) {
            DirectoryChooser dc = new DirectoryChooser("Choose folder with 1.jpg, 2.jpg, ...");
            String d = dc.getDirectory();
            if (d == null) return;
            opts.dir = d;

            GenericDialog gd = new GenericDialog("FIBA All (MATLAB-modeled)");
            gd.addNumericField("How many overlapping cropbox? (count)", opts.count, 0);
            gd.addCheckbox("Wrap pAng to [-90,90]", opts.wrap90);
            gd.addCheckbox("Save SOL + pAng CSV", opts.saveOutputs);
            gd.addCheckbox("Debug", opts.debug);
            gd.showDialog();
            if (gd.wasCanceled()) return;
            opts.count = (int) Math.round(gd.getNextNumber());
            opts.wrap90 = gd.getNextBoolean();
            opts.saveOutputs = gd.getNextBoolean();
            opts.debug = gd.getNextBoolean();
        }

        if (opts.dir == null || opts.dir.trim().isEmpty()) {
            IJ.error("No directory specified (use dir=...)");
            return;
        }
        if (opts.count <= 0) {
            IJ.error("count must be > 0");
            return;
        }

        File dir = new File(opts.dir);
        if (!dir.exists() || !dir.isDirectory()) {
            IJ.error("Not a directory: " + opts.dir);
            return;
        }

        final double[][] solRows = new double[opts.count][180];
        final double[] pAngList = new double[opts.count];

        for (int i = 1; i <= opts.count; i++) {
            File f = new File(dir, i + ".jpg");
            if (!f.exists()) {
                // Try png as a convenience
                File p = new File(dir, i + ".png");
                if (p.exists()) f = p;
            }
            if (!f.exists()) {
                IJ.log("[FIBA_ALL] Missing: " + f.getAbsolutePath());
                pAngList[i - 1] = Double.NaN;
                for (int k = 0; k < 180; k++) solRows[i - 1][k] = Double.NaN;
                continue;
            }

            ImagePlus imp = IJ.openImage(f.getAbsolutePath());
            if (imp == null) {
                IJ.log("[FIBA_ALL] Failed to open: " + f.getAbsolutePath());
                pAngList[i - 1] = Double.NaN;
                for (int k = 0; k < 180; k++) solRows[i - 1][k] = Double.NaN;
                continue;
            }

            ImageProcessor ip = imp.getProcessor();
            if (!(ip instanceof ByteProcessor)) {
                ip = ip.convertToByte(true);
            }

            // Crop to square (center) like the single-image plugin, in case crops aren't perfectly square.
            final int w0 = ip.getWidth();
            final int h0 = ip.getHeight();
            final int n = Math.min(w0, h0);
            final int x0 = (w0 - n) / 2;
            final int y0 = (h0 - n) / 2;
            ip.setRoi(x0, y0, n, n);
            ip = ip.crop();

            // Per-image defaults
            final int w = n / 2;
            if (params.rmax <= 0 || params.rmax > (w - 1)) params.rmax = w - 1;

            double[][] j = toDouble(ip);

            FibaMatlabProcessor.Result res;
            try {
                res = FibaMatlabProcessor.process(j, params);
            } catch (Throwable t) {
                IJ.log("[FIBA_ALL] Error on " + f.getName() + ": " + t);
                pAngList[i - 1] = Double.NaN;
                for (int k = 0; k < 180; k++) solRows[i - 1][k] = Double.NaN;
                continue;
            }

            double pAng = res.pAng;
            if (opts.wrap90 && pAng > 90) pAng = pAng - 180;
            pAngList[i - 1] = pAng;

            if (res.sol != null && res.sol.length >= 180) {
                for (int k = 0; k < 180; k++) solRows[i - 1][k] = res.sol[k];
            }
        }

        if (opts.saveOutputs) {
            File solFile = new File(dir, "SOL.csv");
            writeSolCsv(solFile, solRows);
            File pangFile = new File(dir, "pAng.csv");
            writeColumnCsv(pangFile, "pAng_adj_deg", pAngList);
            IJ.log("[FIBA_ALL] Wrote: " + solFile.getAbsolutePath());
            IJ.log("[FIBA_ALL] Wrote: " + pangFile.getAbsolutePath());
        }
    }

    private static final class Options {
        String dir;
        int count;
        boolean saveOutputs;
        boolean wrap90;
        boolean debug;
    }

    private static void applyOptions(String opts, Options out, FibaMatlabProcessor.Params params) {
        final Map<String, String> kv = parseKeyValueOptions(opts);

        final String dir = kv.get("dir");
        if (dir != null && dir.trim().length() > 0) out.dir = dir.trim();

        out.count = parseInt(kv.get("count"), out.count);
        out.saveOutputs = parseBoolean(kv.get("save"), out.saveOutputs);
        out.wrap90 = parseBoolean(kv.get("wrap90"), out.wrap90);
        out.debug = parseBoolean(kv.get("debug"), out.debug);

        params.rmin = parseInt(kv.get("rmin"), params.rmin);
        params.rmax = parseInt(kv.get("rmax"), params.rmax);
        params.alpha = parseDouble(kv.get("alpha"), params.alpha);
        params.beta = parseDouble(kv.get("beta"), params.beta);
        params.gamma = parseDouble(kv.get("gamma"), params.gamma);

        // Artifact suppression at an exact angle before peak+mask.
        params.suppressAngleSpike = parseBoolean(firstNonNull(kv.get("suppressanglespike"), kv.get("suppress90"), kv.get("removeninety")), params.suppressAngleSpike);
        params.suppressAngleDeg = parseInt(firstNonNull(kv.get("suppressangledeg"), kv.get("suppressangle"), kv.get("suppresstheta")), params.suppressAngleDeg);
        params.suppressHalfWidthDeg = parseInt(firstNonNull(kv.get("suppresshalfwidthdeg"), kv.get("suppresswidth"), kv.get("suppressw")), params.suppressHalfWidthDeg);
        params.suppressIfOverMedianRatio = parseDouble(firstNonNull(kv.get("suppressifovermedianratio"), kv.get("suppressratio"), kv.get("suppressr")), params.suppressIfOverMedianRatio);
    }

    private static String firstNonNull(String a, String b, String c) {
        if (a != null) return a;
        if (b != null) return b;
        return c;
    }

    private static Map<String, String> parseKeyValueOptions(String opts) {
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

    private static void writeSolCsv(File file, double[][] solRows) {
        FileWriter fw = null;
        try {
            fw = new FileWriter(file, false);
            // First row: 0..179
            for (int k = 0; k < 180; k++) {
                if (k > 0) fw.write(",");
                fw.write(Integer.toString(k));
            }
            fw.write(System.lineSeparator());

            // One row per image
            for (int i = 0; i < solRows.length; i++) {
                final double[] row = solRows[i];
                for (int k = 0; k < 180; k++) {
                    if (k > 0) fw.write(",");
                    fw.write(Double.toString(row[k]));
                }
                fw.write(System.lineSeparator());
            }
        } catch (Exception ignore) {
            // ignore
        } finally {
            if (fw != null) {
                try { fw.close(); } catch (Exception ignore) { /* ignore */ }
            }
        }
    }

    private static void writeColumnCsv(File file, String header, double[] values) {
        FileWriter fw = null;
        try {
            fw = new FileWriter(file, false);
            fw.write(header == null ? "value" : header);
            fw.write(System.lineSeparator());
            for (int i = 0; i < values.length; i++) {
                fw.write(Double.toString(values[i]));
                fw.write(System.lineSeparator());
            }
        } catch (Exception ignore) {
            // ignore
        } finally {
            if (fw != null) {
                try { fw.close(); } catch (Exception ignore) { /* ignore */ }
            }
        }
    }
}
