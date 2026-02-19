package fftanalysis.imagej;

import ij.IJ;
import ij.ImagePlus;
import ij.Macro;
import ij.io.FileInfo;
import ij.io.DirectoryChooser;
import ij.io.FileSaver;
import ij.measure.ResultsTable;
import ij.plugin.filter.PlugInFilter;
import ij.process.ByteProcessor;
import ij.process.ColorProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

import java.awt.Color;
import java.awt.GraphicsEnvironment;
import java.io.File;
import java.io.FileWriter;
import java.util.HashMap;
import java.util.Map;

/**
 * ImageJ plugin: sliding-window (overlapping boxes) orientation profiling.
 *
 * Produces a 2D (x,y) vector-field style plot of local orientations,
 * similar in spirit to MATLAB fibaall.m, by running the same MATLAB-modeled
 * core processor on overlapping square tiles.
 */
public class FIBA_Orientation_Profile implements PlugInFilter {

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
        // Batch-safe: avoid any UI prompts unless we have an ImageJ instance.
        final boolean headless = GraphicsEnvironment.isHeadless();
        final boolean interactiveUI = !headless && IJ.getInstance() != null;

        OutputOptions outOpts = new OutputOptions();
        outOpts.saveOutputs = true;
        outOpts.showVectorPlot = true;
        outOpts.showAngleMap = false;
        outOpts.overlayOnOriginal = true;
        outOpts.exportSOL = true;
        outOpts.wrap90 = true;

        TileOptions tileOpts = new TileOptions();
        // Default behavior requested: square tiles whose size equals the original image width.
        // This yields 1 tile across X and a sliding profile along Y if the image is taller than it is wide.
        // If the image is shorter than it is wide, the tile will be clamped to fit.
        tileOpts.tileIsImageWidth = true;
        tileOpts.tileSize = 0; // ignored when tileIsImageWidth=true
        tileOpts.step = 0; // 0 => auto (50% overlap)
        // Requested behavior: find image edges and place a fixed number of boxes from top->bottom.
        // Default to 6 tiles in Y.
        tileOpts.tilesY = 6;
        tileOpts.distributeYToEdges = true;
        tileOpts.vectorScale = 0.9; // relative to tileSize
        tileOpts.minStrength = 0.0;

        FibaMatlabProcessor.Params params = new FibaMatlabProcessor.Params();
        // These will be re-derived per tile based on tile size.
        params.rmin = 4;
        params.rmax = 0;
        params.alpha = 0.4;
        params.beta = 0.3;
        params.gamma = 0.3;

        String macroOpts = (argOptions != null && argOptions.trim().length() > 0)
                ? argOptions
                : Macro.getOptions();
        if (macroOpts == null || macroOpts.trim().isEmpty()) {
            final String sysOpts = System.getProperty("fiba.options");
            if (sysOpts != null && sysOpts.trim().length() > 0) macroOpts = sysOpts;
        }
        if (macroOpts != null && macroOpts.trim().length() > 0) {
            applyOptions(macroOpts, params, tileOpts, outOpts);
        }

        if (imp == null || ip == null) {
            IJ.error("No image");
            return;
        }

        // Determine output directory early (also used for debug logging)
        final String baseName = stripExtension(imp.getTitle());
        final String outputDir = determineOutputDir(imp, outOpts, interactiveUI);
        if (outOpts.debug) {
            debugToFile(outOpts, outputDir, "[FIBA_PROFILE] start interactiveUI=" + interactiveUI + " headless=" + headless);
            debugToFile(outOpts, outputDir, "[FIBA_PROFILE] opts='" + (macroOpts == null ? "" : macroOpts) + "'");
        }

        if (outOpts.saveOutputs && (outputDir == null || outputDir.trim().isEmpty())) {
            IJ.error("No output directory available (set outputDir=... in options)");
            debugToFile(outOpts, outputDir, "[FIBA_PROFILE] ERROR: outputDir is null/empty");
            return;
        }

        // Ensure grayscale for analysis
        ImageProcessor gray = ip;
        if (!(gray instanceof ByteProcessor)) {
            gray = gray.convertToByte(true);
        }

        final int width = gray.getWidth();
        final int height = gray.getHeight();

        final int tileRequested = tileOpts.tileIsImageWidth ? width : tileOpts.tileSize;
        final int tile = clampInt(tileRequested, 16, Math.min(width, height));
        final int stepAuto = Math.max(1, tile / 2);
        final int stepRequested = tileOpts.maxTilesY ? 1 : ((tileOpts.step <= 0) ? stepAuto : tileOpts.step);
        final int step = clampInt(stepRequested, 1, tile);

        final int[] xPositions = buildPositions1D(width, tile, step, 0, false);
        final int[] yPositions;
        if (tileOpts.tilesY > 0) {
            yPositions = buildPositions1D(height, tile, step, tileOpts.tilesY, tileOpts.distributeYToEdges);
        } else {
            yPositions = buildPositions1D(height, tile, step, 0, false);
        }

        final int nx = xPositions.length;
        final int ny = yPositions.length;

        if (outOpts.debug) {
            if (tileRequested != tile) {
                debugToFile(outOpts, outputDir, "[FIBA_PROFILE] NOTE: requested tile=" + tileRequested + " clampedTo=" + tile + " (image w=" + width + " h=" + height + ")");
            }
            if (tileOpts.maxTilesY) {
                debugToFile(outOpts, outputDir, "[FIBA_PROFILE] NOTE: maxTilesY=true => step forced to 1 (maximum overlap; may be slow)");
            }
            if (tileOpts.tilesY > 0) {
                debugToFile(outOpts, outputDir, "[FIBA_PROFILE] NOTE: tilesY=" + tileOpts.tilesY + " distributeYToEdges=" + tileOpts.distributeYToEdges);
            }
            debugToFile(outOpts, outputDir, "[FIBA_PROFILE] tile=" + tile + " step=" + step + " nx=" + nx + " ny=" + ny);
        }

        // Vector plot (either over original, or on a blank canvas)
        final ColorProcessor vectorCanvas;
        if (outOpts.overlayOnOriginal) {
            vectorCanvas = (ColorProcessor) ip.convertToRGB();
        } else {
            vectorCanvas = new ColorProcessor(width, height);
            vectorCanvas.setColor(Color.white);
            vectorCanvas.fill();
        }

        // Map images (grid resolution)
        final float[] angleGrid = new float[nx * ny];
        final float[] strengthGrid = new float[nx * ny];

        // MATLAB-like outputs (fibaall.m): pAng list and SOL matrix.
        // SOL layout: first row is 0..179, then one row per tile with sol[0..179].
        final double[] pAngList = new double[nx * ny];
        final double[][] solRows = new double[nx * ny][180];

        final ResultsTable rt = ResultsTable.getResultsTable();

        int idx = 0;
        for (int gy = 0; gy < ny; gy++) {
            final int y0 = yPositions[gy];
            for (int gx = 0; gx < nx; gx++) {
                final int x0 = xPositions[gx];

                gray.setRoi(x0, y0, tile, tile);
                ImageProcessor tileIp = gray.crop();

                final double[][] j = toDouble(tileIp);

                // Derive per-tile defaults.
                final int w = tile / 2;
                params.rmax = (params.rmax <= 0) ? (w - 1) : Math.min(params.rmax, w - 1);

                FibaMatlabProcessor.Result res;
                try {
                    res = FibaMatlabProcessor.process(j, params);
                } catch (Throwable t) {
                    // Keep going; mark as NaN.
                    angleGrid[idx] = Float.NaN;
                    strengthGrid[idx] = Float.NaN;
                    if (outOpts.debug) debugToFile(outOpts, outputDir, "[FIBA_PROFILE] EXCEPTION at gx=" + gx + " gy=" + gy + ": " + t);
                    idx++;
                    continue;
                }

                final float ang = (float) res.pAng;
                double pAngAdj = res.pAng;
                if (outOpts.wrap90 && pAngAdj > 90) {
                    pAngAdj = pAngAdj - 180;
                }
                final float strength = (float) res.bandStrength;
                angleGrid[idx] = (float) pAngAdj;
                strengthGrid[idx] = strength;
                pAngList[idx] = pAngAdj;
                if (res.sol != null && res.sol.length >= 180) {
                    for (int k = 0; k < 180; k++) {
                        solRows[idx][k] = res.sol[k];
                    }
                }

                // Results table row
                rt.incrementCounter();
                rt.addValue("x0", x0);
                rt.addValue("y0", y0);
                rt.addValue("tile", tile);
                rt.addValue("pAng_deg", res.pAng);
                rt.addValue("pAng_adj_deg", pAngAdj);
                rt.addValue("bandStrength", res.bandStrength);
                rt.addValue("spWid_deg", res.spWid);
                rt.addValue("pWidth_deg", res.pWidth);
                rt.addValue("warnPk", res.warnPk);

                // Draw a vector at the tile center
                if (strength >= tileOpts.minStrength) {
                    final int cx = x0 + tile / 2;
                    final int cy = y0 + tile / 2;

                    final double theta = ang * Math.PI / 180.0;
                    final double baseLen = (tile * 0.5) * tileOpts.vectorScale;
                    final double len = baseLen * clamp01(strength);

                    final double dx = Math.cos(theta) * (len / 2.0);
                    final double dy = -Math.sin(theta) * (len / 2.0); // y down in image coords

                    final int x1 = (int) Math.round(cx - dx);
                    final int y1 = (int) Math.round(cy - dy);
                    final int x2 = (int) Math.round(cx + dx);
                    final int y2 = (int) Math.round(cy + dy);

                    final Color c = Color.getHSBColor((float) (clamp01(ang / 180.0)), 1.0f, 1.0f);
                    vectorCanvas.setColor(c);
                    vectorCanvas.setLineWidth(2);
                    vectorCanvas.drawLine(x1, y1, x2, y2);
                }

                idx++;
            }
        }

        if (interactiveUI) {
            rt.show("Results");
        }

        ImagePlus vecImp = null;
        if (outOpts.showVectorPlot || outOpts.saveOutputs) {
            vecImp = new ImagePlus(baseName + "_profile_vec", vectorCanvas);
            if (outOpts.showVectorPlot) vecImp.show();
        }

        ImagePlus angMapImp = null;
        if (outOpts.showAngleMap || outOpts.saveOutputs) {
            // Simple angle map at grid resolution (ny rows x nx cols)
            FloatProcessor fp = new FloatProcessor(nx, ny, angleGrid);
            angMapImp = new ImagePlus(baseName + "_profile_ang", fp);
            // For viewing: set display range 0..180
            angMapImp.getProcessor().setMinAndMax(0, 180);
            if (outOpts.showAngleMap) angMapImp.show();
        }

        if (outOpts.saveOutputs) {
            ensureDirExists(outputDir);
            // Save vector plot as JPG
            if (vecImp != null) {
                final File outFile = new File(outputDir, baseName + "_profile_vec.jpg");
                new FileSaver(vecImp).saveAsJpeg(outFile.getAbsolutePath());
                debugToFile(outOpts, outputDir, "[FIBA_PROFILE] saved: " + outFile.getAbsolutePath());
            }

            // Save angle map as TIFF (keeps floats) + a JPG preview for convenience
            if (angMapImp != null) {
                final File outFile = new File(outputDir, baseName + "_profile_ang.tif");
                new FileSaver(angMapImp).saveAsTiff(outFile.getAbsolutePath());
                debugToFile(outOpts, outputDir, "[FIBA_PROFILE] saved: " + outFile.getAbsolutePath());

                final File outJpg = new File(outputDir, baseName + "_profile_ang.jpg");
                ImagePlus preview = new ImagePlus("preview", angMapImp.getProcessor().convertToByte(true));
                new FileSaver(preview).saveAsJpeg(outJpg.getAbsolutePath());
                debugToFile(outOpts, outputDir, "[FIBA_PROFILE] saved: " + outJpg.getAbsolutePath());
            }

            if (outOpts.exportSOL) {
                final File solFile = new File(outputDir, baseName + "_SOL.csv");
                writeSolCsv(solFile, solRows);
                debugToFile(outOpts, outputDir, "[FIBA_PROFILE] saved: " + solFile.getAbsolutePath());

                final File pangFile = new File(outputDir, baseName + "_pAng.csv");
                writeColumnCsv(pangFile, "pAng_adj_deg", pAngList);
                debugToFile(outOpts, outputDir, "[FIBA_PROFILE] saved: " + pangFile.getAbsolutePath());
            }
        }

        debugToFile(outOpts, outputDir, "[FIBA_PROFILE] done");
    }

    private static final class OutputOptions {
        boolean saveOutputs;
        boolean showVectorPlot;
        boolean showAngleMap;
        boolean overlayOnOriginal;
        boolean exportSOL;
        boolean wrap90;
        String outputDirOverride;
        boolean debug;
    }

    private static final class TileOptions {
        boolean tileIsImageWidth;
        boolean maxTilesY;
        int tilesY;
        boolean distributeYToEdges;
        int tileSize;
        int step;
        double vectorScale;
        double minStrength;
    }

    private static void applyOptions(String opts, FibaMatlabProcessor.Params params, TileOptions tile, OutputOptions out) {
        final Map<String, String> kv = parseKeyValueOptions(opts);

        // Analysis parameters
        params.rmin = parseInt(kv.get("rmin"), params.rmin);
        params.rmax = parseInt(kv.get("rmax"), params.rmax);
        params.alpha = parseDouble(kv.get("alpha"), params.alpha);
        params.beta = parseDouble(kv.get("beta"), params.beta);
        params.gamma = parseDouble(kv.get("gamma"), params.gamma);

        // Tiling / overlap parameters
        tile.tileIsImageWidth = parseBoolean(kv.get("tilewidth"), tile.tileIsImageWidth);
        // If true, use the maximum number of overlapping boxes from top to bottom (step=1).
        // Supported aliases: maxTilesY / maxOverlap / maxY
        tile.maxTilesY = parseBoolean(firstNonNull(kv.get("maxtilesy"), kv.get("maxoverlap"), kv.get("maxy")), tile.maxTilesY);
        tile.tilesY = parseInt(firstNonNull(kv.get("tilesy"), kv.get("ny"), kv.get("boxesy")), tile.tilesY);
        tile.distributeYToEdges = parseBoolean(firstNonNull(kv.get("distributey"), kv.get("edgestoy"), kv.get("edgesy")), tile.distributeYToEdges);
        tile.tileSize = parseInt(kv.get("tile"), tile.tileSize);
        tile.step = parseInt(kv.get("step"), tile.step);
        tile.vectorScale = parseDouble(kv.get("vectorscale"), tile.vectorScale);
        tile.minStrength = parseDouble(kv.get("minstrength"), tile.minStrength);

        // Outputs
        out.saveOutputs = parseBoolean(kv.get("save"), out.saveOutputs);
        out.showVectorPlot = parseBoolean(kv.get("showvec"), out.showVectorPlot);
        out.showAngleMap = parseBoolean(kv.get("showang"), out.showAngleMap);
        out.overlayOnOriginal = parseBoolean(kv.get("overlay"), out.overlayOnOriginal);
        out.exportSOL = parseBoolean(kv.get("exports"), out.exportSOL);
        out.wrap90 = parseBoolean(kv.get("wrap90"), out.wrap90);
        out.debug = parseBoolean(kv.get("debug"), out.debug);

        final String outDir = kv.get("outputdir");
        if (outDir != null && outDir.trim().length() > 0) {
            out.outputDirOverride = outDir.trim();
        }
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

    private static String firstNonNull(String a, String b, String c) {
        if (a != null) return a;
        if (b != null) return b;
        return c;
    }

    /**
     * Build top-left positions for placing tiles along one dimension.
     *
     * If fixedCount>0 and distributeToEdges=true, positions are spaced so the first tile starts at 0
     * and the last tile starts at (size - tile). This matches "find the edge" behavior.
     */
    private static int[] buildPositions1D(int size, int tile, int step, int fixedCount, boolean distributeToEdges) {
        final int maxStart = Math.max(0, size - tile);
        if (fixedCount > 0) {
            final int n = Math.max(1, fixedCount);
            final int[] pos = new int[n];
            if (n == 1) {
                pos[0] = 0;
                return pos;
            }
            if (distributeToEdges) {
                for (int i = 0; i < n; i++) {
                    final double t = (double) i / (double) (n - 1);
                    final int p = (int) Math.round(t * maxStart);
                    pos[i] = clampInt(p, 0, maxStart);
                }
            } else {
                // Fixed count but not forced to edges: just step down.
                for (int i = 0; i < n; i++) {
                    pos[i] = clampInt(i * step, 0, maxStart);
                }
            }
            return pos;
        }

        // Default stepping.
        final int n = 1 + Math.max(0, maxStart / Math.max(1, step));
        final int[] pos = new int[n];
        for (int i = 0; i < n; i++) {
            pos[i] = clampInt(i * step, 0, maxStart);
        }
        return pos;
    }

    private static int clampInt(int v, int lo, int hi) {
        if (v < lo) return lo;
        if (v > hi) return hi;
        return v;
    }

    private static double clamp01(double v) {
        if (Double.isNaN(v)) return 0;
        if (v < 0) return 0;
        if (v > 1) return 1;
        return v;
    }

    private static String determineOutputDir(ImagePlus imp, OutputOptions opts, boolean interactiveUI) {
        if (opts.outputDirOverride != null && opts.outputDirOverride.trim().length() > 0) {
            return opts.outputDirOverride;
        }

        // Requested behavior: default to the user's Downloads folder.
        final String downloads = defaultDownloadsDir();
        if (downloads != null && downloads.trim().length() > 0) {
            return downloads;
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

    private static String defaultDownloadsDir() {
        final String home = System.getProperty("user.home");
        if (home == null || home.trim().isEmpty()) return null;
        final File d = new File(home, "Downloads");
        return d.getAbsolutePath();
    }

    private static void ensureDirExists(String dir) {
        if (dir == null || dir.trim().isEmpty()) return;
        try {
            new File(dir).mkdirs();
        } catch (Exception ignore) {
            // ignore
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

    private static void debugToFile(OutputOptions out, String outputDir, String message) {
        if (out == null || !out.debug) return;
        if (outputDir == null || outputDir.trim().isEmpty()) return;
        FileWriter fw = null;
        try {
            final File f = new File(outputDir, "fiba_profile_debug.txt");
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

            // One row per tile
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
