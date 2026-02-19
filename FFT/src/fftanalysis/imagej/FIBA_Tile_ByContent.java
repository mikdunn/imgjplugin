package fftanalysis.imagej;

import ij.IJ;
import ij.ImagePlus;
import ij.Macro;
import ij.io.FileSaver;
import ij.plugin.filter.PlugInFilter;
import ij.process.ByteProcessor;
import ij.process.ColorProcessor;
import ij.process.ImageProcessor;

import java.awt.Color;
import java.awt.Font;
import java.awt.GraphicsEnvironment;
import java.io.File;
import java.io.FileWriter;
import java.util.HashMap;
import java.util.Map;

/**
 * Automatically tiles the current image into vertically-overlapping square crops.
 *
 * Logic:
 *  - find the "content" region based on non-black pixels
 *  - set the tile width to the content width (left/right bounds), forced to even
 *  - generate the minimal number of overlapping tiles that reach the bottom
 *  - write numbered tiles (1.jpg, 2.jpg, ...) to an output folder
 */
public class FIBA_Tile_ByContent implements PlugInFilter {

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
        if (imp == null || ip == null) {
            IJ.error("No image");
            return;
        }

        final boolean headless = GraphicsEnvironment.isHeadless();
        final boolean interactiveUI = !headless && IJ.getInstance() != null;

        // Ensure 8-bit grayscale.
        ImageProcessor gray = ip;
        if (!(gray instanceof ByteProcessor)) {
            gray = gray.convertToByte(true);
        }

        Options opts = new Options();
        opts.threshold = 0;         // treat >0 as non-black by default
        opts.colFrac = 0.10;        // columns with >=10% of max non-black count define content
        opts.rowFrac = 0.10;        // rows with >=10% of max non-black count define content
        opts.overwrite = true;
        opts.debug = false;
        // Visual debug/validation output: draw tile boxes on the original image.
        opts.saveOverlay = true;
        opts.labelOverlay = true;
        opts.overlayLineWidth = 2;

        // Options precedence: explicit arg -> Macro options -> system property.
        String macroOpts = (argOptions != null && argOptions.trim().length() > 0)
                ? argOptions
                : Macro.getOptions();
        if (macroOpts == null || macroOpts.trim().isEmpty()) {
            String sys = System.getProperty("fiba.tile.options");
            if (sys != null && sys.trim().length() > 0) macroOpts = sys;
        }
        if (macroOpts == null || macroOpts.trim().isEmpty()) {
            // Back-compat convenience: allow reusing -Dfiba.options=...
            String sys = System.getProperty("fiba.options");
            if (sys != null && sys.trim().length() > 0) macroOpts = sys;
        }

        if (macroOpts != null && macroOpts.trim().length() > 0) {
            applyOptions(macroOpts, opts);
        }

        final String baseName = stripExtension(imp.getTitle());

        final String outputDir = (opts.outputDirOverride != null && opts.outputDirOverride.trim().length() > 0)
                ? opts.outputDirOverride.trim()
                : defaultTilesDir(baseName);

        final File outDir = new File(outputDir);
        if (!outDir.exists()) {
            if (!outDir.mkdirs()) {
                IJ.error("Failed to create output directory: " + outDir.getAbsolutePath());
                return;
            }
        }
        if (!outDir.isDirectory()) {
            IJ.error("Output is not a directory: " + outDir.getAbsolutePath());
            return;
        }

        if (opts.overwrite) {
            clearDirectoryFiles(outDir);
        }

        final int width = gray.getWidth();
        final int height = gray.getHeight();

        final Object pixObj = gray.getPixels();
        if (!(pixObj instanceof byte[])) {
            IJ.error("Expected 8-bit pixels");
            return;
        }
        final byte[] pixels = (byte[]) pixObj;

        TileUtils.ContentBox box = TileUtils.findContentBox(pixels, width, height, opts.threshold, opts.colFrac, opts.rowFrac);

        int tileSize = TileUtils.forceEvenSize(box.width());
        if (tileSize < 8) {
            IJ.error("Tile size too small after content detection: " + tileSize);
            return;
        }

        // Ensure the square fits vertically.
        if (tileSize > height) {
            tileSize = TileUtils.forceEvenSize(height);
        }
        if (tileSize <= 0) {
            IJ.error("Invalid tile size");
            return;
        }

        int left = box.left;
        // Clamp left so [left, left+tileSize) stays in image.
        left = TileUtils.clamp(left, 0, Math.max(0, width - tileSize));

        int[] tops;
        try {
            tops = TileUtils.computeTileTops(box.top, height, tileSize);
        } catch (IllegalArgumentException e) {
            IJ.error("Cannot tile image: " + e.getMessage());
            return;
        }

        IJ.log("[TILE] content left=" + box.left + " right=" + box.right + " top=" + box.top);
        IJ.log("[TILE] tileSize=" + tileSize + " count=" + tops.length + " outputDir=" + outDir.getAbsolutePath());

        // Persist a small manifest so subsequent steps (e.g., FIBA_All_FromFolder) can
        // deterministically know how many tiles were produced and where they came from.
        writeTileManifests(outDir, box, left, tileSize, width, height, opts, tops);

        // Save an overlay image with the tiling boxes drawn, for visual verification.
        if (opts.saveOverlay) {
            try {
                final ColorProcessor overlay = (ColorProcessor) ip.convertToRGB();
                overlay.setColor(Color.red);
                overlay.setLineWidth(Math.max(1, opts.overlayLineWidth));
                if (opts.labelOverlay) {
                    overlay.setFont(new Font("SansSerif", Font.BOLD, 14));
                }

                for (int i = 0; i < tops.length; i++) {
                    final int top = tops[i];
                    overlay.drawRect(left, top, tileSize, tileSize);
                    if (opts.labelOverlay) {
                        overlay.drawString(Integer.toString(i + 1), left + 3, top + 16);
                    }
                }

                final ImagePlus overlayImp = new ImagePlus(baseName + "_tile_boxes", overlay);
                final File outFile = new File(outDir, baseName + "_tile_boxes.jpg");
                new FileSaver(overlayImp).saveAsJpeg(outFile.getAbsolutePath());
            } catch (Exception e) {
                IJ.log("[TILE] WARNING: failed to write tile box overlay: " + e);
            }
        }

        // Save tiles in top-to-bottom order as 1.jpg, 2.jpg, ...
        for (int i = 0; i < tops.length; i++) {
            int top = tops[i];
            ImageProcessor dup = gray.duplicate();
            dup.setRoi(left, top, tileSize, tileSize);
            ImageProcessor crop = dup.crop();
            ImagePlus tile = new ImagePlus(Integer.toString(i + 1), crop);

            File outFile = new File(outDir, (i + 1) + ".jpg");
            boolean ok = new FileSaver(tile).saveAsJpeg(outFile.getAbsolutePath());
            if (!ok) {
                IJ.log("[TILE] Failed to write: " + outFile.getAbsolutePath());
            }
        }

        if (interactiveUI) {
            IJ.showStatus("Tiled " + tops.length + " crops to " + outDir.getAbsolutePath());
        }
    }

    private static final class Options {
        String outputDirOverride;
        int threshold;
        double colFrac;
        double rowFrac;
        boolean overwrite;
        boolean debug;
        boolean saveOverlay;
        boolean labelOverlay;
        int overlayLineWidth;
    }

    private static void applyOptions(String opts, Options out) {
        final Map<String, String> kv = parseKeyValueOptions(opts);

        final String outDir = kv.get("outputdir");
        if (outDir != null && outDir.trim().length() > 0) out.outputDirOverride = outDir.trim();

        out.threshold = parseInt(kv.get("threshold"), out.threshold);
        out.colFrac = parseDouble(kv.get("colfrac"), out.colFrac);
        out.rowFrac = parseDouble(kv.get("rowfrac"), out.rowFrac);

        out.overwrite = parseBoolean(kv.get("overwrite"), out.overwrite);
        out.debug = parseBoolean(kv.get("debug"), out.debug);

        // Visual validation outputs
        out.saveOverlay = parseBoolean(firstNonNull(kv.get("saveoverlay"), kv.get("overlay"), kv.get("drawboxes")), out.saveOverlay);
        out.labelOverlay = parseBoolean(firstNonNull(kv.get("labeloverlay"), kv.get("labelboxes"), kv.get("labels")), out.labelOverlay);
        out.overlayLineWidth = parseInt(firstNonNull(kv.get("overlaylinewidth"), kv.get("boxlinewidth"), kv.get("linewidth")), out.overlayLineWidth);
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

    private static String defaultTilesDir(String baseName) {
        // Windows-friendly default: %USERPROFILE%\Downloads\<baseName>_tiles
        String home = System.getProperty("user.home");
        if (home == null || home.trim().isEmpty()) {
            // Fallback: current directory
            return new File(baseName + "_tiles").getAbsolutePath();
        }
        File downloads = new File(home, "Downloads");
        return new File(downloads, baseName + "_tiles").getAbsolutePath();
    }

    private static void clearDirectoryFiles(File dir) {
        File[] files = dir.listFiles();
        if (files == null) return;
        for (int i = 0; i < files.length; i++) {
            File f = files[i];
            if (f.isFile()) {
                // best effort
                try {
                    //noinspection ResultOfMethodCallIgnored
                    f.delete();
                } catch (Exception ignore) {
                    // ignore
                }
            }
        }
    }

    private static void writeTileManifests(
            File outDir,
            TileUtils.ContentBox box,
            int usedLeft,
            int tileSize,
            int imageWidth,
            int imageHeight,
            Options opts,
            int[] tops
    ) {
        // Best-effort; never fail the tiling operation due to manifest issues.
        FileWriter fw = null;
        try {
            File info = new File(outDir, "tiles_info.txt");
            fw = new FileWriter(info, false);
            fw.write("content_left=" + box.left + System.lineSeparator());
            fw.write("content_right=" + box.right + System.lineSeparator());
            fw.write("content_top=" + box.top + System.lineSeparator());
            fw.write("used_left=" + usedLeft + System.lineSeparator());
            fw.write("tile_size=" + tileSize + System.lineSeparator());
            fw.write("tile_count=" + (tops == null ? 0 : tops.length) + System.lineSeparator());
            fw.write("image_width=" + imageWidth + System.lineSeparator());
            fw.write("image_height=" + imageHeight + System.lineSeparator());
            fw.write("threshold=" + opts.threshold + System.lineSeparator());
            fw.write("colFrac=" + opts.colFrac + System.lineSeparator());
            fw.write("rowFrac=" + opts.rowFrac + System.lineSeparator());
        } catch (Exception ignore) {
            // ignore
        } finally {
            if (fw != null) {
                try { fw.close(); } catch (Exception ignore) { /* ignore */ }
            }
        }

        try {
            File csv = new File(outDir, "tiles.csv");
            fw = new FileWriter(csv, false);
            fw.write("index,left,top,size" + System.lineSeparator());
            for (int i = 0; i < tops.length; i++) {
                fw.write((i + 1) + "," + usedLeft + "," + tops[i] + "," + tileSize + System.lineSeparator());
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
