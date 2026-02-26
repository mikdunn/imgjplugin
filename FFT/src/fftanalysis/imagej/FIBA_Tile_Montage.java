package fftanalysis.imagej;

import ij.IJ;
import ij.ImagePlus;
import ij.Macro;
import ij.gui.Plot;
import ij.io.DirectoryChooser;
import ij.io.FileInfo;
import ij.io.FileSaver;
import ij.plugin.filter.PlugInFilter;
import ij.process.ByteProcessor;
import ij.process.ColorProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import java.awt.Color;
import java.awt.Font;
import java.awt.GraphicsEnvironment;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.OutputStreamWriter;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

/**
 * Tiled FIBA workflow (no full-image processing):
 * <ul>
 *   <li>Detect content region (left/right/top) and compute a square tile size.</li>
 *   <li>Place a fixed number of square crops along Y (optionally distributed to edges).</li>
 *   <li>Run the FIBA FFT+mask reconstruction on each crop.</li>
 *   <li>Save: tile-box overlay, a montage (crop/FFT/mask/recon per tile), and an orientation-vs-tile plot.</li>
 * </ul>
 *
 * This plugin is intended for already-averaged Z-stacks (or other preprocessed images).
 * It will NOT compute any projection/averaging.
 */
public class FIBA_Tile_Montage implements PlugInFilter {

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

        Options opts = new Options();
        // IMPORTANT: keep output location stable for troubleshooting.
        // Default to the user's Downloads folder unless outputDir=... is explicitly provided.
        opts.outputDirOverride = defaultDownloadsDir();
        opts.threshold = 0;
        opts.colFrac = 0.10;
        opts.rowFrac = 0.10;
        opts.tilesY = 10;
        opts.distributeYToEdges = true;
        // Safety: never delete arbitrary files by default.
        opts.overwrite = false;
        opts.saveOverlay = true;
        opts.labelOverlay = true;
        opts.overlayLineWidth = 2;
        opts.saveMontage = true;
        opts.savePlot = true;
        opts.savePerTilePanels = true;
        opts.saveSolPlots = true;
        opts.saveCsv = true;
        opts.wrap90 = false;
        // Report/display fiber-axis direction by default: pAng_fiber = (pAng_adj + 90) mod 180.
        // This aligns reported orientation with the visible fibril direction convention.
        opts.reportFiberAxis = true;
        // Draw axes/angle labels onto the saved per-tile polar image.
        // (The montage uses a compact grayscale polar panel without labels.)
        opts.polarAxes = true;
        opts.debug = false;
        // Track specimen position in X for each tile by measuring content in that Y-band.
        // This prevents sampling the same vertical strip repeatedly when the specimen drifts/curves.
        opts.trackXPerTile = true;
        // Draw horizontal full-width section boxes (like the reference diagram).
        opts.overlayFullWidthBands = true;
        // Use grayscale jump/drop (profile gradient) to estimate left/right edges per band.
        opts.useGradientEdges = true;
        opts.edgeSmoothRadius = 6;
        // When we can estimate specimen left/right edges, keep the square slightly inside the band.
        // 1.0 = full edge-to-edge width, 0.9 = 10% margin.
        opts.tileWidthFrac = 0.90;
        // Prevent edge-based refinement from collapsing to tiny tiles.
        // (Still clamped to the image size.)
        opts.minTileSize = 128;
        // Optional explicit override (0 = auto).
        opts.tileSizeOverride = 0;

        final FibaMatlabProcessor.Params params = new FibaMatlabProcessor.Params();
        params.rmin = 4;
        params.rmax = 0;
        params.alpha = 0.4;
        params.beta = 0.3;
        params.gamma = 0.3;
        // Artifact suppression is OFF by default. We will deal with artifacts ONLY with alpha values.
        params.suppressAngleSpike = false;
        params.removeFullHeightVerticalLine = false;
        params.removeFullWidthHorizontalLine = false;

        // Options precedence: explicit arg -> Macro options -> system property.
        String macroOpts = (argOptions != null && argOptions.trim().length() > 0)
                ? argOptions
                : Macro.getOptions();
        if (macroOpts == null || macroOpts.trim().isEmpty()) {
            final String sys = System.getProperty("fiba.tilemontage.options");
            if (sys != null && sys.trim().length() > 0) macroOpts = sys;
        }
        if (macroOpts == null || macroOpts.trim().isEmpty()) {
            final String sys = System.getProperty("fiba.options");
            if (sys != null && sys.trim().length() > 0) macroOpts = sys;
        }
        if (macroOpts != null && macroOpts.trim().length() > 0) {
            applyOptions(macroOpts, opts, params);
        }

        // Ensure 8-bit grayscale for content detection and cropping.
        ImageProcessor gray = ip;
        if (!(gray instanceof ByteProcessor)) {
            gray = gray.convertToByte(true);
        }

        final String baseName = stripExtension(imp.getTitle());
        final String outputDir = determineOutputDir(imp, opts, interactiveUI);
        if (opts.saveAny() && (outputDir == null || outputDir.trim().isEmpty())) {
            IJ.error("No output directory available (set outputDir=... in options)");
            return;
        }

        final File outDir = new File(outputDir);
        if (!outDir.exists() && !outDir.mkdirs()) {
            IJ.error("Failed to create output directory: " + outDir.getAbsolutePath());
            return;
        }
        if (!outDir.isDirectory()) {
            IJ.error("Output is not a directory: " + outDir.getAbsolutePath());
            return;
        }

        IJ.log("[TILE_MONTAGE] Writing outputs to: " + outDir.getAbsolutePath());

        if (opts.overwrite) {
            // Delete ONLY files produced by this plugin run (baseName-prefixed), not the entire directory.
            clearOutputsForBaseName(outDir, baseName);
        }

        // Content detection for choosing an initial tile width and anchors.
        final int width = gray.getWidth();
        final int height = gray.getHeight();
        final Object pixObj = gray.getPixels();
        if (!(pixObj instanceof byte[])) {
            IJ.error("Expected 8-bit pixels");
            return;
        }
        final byte[] pixels = (byte[]) pixObj;
        final TileUtils.ContentBox box = TileUtils.findContentBox(pixels, width, height, opts.threshold, opts.colFrac, opts.rowFrac);

        int tileSizeGuess = TileUtils.forceEvenSize((opts.tileSizeOverride > 0) ? opts.tileSizeOverride : box.width());
        if (tileSizeGuess < 16) {
            IJ.error("Tile size too small after content detection: " + tileSizeGuess);
            return;
        }
        // Ensure the square fits in both dimensions.
        tileSizeGuess = Math.min(tileSizeGuess, Math.min(width, height));
        tileSizeGuess = TileUtils.forceEvenSize(tileSizeGuess);

        // Enforce a minimum tile size (unless image is smaller).
        if (opts.minTileSize > 0) {
            final int minEven = TileUtils.forceEvenSize(Math.min(opts.minTileSize, Math.min(width, height)));
            if (minEven >= 16) {
                tileSizeGuess = Math.max(tileSizeGuess, minEven);
            }
        }

        // Provisional Y-band placement based on the initial guess.
        int[] topsGuess = computeTopsFixedCount(box.top, height, tileSizeGuess, opts.tilesY, opts.distributeYToEdges);

        // Optional: refine tileSize using per-band gradient-detected edges.
        int tileSize = tileSizeGuess;
        if (opts.useGradientEdges && (opts.tileWidthFrac > 0) && (opts.tileWidthFrac <= 1.0) && topsGuess.length > 0) {
            final int bandH = Math.max(16,
                    Math.min(tileSizeGuess,
                            Math.max(32, (int) Math.round(height / (double) Math.max(1, opts.tilesY)))));
            final int[] widths = new int[topsGuess.length];
            int wCount = 0;
            for (int i = 0; i < topsGuess.length; i++) {
                final int[] e = TileUtils.findBandEdgesByGradient(pixels, width, height, topsGuess[i], bandH, opts.edgeSmoothRadius);
                if (e[0] >= 0 && e[1] > e[0]) {
                    final int w = (e[1] - e[0] + 1);
                    if (w >= 16) widths[wCount++] = w;
                }
            }
            if (wCount >= 3) {
                final int medW = median(widths, wCount);
                final int refined = TileUtils.forceEvenSize((int) Math.round(medW * opts.tileWidthFrac));
                if (refined >= 16) tileSize = refined;
            }
        }

        // Enforce minimum after refinement too.
        if (opts.minTileSize > 0) {
            final int minEven = TileUtils.forceEvenSize(Math.min(opts.minTileSize, Math.min(width, height)));
            if (minEven >= 16) tileSize = Math.max(tileSize, minEven);
        }
        tileSize = Math.min(tileSize, Math.min(width, height));
        tileSize = TileUtils.forceEvenSize(tileSize);

        final int[] tops = computeTopsFixedCount(box.top, height, tileSize, opts.tilesY, opts.distributeYToEdges);

        // Compute per-tile left positions.
        final int[] lefts = new int[tops.length];
        final int defaultCenterX = (box.left + box.right) / 2;
        for (int i = 0; i < tops.length; i++) {
            if (opts.trackXPerTile) {
                final int[] b = (opts.useGradientEdges)
                        ? TileUtils.findBandEdgesByGradient(pixels, width, height, tops[i], tileSize, opts.edgeSmoothRadius)
                        : TileUtils.findBandContentBounds(pixels, width, height, tops[i], tileSize, opts.threshold, opts.colFrac);
                final int cx = (b[0] >= 0) ? b[2] : defaultCenterX;
                int left = TileUtils.clamp(cx - tileSize / 2, 0, Math.max(0, width - tileSize));
                if (opts.useGradientEdges && b[0] >= 0 && b[1] > b[0]) {
                    final int minLeft = TileUtils.clamp(b[0], 0, Math.max(0, width - tileSize));
                    final int maxLeft = TileUtils.clamp(b[1] - tileSize, 0, Math.max(0, width - tileSize));
                    if (maxLeft >= minLeft) {
                        left = TileUtils.clamp(cx - tileSize / 2, minLeft, maxLeft);
                    }
                }
                lefts[i] = left;
            } else {
                lefts[i] = TileUtils.clamp(box.left, 0, Math.max(0, width - tileSize));
            }
        }

        IJ.log("[TILE_MONTAGE] tileSize=" + tileSize + " tilesY=" + tops.length + " outDir=" + outDir.getAbsolutePath());

        // Save overlay (tile boxes on original) for validation.
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
                    final int left = lefts[i];
                    if (opts.overlayFullWidthBands) {
                        overlay.drawRect(0, top, width - 1, tileSize);
                    }
                    overlay.drawRect(left, top, tileSize, tileSize);
                    if (opts.labelOverlay) {
                        overlay.drawString(Integer.toString(i + 1), left + 3, top + 16);
                    }
                }
                final ImagePlus overlayImp = new ImagePlus(baseName + "_tile_boxes", overlay);
                final File outFile = new File(outDir, baseName + "_tile_boxes.jpg");
                new FileSaver(overlayImp).saveAsJpeg(outFile.getAbsolutePath());
            } catch (Exception e) {
                IJ.log("[TILE_MONTAGE] WARNING: failed to write overlay: " + e);
            }
        }

        // Per-tile processing, plus build montage and orientation plot.
        final int nTiles = tops.length;
        final double[] pangAdj = new double[nTiles];

        final String tilePrefix = (opts.tilePrefix == null || opts.tilePrefix.trim().isEmpty())
                ? (baseName + "_tile")
                : opts.tilePrefix.trim();

        CsvWriter csv = null;
        if (opts.saveCsv) {
            try {
                final File csvFile = new File(outDir, baseName + "_tile_results.csv");
                csv = new CsvWriter(csvFile);
                csv.writeLine("tile_id,left,top,size,pAng,pAng_adj,pAng_fiber_axis");
            } catch (Exception e) {
                IJ.log("[TILE_MONTAGE] WARNING: failed to create CSV: " + e);
            }
        }

        // Montage canvas: nTiles rows x 4 columns (crop, FFT, mask, recon)
        final int gap = 5;
        final int cols = 4;
        final int montageW = cols * tileSize + (cols - 1) * gap;
        final int montageH = nTiles * tileSize + (nTiles - 1) * gap;
        final ColorProcessor montage = new ColorProcessor(montageW, montageH);
        montage.setColor(Color.white);
        montage.fill();

        for (int i = 0; i < nTiles; i++) {
            final int top = tops[i];
            final int left = lefts[i];

            // Crop tile
            ImageProcessor dup = gray.duplicate();
            dup.setRoi(left, top, tileSize, tileSize);
            final ImageProcessor crop = dup.crop();

            final double[][] j = toDouble(crop);

            // Per-tile derived defaults.
            params.requireSquare = true;
            final int w = tileSize / 2;
            params.rmax = (params.rmax <= 0) ? (w - 1) : Math.min(params.rmax, w - 1);

            final FibaMatlabProcessor.Result res;
            try {
                res = FibaMatlabProcessor.process(j, params);
            } catch (Throwable t) {
                IJ.log("[TILE_MONTAGE] ERROR tile " + (i + 1) + ": " + t);
                pangAdj[i] = Double.NaN;
                continue;
            }

            double p = res.pAng;
            if (opts.wrap90 && p > 90) p = p - 180;
            final double pFiber = normalizeAxis180(p + 90.0);
            final double pReported = opts.reportFiberAxis ? pFiber : p;
            pangAdj[i] = pReported;

            if (csv != null) {
                try {
                    csv.writeLine(
                            (i + 1) + "," + left + "," + top + "," + tileSize + "," + res.pAng + "," + p + "," + pFiber
                    );
                } catch (Exception e) {
                    IJ.log("[TILE_MONTAGE] WARNING: failed to write CSV row: " + e);
                }
            }

            if (opts.savePerTilePanels || opts.savePerTileFftTiff32) {
                try {
                    if (opts.savePerTilePanels) {
                        savePanelJpeg(outDir, tilePrefix + (i + 1) + "_crop.jpg", res.origNorm);
                        savePanelJpeg(outDir, tilePrefix + (i + 1) + "_fft.jpg", res.imgFDisp);
                    }

                    // Optional lossless output for diagnosing JPEG/8-bit artifacts in the FFT display.
                    if (opts.savePerTileFftTiff32) {
                        savePanelTiff32(outDir, tilePrefix + (i + 1) + "_fft.tif", res.imgFDisp);
                    }

                    if (res.reconMask != null) {
                        // Frequency-domain reconstruction mask (cartesian FFT plane)
                        savePanelJpeg(outDir, tilePrefix + (i + 1) + "_mask.jpg", normalize01ByMax(res.reconMask));
                    }
                    // Polar-coordinate visualization of SOL with selected peak-band highlighted (Figure-4-style)
                    if (res.sol != null) {
                        savePolarSolJpeg(outDir, tilePrefix + (i + 1) + "_polar.jpg",
                                res.sol, res.ang1, res.ang2, tileSize, opts.polarAxes);
                    }
                    savePanelJpeg(outDir, tilePrefix + (i + 1) + "_rec.jpg", res.imgR2);
                } catch (Exception e) {
                    IJ.log("[TILE_MONTAGE] WARNING: failed to save per-tile panels: " + e);
                }
            }

            if (opts.saveSolPlots && res.sol != null) {
                try {
                    final ImagePlus solPlot = buildSolPlot(tilePrefix + (i + 1) + "_sol", res.sol, res.ang1, res.ang2);
                    final File outFile = new File(outDir, tilePrefix + (i + 1) + "_sol.jpg");
                    new FileSaver(solPlot).saveAsJpeg(outFile.getAbsolutePath());
                } catch (Exception e) {
                    IJ.log("[TILE_MONTAGE] WARNING: failed to save per-tile SOL plot: " + e);
                }
            }

            // Compute montage placement
            final int y0 = i * (tileSize + gap);

            // Column 0: original tile (normalized)
            blitGray(montage, res.origNorm, 0, y0);

            // Column 1: FFT display
            blitGray(montage, res.imgFDisp, (tileSize + gap) * 1, y0);

            // Column 2: polar-coordinate SOL with selected peak-band highlighted (matches Figure 4 "c" panels)
            if (res.sol != null) {
                blitGray(montage, renderPolarSol(res.sol, res.ang1, res.ang2, tileSize), (tileSize + gap) * 2, y0);
            } else if (res.reconMask != null) {
                // Fallback: show frequency-domain mask if SOL isn't available
                blitGray(montage, normalize01ByMax(res.reconMask), (tileSize + gap) * 2, y0);
            }

            // Column 3: reconstructed mask image
            blitGray(montage, res.imgR2, (tileSize + gap) * 3, y0);

            if (opts.debug) {
                debugToFile(outDir, "[TILE_MONTAGE] tile=" + (i + 1) + " y=" + top + " pAng=" + res.pAng + " pAdj=" + p + " pFiber=" + pFiber + " pReported=" + pReported);
            }
        }

        if (opts.saveMontage) {
            final ImagePlus montageImp = new ImagePlus(baseName + "_tile_montage", montage);
            final File outFile = new File(outDir, baseName + "_tile_montage.jpg");
            new FileSaver(montageImp).saveAsJpeg(outFile.getAbsolutePath());
        }

        if (opts.savePlot) {
            final ImagePlus plotImp = buildTilePlot(baseName + "_tile_profile", pangAdj, opts.wrap90, opts.reportFiberAxis);
            final File outFile = new File(outDir, baseName + "_tile_profile.jpg");
            new FileSaver(plotImp).saveAsJpeg(outFile.getAbsolutePath());
        }

        if (csv != null) {
            try {
                csv.close();
            } catch (Exception ignore) {
                // ignore
            }
        }

        if (interactiveUI) {
            IJ.showStatus("Done: " + nTiles + " tiles");
        }
    }

    private static final class Options {
        String outputDirOverride;
        String tilePrefix;
        int threshold;
        double colFrac;
        double rowFrac;
        int tilesY;
        boolean distributeYToEdges;
        boolean trackXPerTile;
        boolean overlayFullWidthBands;
        boolean useGradientEdges;
        int edgeSmoothRadius;
        double tileWidthFrac;
        int minTileSize;
        int tileSizeOverride;
        boolean overwrite;
        boolean saveOverlay;
        boolean labelOverlay;
        int overlayLineWidth;
        boolean saveMontage;
        boolean savePlot;
        boolean savePerTilePanels;
        boolean savePerTileFftTiff32;
        boolean polarAxes;
        boolean saveSolPlots;
        boolean saveCsv;
        boolean wrap90;
        boolean reportFiberAxis;
        boolean debug;

        boolean saveAny() {
            return saveOverlay || saveMontage || savePlot || savePerTilePanels || savePerTileFftTiff32 || saveSolPlots || saveCsv;
        }
    }

    private static void applyOptions(String opts, Options out, FibaMatlabProcessor.Params params) {
        final Map<String, String> kv = parseKeyValueOptions(opts);

        final String outDir = kv.get("outputdir");
        if (outDir != null && outDir.trim().length() > 0) out.outputDirOverride = outDir.trim();

        final String tilePrefix = firstNonNull(kv.get("tileprefix"), kv.get("prefix"), kv.get("tilebasename"));
        if (tilePrefix != null && tilePrefix.trim().length() > 0) out.tilePrefix = tilePrefix.trim();

        out.threshold = parseInt(kv.get("threshold"), out.threshold);
        out.colFrac = parseDouble(kv.get("colfrac"), out.colFrac);
        out.rowFrac = parseDouble(kv.get("rowfrac"), out.rowFrac);
        out.tilesY = parseInt(firstNonNull(kv.get("tilesy"), kv.get("count"), kv.get("tiles")), out.tilesY);
        out.distributeYToEdges = parseBoolean(firstNonNull(kv.get("distributetoedges"), kv.get("distributeytoedges"), kv.get("edges")), out.distributeYToEdges);
        out.trackXPerTile = parseBoolean(firstNonNull(kv.get("trackx"), kv.get("trackxpertile"), kv.get("followblob")), out.trackXPerTile);
        out.overlayFullWidthBands = parseBoolean(firstNonNull(kv.get("overlayfullwidthbands"), kv.get("fullwidthbands"), kv.get("hsections")), out.overlayFullWidthBands);
        out.useGradientEdges = parseBoolean(firstNonNull(kv.get("usegradientedges"), kv.get("gradientedges"), kv.get("usegradient")), out.useGradientEdges);
        out.edgeSmoothRadius = parseInt(firstNonNull(kv.get("edgesmoothradius"), kv.get("smoothradius"), kv.get("edgesmooth")), out.edgeSmoothRadius);
        out.tileWidthFrac = parseDouble(firstNonNull(kv.get("tilewidthfrac"), kv.get("tilefrac"), kv.get("widthfrac")), out.tileWidthFrac);
        out.minTileSize = parseInt(firstNonNull(kv.get("mintilesize"), kv.get("mintile"), kv.get("minsize")), out.minTileSize);
        out.tileSizeOverride = parseInt(firstNonNull(kv.get("tilesize"), kv.get("tileside"), kv.get("size")), out.tileSizeOverride);

        out.overwrite = parseBoolean(kv.get("overwrite"), out.overwrite);
        out.saveOverlay = parseBoolean(firstNonNull(kv.get("saveoverlay"), kv.get("overlay"), kv.get("drawboxes")), out.saveOverlay);
        out.labelOverlay = parseBoolean(firstNonNull(kv.get("labeloverlay"), kv.get("labelboxes"), kv.get("labels")), out.labelOverlay);
        out.overlayLineWidth = parseInt(firstNonNull(kv.get("overlaylinewidth"), kv.get("boxlinewidth"), kv.get("linewidth")), out.overlayLineWidth);
        out.saveMontage = parseBoolean(firstNonNull(kv.get("savemontage"), kv.get("montage"), kv.get("savemontagejpg")), out.saveMontage);
        out.savePlot = parseBoolean(firstNonNull(kv.get("saveplot"), kv.get("plot"), kv.get("saveprofile")), out.savePlot);
        out.savePerTilePanels = parseBoolean(firstNonNull(kv.get("savepertile"), kv.get("savepanels"), kv.get("savetiles")), out.savePerTilePanels);
        out.savePerTileFftTiff32 = parseBoolean(firstNonNull(
            kv.get("savepertileffttiff32"),
            kv.get("savepertileffttiff"),
            firstNonNull(kv.get("saveffttiff32"), kv.get("saveffttiff"), kv.get("saveffttif"))
        ), out.savePerTileFftTiff32);

        out.polarAxes = parseBoolean(firstNonNull(
            kv.get("polaraxes"),
            kv.get("showpolaraxes"),
            kv.get("labelpolar")
        ), out.polarAxes);

        out.saveSolPlots = parseBoolean(firstNonNull(kv.get("savesol"), kv.get("savesolplots"), kv.get("solplots")), out.saveSolPlots);
        out.saveCsv = parseBoolean(firstNonNull(kv.get("savecsv"), kv.get("csv"), null), out.saveCsv);
        out.wrap90 = parseBoolean(kv.get("wrap90"), out.wrap90);
        out.reportFiberAxis = parseBoolean(firstNonNull(
            kv.get("reportfiberaxis"),
            kv.get("fiberaxis"),
            kv.get("add90forreport")
        ), out.reportFiberAxis);
        out.debug = parseBoolean(kv.get("debug"), out.debug);

        // FIBA params
        params.rmin = parseInt(kv.get("rmin"), params.rmin);
        params.rmax = parseInt(kv.get("rmax"), params.rmax);
        params.alpha = parseDouble(kv.get("alpha"), params.alpha);
        params.beta = parseDouble(kv.get("beta"), params.beta);
        params.gamma = parseDouble(kv.get("gamma"), params.gamma);

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

        params.removeFullWidthHorizontalLine = parseBoolean(
                firstNonNull(kv.get("removefullwidthhorizontalline"), kv.get("removehorizontalfullline"), kv.get("removehorizontalrow")),
                params.removeFullWidthHorizontalLine);
        params.removeFullWidthHorizontalLineMinCoverage = parseDouble(
                firstNonNull(kv.get("removefullwidthhorizontallinemincoverage"), kv.get("removehorizontalmincoverage"), kv.get("horizontalmincoverage")),
                params.removeFullWidthHorizontalLineMinCoverage);
    }

    private static final class CsvWriter {
        private final FileOutputStream fos;
        private final OutputStreamWriter osw;

        CsvWriter(File file) throws Exception {
            this.fos = new FileOutputStream(file);
            this.osw = new OutputStreamWriter(fos, StandardCharsets.UTF_8);
        }

        void writeLine(String line) throws Exception {
            osw.write(line);
            osw.write("\n");
            osw.flush();
        }

        void close() throws Exception {
            try {
                osw.flush();
            } finally {
                osw.close();
                fos.close();
            }
        }
    }

    private static int[] computeTopsFixedCount(int startY, int imageHeight, int tileSize, int tilesY, boolean distributeToEdges) {
        final int endTop = imageHeight - tileSize;
        startY = TileUtils.clamp(startY, 0, endTop);
        if (tilesY <= 1) {
            return new int[]{startY};
        }

        // If distributeToEdges, include endTop; otherwise distribute from startY down by step.
        final int first = startY;
        final int last = distributeToEdges ? endTop : startY + (tilesY - 1) * Math.max(1, tileSize / 2);
        final int clampedLast = TileUtils.clamp(last, 0, endTop);

        final int[] tops = new int[tilesY];
        for (int i = 0; i < tilesY; i++) {
            final double t = (tilesY == 1) ? 0.0 : (i / (double) (tilesY - 1));
            final int y = (int) Math.round(first + t * (clampedLast - first));
            tops[i] = TileUtils.clamp(y, 0, endTop);
        }
        return tops;
    }

    private static int median(int[] values, int length) {
        if (length <= 0) throw new IllegalArgumentException("length must be > 0");
        if (length > values.length) throw new IllegalArgumentException("length exceeds array");
        final int[] tmp = Arrays.copyOf(values, length);
        Arrays.sort(tmp);
        final int mid = length / 2;
        if (length % 2 == 1) return tmp[mid];
        return (tmp[mid - 1] + tmp[mid]) / 2;
    }

    private static ImagePlus buildTilePlot(String title, double[] pangAdj, boolean wrap90, boolean reportFiberAxis) {
        final int n = pangAdj.length;
        final double[] x = new double[n];
        final double[] y = new double[n];
        double sum = 0;
        int count = 0;
        for (int i = 0; i < n; i++) {
            x[i] = i + 1;
            y[i] = pangAdj[i];
            if (!Double.isNaN(y[i]) && !Double.isInfinite(y[i])) {
                sum += y[i];
                count++;
            }
        }
        final double mean = (count == 0) ? Double.NaN : (sum / count);

        final String yLabel;
        if (reportFiberAxis) {
            yLabel = wrap90 ? "Fiber Axis (deg, +90 adj/wrapped)" : "Fiber Axis (deg, +90 adj)";
        } else {
            yLabel = wrap90 ? "Fibril Orientation (deg, adj)" : "Fibril Orientation (deg)";
        }
        final Plot plot = new Plot(title, "Crop Square ID", yLabel, x, y);
        plot.setLineWidth(2);
        plot.setColor(Color.black);

        if (!Double.isNaN(mean)) {
            plot.setColor(Color.blue);
            final double[] my = new double[n];
            for (int i = 0; i < n; i++) my[i] = mean;
            plot.add("line", x, my);
        }

        return plot.getImagePlus();
    }

    private static ImagePlus buildSolPlot(String title, double[] sol180, int ang1, int ang2) {
        final int n = 180;
        final double[] x = new double[n];
        final double[] y = new double[n];
        final double[] yBand = new double[n];
        final boolean[] inBand = computeInBand(ang1, ang2);

        double max = 0;
        for (int i = 0; i < n; i++) {
            x[i] = i;
            y[i] = sol180[i];
            if (y[i] > max) max = y[i];
        }
        if (!(max > 0)) max = 1;

        for (int i = 0; i < n; i++) {
            // Normalize for visibility
            y[i] = y[i] / max;
            yBand[i] = inBand[i] ? y[i] : Double.NaN;
        }

        final Plot plot = new Plot(title, "Angle (deg)", "SOL (norm)", x, y);
        plot.setLineWidth(1);
        plot.setColor(Color.black);
        plot.add("line", x, y);

        plot.setColor(Color.red);
        plot.setLineWidth(4);
        plot.add("line", x, yBand);

        return plot.getImagePlus();
    }

    private static boolean[] computeInBand(int ang1, int ang2) {
        final boolean[] in = new boolean[180];
        final int a1 = mod180(ang1);
        final int a2 = mod180(ang2);
        if (a2 >= a1) {
            for (int a = a1; a <= a2; a++) in[a] = true;
        } else {
            for (int a = a1; a < 180; a++) in[a] = true;
            for (int a = 0; a <= a2; a++) in[a] = true;
        }
        return in;
    }

    /**
     * Render a simple polar plot of SOL in a square grayscale image (0..1) and highlight the selected peak-band.
     *
     * Angle convention matches the processor: theta=0 follows +row (down), theta=90 follows +col (right).
     */
    private static double[][] renderPolarSol(double[] sol180, int ang1, int ang2, int size) {
        if (sol180 == null || sol180.length != 180 || size <= 0) return null;

        double max = 0;
        for (int i = 0; i < 180; i++) if (sol180[i] > max) max = sol180[i];
        if (!(max > 0)) max = 1;

        final int a1 = mod180(ang1);
        final int a2 = mod180(ang2);
        final boolean[] inBand = computeInBand(ang1, ang2);

        final double[][] out = new double[size][size];
        final double cx = (size - 1) / 2.0;
        final double cy = (size - 1) / 2.0;
        final double R = (size / 2.0) - 2.0;

        // Paint a visible wedge background for the selected peak band so highlighting is obvious.
        for (int deg = 0; deg < 360; deg++) {
            final int a = deg % 180;
            if (!isInBandAngle(a, a1, a2)) continue;
            final double theta = Math.toRadians(deg);
            final double dirRow = Math.cos(theta);
            final double dirCol = Math.sin(theta);
            for (int r = 0; r <= (int) Math.round(R); r++) {
                final int yy = (int) Math.round(cy + dirRow * r);
                final int xx = (int) Math.round(cx + dirCol * r);
                if (yy < 0 || yy >= size || xx < 0 || xx >= size) continue;
                if (out[yy][xx] < 0.28) out[yy][xx] = 0.28;
            }
        }

        // Draw SOL rays (0..359) using mirrored SOL (SOL is 180-periodic)
        for (int deg = 0; deg < 360; deg++) {
            final int a = deg % 180;
            final double v = sol180[a] / max;
            final double theta = Math.toRadians(deg);

            // Image coordinates: row increases downward; theta=0 should point down.
            final double dirRow = Math.cos(theta);
            final double dirCol = Math.sin(theta);

            final double rLen = v * R;
            final double intensity = inBand[a] ? 1.0 : 0.10;

            for (int r = 0; r <= (int) Math.round(rLen); r++) {
                final int yy = (int) Math.round(cy + dirRow * r);
                final int xx = (int) Math.round(cx + dirCol * r);
                if (yy < 0 || yy >= size || xx < 0 || xx >= size) continue;
                if (out[yy][xx] < intensity) out[yy][xx] = intensity;
            }
        }

        // Draw a circle outline
        for (int deg = 0; deg < 360; deg++) {
            final double theta = Math.toRadians(deg);
            final int yy = (int) Math.round(cy + Math.cos(theta) * R);
            final int xx = (int) Math.round(cx + Math.sin(theta) * R);
            if (yy < 0 || yy >= size || xx < 0 || xx >= size) continue;
            out[yy][xx] = 1.0;
        }

        // Draw band boundary rays for ang1 and ang2 (and +180)
        drawPolarRay(out, cx, cy, R, a1, 1.0);
        drawPolarRay(out, cx, cy, R, a2, 1.0);
        drawPolarRay(out, cx, cy, R, a1 + 180, 1.0);
        drawPolarRay(out, cx, cy, R, a2 + 180, 1.0);

        return out;
    }

    private static void drawPolarRay(double[][] img, double cx, double cy, double R, int deg, double intensity) {
        final int size = img.length;
        final double theta = Math.toRadians(deg);
        final double dirRow = Math.cos(theta);
        final double dirCol = Math.sin(theta);
        for (int r = 0; r <= (int) Math.round(R); r++) {
            final int yy = (int) Math.round(cy + dirRow * r);
            final int xx = (int) Math.round(cx + dirCol * r);
            if (yy < 0 || yy >= size || xx < 0 || xx >= size) continue;
            if (img[yy][xx] < intensity) img[yy][xx] = intensity;
        }
    }

    private static boolean isInBandAngle(int angle0to179, int a1, int a2) {
        final int a = mod180(angle0to179);
        if (a2 >= a1) return a >= a1 && a <= a2;
        return a >= a1 || a <= a2;
    }

    private static int mod180(int a) {
        int m = a % 180;
        if (m < 0) m += 180;
        return m;
    }

    private static double normalizeAxis180(double a) {
        double m = a % 180.0;
        if (m < 0) m += 180.0;
        return m;
    }

    private static String determineOutputDir(ImagePlus imp, Options opts, boolean interactiveUI) {
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

    /**
     * Best-effort default output directory.
     *
     * For reproducible troubleshooting, prefer a stable location (Downloads) over the image's source directory.
     */
    private static String defaultDownloadsDir() {
        final String home = System.getProperty("user.home");
        if (home == null || home.trim().isEmpty()) return null;
        final File downloads = new File(home, "Downloads");
        if (downloads.exists() && downloads.isDirectory()) {
            return downloads.getAbsolutePath();
        }
        // Fallback: use home directory if Downloads doesn't exist.
        final File h = new File(home);
        if (h.exists() && h.isDirectory()) {
            return h.getAbsolutePath();
        }
        return null;
    }

    private static void clearOutputsForBaseName(File dir, String baseName) {
        if (dir == null || baseName == null || baseName.trim().isEmpty()) return;
        final File[] files = dir.listFiles();
        if (files == null) return;
        final String prefix = baseName + "_";
        for (final File f : files) {
            if (!f.isFile()) continue;
            final String name = f.getName();
            if (name == null) continue;
            if (!name.startsWith(prefix)) continue;
            try {
                //noinspection ResultOfMethodCallIgnored
                f.delete();
            } catch (Exception ignore) {
                // ignore
            }
        }
    }

    private static void debugToFile(File outDir, String message) {
        if (outDir == null) return;
        FileWriter fw = null;
        try {
            final File f = new File(outDir, "fiba_tile_montage_debug.txt");
            fw = new FileWriter(f, true);
            fw.write(message);
            fw.write(System.lineSeparator());
        } catch (java.io.IOException ignore) {
            // ignore
        } finally {
            if (fw != null) {
                try { fw.close(); } catch (java.io.IOException ignore) { /* ignore */ }
            }
        }
    }

    private static void blitGray(ColorProcessor cp, double[][] img01, int x0, int y0) {
        if (img01 == null || img01.length == 0 || img01[0].length == 0) return;
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

    private static void savePanelJpeg(File outDir, String fileName, double[][] img01) {
        if (outDir == null || fileName == null || img01 == null) return;
        final int h = img01.length;
        final int w = img01[0].length;
        final byte[] pix = new byte[w * h];
        int k = 0;
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                pix[k++] = (byte) (int) Math.round(clamp01(img01[y][x]) * 255.0);
            }
        }
        final ByteProcessor bp = new ByteProcessor(w, h, pix);
        final ImagePlus imp = new ImagePlus(fileName, bp);
        final File outFile = new File(outDir, fileName);
        new FileSaver(imp).saveAsJpeg(outFile.getAbsolutePath());
    }

    private static void savePolarSolJpeg(
            File outDir,
            String fileName,
            double[] sol180,
            int ang1,
            int ang2,
            int size,
            boolean drawAxes
    ) {
        if (outDir == null || fileName == null) return;
        final double[][] img01 = renderPolarSol(sol180, ang1, ang2, size);
        if (img01 == null) return;

        final int h = img01.length;
        final int w = img01[0].length;
        final byte[] pix = new byte[w * h];
        int k = 0;
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                pix[k++] = (byte) (int) Math.round(clamp01(img01[y][x]) * 255.0);
            }
        }
        final ByteProcessor bp = new ByteProcessor(w, h, pix);
        ImageProcessor ip = bp;
        if (drawAxes) {
            final ColorProcessor cp = (ColorProcessor) bp.convertToRGB();
            drawPolarAxes(cp);
            ip = cp;
        }
        final ImagePlus imp = new ImagePlus(fileName, ip);
        final File outFile = new File(outDir, fileName);
        new FileSaver(imp).saveAsJpeg(outFile.getAbsolutePath());
    }

    /**
     * Draw simple polar axes and angle labels on a polar SOL image.
     * Angle convention matches the processor: theta=0 follows +row (down), theta=90 follows +col (right).
     */
    private static void drawPolarAxes(ColorProcessor cp) {
        if (cp == null) return;
        final int w = cp.getWidth();
        final int h = cp.getHeight();
        final int cx = (w - 1) / 2;
        final int cy = (h - 1) / 2;

        cp.setLineWidth(1);
        // Light gray so it doesn't overpower the SOL rays.
        cp.setColor(new Color(210, 210, 210));
        cp.drawLine(cx, 0, cx, h - 1);
        cp.drawLine(0, cy, w - 1, cy);

        // Tick marks every 45 degrees, near the perimeter.
        final double R = (Math.min(w, h) / 2.0) - 2.0;
        for (int deg = 0; deg < 360; deg += 45) {
            final double theta = Math.toRadians(deg);
            final double dirRow = Math.cos(theta);
            final double dirCol = Math.sin(theta);
            final int x1 = (int) Math.round(cx + dirCol * (R - 6));
            final int y1 = (int) Math.round(cy + dirRow * (R - 6));
            final int x2 = (int) Math.round(cx + dirCol * (R));
            final int y2 = (int) Math.round(cy + dirRow * (R));
            cp.drawLine(x1, y1, x2, y2);
        }

        cp.setFont(new Font("SansSerif", Font.PLAIN, 12));
        cp.setColor(Color.white);
        // Labels: 0=down, 90=right, 180=up, 270=left.
        cp.drawString("0", cx + 4, h - 6);
        cp.drawString("90", w - 24, cy - 4);
        cp.drawString("180", cx - 10, 14);
        cp.drawString("270", 4, cy - 4);
    }

    /**
     * Save a panel as a lossless 32-bit TIFF.
     *
     * This is intended for debugging: it preserves exact pixel values (no JPEG artifacts, no 8-bit quantization).
     */
    private static void savePanelTiff32(File outDir, String fileName, double[][] img) {
        if (outDir == null || fileName == null || img == null) return;
        final int h = img.length;
        final int w = img[0].length;
        final float[] pix = new float[w * h];
        int k = 0;
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                // Do not clamp here: if something is outside [0,1], we want to preserve that for diagnosis.
                pix[k++] = (float) img[y][x];
            }
        }
        final FloatProcessor fp = new FloatProcessor(w, h, pix);
        final ImagePlus imp = new ImagePlus(fileName, fp);
        final File outFile = new File(outDir, fileName);
        new FileSaver(imp).saveAsTiff(outFile.getAbsolutePath());
    }

    private static double[][] normalize01ByMax(double[][] in) {
        if (in == null || in.length == 0 || in[0].length == 0) return in;
        final int h = in.length;
        final int w = in[0].length;
        double max = 0;
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                final double v = in[y][x];
                if (v > max) max = v;
            }
        }
        if (!(max > 0)) return in;
        final double inv = 1.0 / max;
        final double[][] out = new double[h][w];
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                out[y][x] = in[y][x] * inv;
            }
        }
        return out;
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

    private static double clamp01(double v) {
        if (v < 0) return 0;
        if (v > 1) return 1;
        return v;
    }

    private static String stripExtension(String name) {
        if (name == null) return "image";
        int dot = name.lastIndexOf('.');
        if (dot > 0) return name.substring(0, dot);
        return name;
    }

    private static String firstNonNull(String a, String b, String c) {
        if (a != null) return a;
        if (b != null) return b;
        return c;
    }

    private static Map<String, String> parseKeyValueOptions(String opts) {
        final Map<String, String> map = new HashMap<>();
        if (opts == null) return map;
        final String[] tokens = opts.trim().split("\\s+");
        for (final String t : tokens) {
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
        } catch (NumberFormatException e) {
            return def;
        }
    }

    private static double parseDouble(String s, double def) {
        if (s == null) return def;
        try {
            return Double.parseDouble(s);
        } catch (NumberFormatException e) {
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
}
