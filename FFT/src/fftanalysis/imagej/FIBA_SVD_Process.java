package fftanalysis.imagej;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.Plot;
import ij.io.FileInfo;
import ij.io.FileSaver;
import ij.io.OpenDialog;
import ij.plugin.filter.PlugInFilter;
import ij.process.ByteProcessor;
import ij.process.ColorProcessor;
import ij.process.ImageProcessor;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;

import java.awt.Color;
import java.awt.GraphicsEnvironment;
import java.io.File;
import java.text.SimpleDateFormat;
import java.util.Arrays;
import java.util.Date;

/**
 * SVD denoising/reconstruction plugin for ImageJ.
 *
 * Produces rank-20 and energy-99% reconstructions, shows/saves panel + singular-value plot,
 * and writes to timestamped run folders beside input image.
 */
public class FIBA_SVD_Process implements PlugInFilter {

    private ImagePlus imp;

    @Override
    public int setup(String arg, ImagePlus imp) {
        this.imp = imp;
        return DOES_ALL + NO_IMAGE_REQUIRED;
    }

    @Override
    public void run(ImageProcessor ip) {
        if (imp == null || ip == null) {
            final ImagePlus opened = promptForInputImage();
            if (opened == null) {
                IJ.error("No image selected");
                return;
            }
            imp = opened;
            ip = imp.getProcessor();
        }

        final boolean interactiveUI = !GraphicsEnvironment.isHeadless() && IJ.getInstance() != null;

        ImageProcessor gray = ip;
        if (!(gray instanceof ByteProcessor)) {
            gray = gray.convertToByte(true);
        }

        final String baseName = stripExtension(imp.getTitle());
        final String outputDir = determineOutputDir(imp, baseName, interactiveUI);
        if (outputDir == null || outputDir.trim().isEmpty()) {
            IJ.error("No output directory available");
            return;
        }

        final File outDir = new File(outputDir);
        if (!outDir.exists() && !outDir.mkdirs()) {
            IJ.error("Failed to create output directory: " + outDir.getAbsolutePath());
            return;
        }

        final double[][] img = toDouble(gray);
        final SVDComputeResult svd = computeSvdVariants(img, 20, 0.99);

        final ImagePlus rankImp = asGrayImage(baseName + "_svd_rank20", svd.rankRecon01);
        final ImagePlus energyImp = asGrayImage(baseName + "_svd_energy99", svd.energyRecon01);
        final ImagePlus spectrumImp = buildSpectrumPlot(baseName + "_svd_singular_values", svd.singularValues);
        final ImagePlus panelImp = buildPanel(baseName + "_svd_panel", gray, svd.rankRecon01, svd.energyRecon01);

        new FileSaver(rankImp).saveAsJpeg(new File(outDir, baseName + "_svd_rank20.jpg").getAbsolutePath());
        new FileSaver(energyImp).saveAsJpeg(new File(outDir, baseName + "_svd_energy99.jpg").getAbsolutePath());
        new FileSaver(spectrumImp).saveAsJpeg(new File(outDir, baseName + "_svd_singular_values.jpg").getAbsolutePath());
        new FileSaver(panelImp).saveAsJpeg(new File(outDir, baseName + "_svd_panel.jpg").getAbsolutePath());

        IJ.log("[SVD] Writing outputs to: " + outDir.getAbsolutePath());
        IJ.log("[SVD] rank20 energy=" + svd.rank20Energy + " energy99 rank=" + svd.energyRank);

        if (interactiveUI) {
            rankImp.show();
            energyImp.show();
            spectrumImp.show();
            panelImp.show();
            IJ.showStatus("SVD process done");
        }
    }

    private static final class SVDComputeResult {
        final double[][] rankRecon01;
        final double[][] energyRecon01;
        final double[] singularValues;
        final int energyRank;
        final double rank20Energy;

        SVDComputeResult(double[][] rankRecon01, double[][] energyRecon01, double[] singularValues, int energyRank, double rank20Energy) {
            this.rankRecon01 = rankRecon01;
            this.energyRecon01 = energyRecon01;
            this.singularValues = singularValues;
            this.energyRank = energyRank;
            this.rank20Energy = rank20Energy;
        }
    }

    private static SVDComputeResult computeSvdVariants(double[][] img, int fixedRank, double energyCutoff) {
        final RealMatrix X = new Array2DRowRealMatrix(img, true);
        final SingularValueDecomposition svd = new SingularValueDecomposition(X);
        final double[] s = svd.getSingularValues();

        final int rank20 = Math.max(1, Math.min(fixedRank, s.length));
        final int energyRank = findRankForEnergy(s, energyCutoff);

        final RealMatrix rankRecon = reconstruct(svd, rank20);
        final RealMatrix energyRecon = reconstruct(svd, energyRank);

        final double rank20Energy = cumulativeEnergyAtRank(s, rank20);

        final double[][] r20 = rankRecon.getData();
        final double[][] e99 = energyRecon.getData();
        normalizeTo01InPlace(r20);
        normalizeTo01InPlace(e99);

        return new SVDComputeResult(r20, e99, Arrays.copyOf(s, s.length), energyRank, rank20Energy);
    }

    private static RealMatrix reconstruct(SingularValueDecomposition svd, int rank) {
        final RealMatrix U = svd.getU();
        final RealMatrix V = svd.getV();
        final double[] s = svd.getSingularValues();

        final int r = Math.max(1, Math.min(rank, s.length));
        final RealMatrix Usub = U.getSubMatrix(0, U.getRowDimension() - 1, 0, r - 1);
        final RealMatrix Vsub = V.getSubMatrix(0, V.getRowDimension() - 1, 0, r - 1);
        final double[] sr = Arrays.copyOf(s, r);
        final RealMatrix Sdiag = MatrixUtils.createRealDiagonalMatrix(sr);

        return Usub.multiply(Sdiag).multiply(Vsub.transpose());
    }

    private static int findRankForEnergy(double[] s, double cutoff) {
        if (s == null || s.length == 0) return 1;
        double total = 0;
        for (int i = 0; i < s.length; i++) total += s[i] * s[i];
        if (!(total > 0)) return 1;
        double c = 0;
        for (int i = 0; i < s.length; i++) {
            c += s[i] * s[i];
            if (c / total >= cutoff) return i + 1;
        }
        return s.length;
    }

    private static double cumulativeEnergyAtRank(double[] s, int rank) {
        if (s == null || s.length == 0) return 1.0;
        double total = 0;
        for (int i = 0; i < s.length; i++) total += s[i] * s[i];
        if (!(total > 0)) return 1.0;
        final int r = Math.max(1, Math.min(rank, s.length));
        double c = 0;
        for (int i = 0; i < r; i++) c += s[i] * s[i];
        return c / total;
    }

    private static ImagePlus asGrayImage(String title, double[][] img01) {
        final int h = img01.length;
        final int w = img01[0].length;
        final byte[] pix = new byte[w * h];
        int k = 0;
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                pix[k++] = (byte) (int) Math.round(clamp01(img01[y][x]) * 255.0);
            }
        }
        return new ImagePlus(title, new ByteProcessor(w, h, pix));
    }

    private static ImagePlus buildPanel(String title, ImageProcessor gray, double[][] rank20, double[][] energy99) {
        final int w = gray.getWidth();
        final int h = gray.getHeight();
        final ColorProcessor panel = new ColorProcessor(w * 3, h);

        final ImagePlus origImp = new ImagePlus("orig", gray.duplicate().convertToByte(true));
        final ImagePlus rImp = asGrayImage("rank20", rank20);
        final ImagePlus eImp = asGrayImage("energy99", energy99);

        blit(panel, origImp.getProcessor().convertToRGB(), 0, 0);
        blit(panel, rImp.getProcessor().convertToRGB(), w, 0);
        blit(panel, eImp.getProcessor().convertToRGB(), 2 * w, 0);

        return new ImagePlus(title, panel);
    }

    private static void blit(ColorProcessor dst, ImageProcessor src, int x0, int y0) {
        final int w = src.getWidth();
        final int h = src.getHeight();
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                dst.set(x0 + x, y0 + y, src.getPixel(x, y));
            }
        }
    }

    private static ImagePlus buildSpectrumPlot(String title, double[] s) {
        final int n = Math.min(120, s.length);
        final double[] x = new double[n];
        final double[] y = new double[n];
        final double max = (s.length == 0) ? 1.0 : s[0];
        for (int i = 0; i < n; i++) {
            x[i] = i + 1;
            y[i] = (max > 0) ? (s[i] / max) : 0;
        }
        final Plot p = new Plot(title, "Index", "Normalized singular value", x, y);
        p.setColor(Color.magenta);
        p.setLineWidth(2);
        p.add("line", x, y);
        return p.getImagePlus();
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

    private static void normalizeTo01InPlace(double[][] arr) {
        double min = Double.POSITIVE_INFINITY;
        double max = Double.NEGATIVE_INFINITY;
        for (int y = 0; y < arr.length; y++) {
            for (int x = 0; x < arr[0].length; x++) {
                final double v = arr[y][x];
                if (v < min) min = v;
                if (v > max) max = v;
            }
        }
        final double range = max - min;
        if (!(range > 0)) {
            for (int y = 0; y < arr.length; y++) {
                Arrays.fill(arr[y], 0);
            }
            return;
        }
        for (int y = 0; y < arr.length; y++) {
            for (int x = 0; x < arr[0].length; x++) {
                arr[y][x] = (arr[y][x] - min) / range;
            }
        }
    }

    private static double clamp01(double v) {
        if (v < 0) return 0;
        if (v > 1) return 1;
        return v;
    }

    private static ImagePlus promptForInputImage() {
        final OpenDialog od = new OpenDialog("Select image to analyze (SVD)", null);
        final String dir = od.getDirectory();
        final String name = od.getFileName();
        if (dir == null || name == null) return null;
        return IJ.openImage(new File(dir, name).getAbsolutePath());
    }

    private static String determineOutputDir(ImagePlus imp, String baseName, boolean interactiveUI) {
        final FileInfo fi = imp.getOriginalFileInfo();
        if (fi != null && fi.directory != null) {
            final File rootOut = new File(new File(fi.directory), baseName + "_fiba_svd_process");
            return appendRunTimestampDir(rootOut.getAbsolutePath());
        }
        if (!interactiveUI) return null;
        final OpenDialog od = new OpenDialog("Select image again to infer output location", null);
        final String d = od.getDirectory();
        if (d == null) return null;
        return appendRunTimestampDir(new File(d, baseName + "_fiba_svd_process").getAbsolutePath());
    }

    private static String appendRunTimestampDir(String rootOutputDir) {
        if (rootOutputDir == null || rootOutputDir.trim().isEmpty()) return rootOutputDir;
        final String stamp = new SimpleDateFormat("yyyyMMdd_HHmmss").format(new Date());
        return new File(rootOutputDir, stamp).getAbsolutePath();
    }

    private static String stripExtension(String name) {
        if (name == null) return "image";
        final int dot = name.lastIndexOf('.');
        if (dot > 0) return name.substring(0, dot);
        return name;
    }
}
