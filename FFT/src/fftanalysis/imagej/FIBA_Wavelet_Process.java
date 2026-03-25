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

import java.awt.Color;
import java.awt.GraphicsEnvironment;
import java.io.File;
import java.text.SimpleDateFormat;
import java.util.Date;

/**
 * Wavelet-style local frequency analysis plugin (Morlet/Gabor-like oriented bank).
 *
 * Behavior mirrors the tile montage UX:
 *  - Works with no pre-opened image (file picker)
 *  - Writes outputs beside input image in timestamped run folders
 *  - Displays key output images in interactive ImageJ
 */
public class FIBA_Wavelet_Process implements PlugInFilter {

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

        final WaveletResult wr = computeWaveletEnergy(gray);
        final ImagePlus energyImp = asGrayImage(baseName + "_wavelet_energy", wr.energyNorm);
        final ImagePlus angleImp = asAngleColorImage(baseName + "_wavelet_angle", wr.energyNorm, wr.bestAngleDeg);
        final ImagePlus histImp = buildOrientationHistogram(baseName + "_wavelet_orientation_hist", wr.bestAngleDeg, wr.energyNorm);
        final ImagePlus panelImp = buildPanel(baseName + "_wavelet_panel", gray, wr.energyNorm, wr.bestAngleDeg);

        new FileSaver(energyImp).saveAsJpeg(new File(outDir, baseName + "_wavelet_energy.jpg").getAbsolutePath());
        new FileSaver(angleImp).saveAsJpeg(new File(outDir, baseName + "_wavelet_angle.jpg").getAbsolutePath());
        new FileSaver(histImp).saveAsJpeg(new File(outDir, baseName + "_wavelet_orientation_hist.jpg").getAbsolutePath());
        new FileSaver(panelImp).saveAsJpeg(new File(outDir, baseName + "_wavelet_panel.jpg").getAbsolutePath());

        IJ.log("[WAVELET] Writing outputs to: " + outDir.getAbsolutePath());

        if (interactiveUI) {
            energyImp.show();
            angleImp.show();
            histImp.show();
            panelImp.show();
            IJ.showStatus("Wavelet process done");
        }
    }

    private static final class WaveletResult {
        final double[][] energyNorm;
        final double[][] bestAngleDeg;

        WaveletResult(double[][] energyNorm, double[][] bestAngleDeg) {
            this.energyNorm = energyNorm;
            this.bestAngleDeg = bestAngleDeg;
        }
    }

    private static WaveletResult computeWaveletEnergy(ImageProcessor gray) {
        final int w = gray.getWidth();
        final int h = gray.getHeight();

        final double[][] img = toDouble(gray);
        final double[][] bestEnergy = new double[h][w];
        final double[][] bestAngle = new double[h][w];

        final int kernelSize = 21;
        final double sigma = 3.0;
        final double wavelength = 8.0;
        final double gamma = 0.6;

        for (int angle = 0; angle < 180; angle += 10) {
            final double[][] kernel = makeGaborKernel(kernelSize, sigma, wavelength, gamma, angle);
            final double[][] resp = convolveAbs(img, kernel);
            for (int y = 0; y < h; y++) {
                for (int x = 0; x < w; x++) {
                    if (resp[y][x] > bestEnergy[y][x]) {
                        bestEnergy[y][x] = resp[y][x];
                        bestAngle[y][x] = angle;
                    }
                }
            }
        }

        normalizeInPlace(bestEnergy);
        return new WaveletResult(bestEnergy, bestAngle);
    }

    private static double[][] makeGaborKernel(int size, double sigma, double wavelength, double gamma, double thetaDeg) {
        final double[][] k = new double[size][size];
        final int c = size / 2;
        final double th = Math.toRadians(thetaDeg);
        final double ct = Math.cos(th);
        final double st = Math.sin(th);

        double mean = 0;
        for (int y = 0; y < size; y++) {
            for (int x = 0; x < size; x++) {
                final double dx = x - c;
                final double dy = y - c;
                final double xr = dx * ct + dy * st;
                final double yr = -dx * st + dy * ct;
                final double gauss = Math.exp(-(xr * xr + (gamma * gamma) * yr * yr) / (2.0 * sigma * sigma));
                final double carrier = Math.cos(2.0 * Math.PI * xr / wavelength);
                final double v = gauss * carrier;
                k[y][x] = v;
                mean += v;
            }
        }

        mean /= (size * size);
        for (int y = 0; y < size; y++) {
            for (int x = 0; x < size; x++) {
                k[y][x] -= mean;
            }
        }
        return k;
    }

    private static double[][] convolveAbs(double[][] img, double[][] kernel) {
        final int h = img.length;
        final int w = img[0].length;
        final int kh = kernel.length;
        final int kw = kernel[0].length;
        final int cy = kh / 2;
        final int cx = kw / 2;

        final double[][] out = new double[h][w];
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                double sum = 0;
                for (int ky = 0; ky < kh; ky++) {
                    for (int kx = 0; kx < kw; kx++) {
                        int yy = y + ky - cy;
                        int xx = x + kx - cx;
                        if (yy < 0) yy = 0;
                        if (yy >= h) yy = h - 1;
                        if (xx < 0) xx = 0;
                        if (xx >= w) xx = w - 1;
                        sum += img[yy][xx] * kernel[ky][kx];
                    }
                }
                out[y][x] = Math.abs(sum);
            }
        }
        return out;
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

    private static ImagePlus asAngleColorImage(String title, double[][] energy01, double[][] angleDeg) {
        final int h = energy01.length;
        final int w = energy01[0].length;
        final ColorProcessor cp = new ColorProcessor(w, h);
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                final float hue = (float) (angleDeg[y][x] / 180.0);
                final float sat = 1.0f;
                final float bri = (float) clamp01(energy01[y][x]);
                cp.set(x, y, Color.HSBtoRGB(hue, sat, bri));
            }
        }
        return new ImagePlus(title, cp);
    }

    private static ImagePlus buildOrientationHistogram(String title, double[][] angleDeg, double[][] energy01) {
        final int bins = 18;
        final double[] x = new double[bins];
        final double[] y = new double[bins];

        for (int i = 0; i < bins; i++) {
            x[i] = i * 10.0 + 5.0;
        }

        for (int r = 0; r < angleDeg.length; r++) {
            for (int c = 0; c < angleDeg[0].length; c++) {
                int b = (int) Math.floor(angleDeg[r][c] / 10.0);
                if (b < 0) b = 0;
                if (b >= bins) b = bins - 1;
                y[b] += clamp01(energy01[r][c]);
            }
        }

        final Plot p = new Plot(title, "Angle (deg)", "Weighted count", x, y);
        p.setColor(Color.blue);
        p.add("bar", x, y);
        return p.getImagePlus();
    }

    private static ImagePlus buildPanel(String title, ImageProcessor gray, double[][] energy01, double[][] angleDeg) {
        final int w = gray.getWidth();
        final int h = gray.getHeight();
        final ColorProcessor panel = new ColorProcessor(w * 3, h);

        final ImagePlus origImp = new ImagePlus("orig", gray.duplicate().convertToByte(true));
        final ImagePlus enImp = asGrayImage("energy", energy01);
        final ImagePlus angImp = asAngleColorImage("angle", energy01, angleDeg);

        blit(panel, origImp.getProcessor().convertToRGB(), 0, 0);
        blit(panel, enImp.getProcessor().convertToRGB(), w, 0);
        blit(panel, angImp.getProcessor().convertToRGB(), 2 * w, 0);

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

    private static void normalizeInPlace(double[][] arr) {
        double max = 0;
        for (int y = 0; y < arr.length; y++) {
            for (int x = 0; x < arr[0].length; x++) {
                if (arr[y][x] > max) max = arr[y][x];
            }
        }
        if (!(max > 0)) return;
        final double inv = 1.0 / max;
        for (int y = 0; y < arr.length; y++) {
            for (int x = 0; x < arr[0].length; x++) {
                arr[y][x] *= inv;
            }
        }
    }

    private static double clamp01(double v) {
        if (v < 0) return 0;
        if (v > 1) return 1;
        return v;
    }

    private static ImagePlus promptForInputImage() {
        final OpenDialog od = new OpenDialog("Select image to analyze (wavelet)", null);
        final String dir = od.getDirectory();
        final String name = od.getFileName();
        if (dir == null || name == null) return null;
        return IJ.openImage(new File(dir, name).getAbsolutePath());
    }

    private static String determineOutputDir(ImagePlus imp, String baseName, boolean interactiveUI) {
        final FileInfo fi = imp.getOriginalFileInfo();
        if (fi != null && fi.directory != null) {
            final File rootOut = new File(new File(fi.directory), baseName + "_fiba_wavelet_process");
            return appendRunTimestampDir(rootOut.getAbsolutePath());
        }
        if (!interactiveUI) return null;
        final OpenDialog od = new OpenDialog("Select image again to infer output location", null);
        final String d = od.getDirectory();
        if (d == null) return null;
        return appendRunTimestampDir(new File(d, baseName + "_fiba_wavelet_process").getAbsolutePath());
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
