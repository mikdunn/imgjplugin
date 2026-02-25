package fftanalysis.imagej;

import org.jtransforms.fft.DoubleFFT_2D;

import java.util.ArrayList;
import java.util.List;

/**
 * MATLAB-faithful implementation of the core logic from fiba.m:
 * <ul>
 *   <li>imadjust + Tukey window on the input image</li>
 *   <li>2D FFT + fftshift magnitude</li>
 *   <li>orientation signal extraction SOL(theta)</li>
 *   <li>statistically significant peak detection + weighted median peak angle</li>
 *   <li>band-limited inverse FFT reconstruction with blended mask edges</li>
 * </ul>
 *
 * This class intentionally uses plain double arrays so it can be used from
 * ImageJ UI code, batch processing, tests, etc.
 */
public final class FibaMatlabProcessor {

    private FibaMatlabProcessor() {
    }

    public static final class Params {
        /** inner radius of usable frequency band */
        public int rmin = 4;
        /** outer radius of usable frequency band; if <= 0, defaults to w-1 */
        public int rmax = -1;

        /** size of Tukey window on original sample */
        public double alpha = 0.4;
        /** size of Tukey window on FFT result (r-direction) */
        public double beta = 0.3;
        /** size of Tukey window on FFT result (theta-direction) */
        public double gamma = 0.3;

        /** expected image is square; if not, caller should crop first */
        public boolean requireSquare = true;

        /**
         * Optional artifact suppression: remove any perfectly-vertical line that spans (nearly) the full image height
         * in the reconstructed mask (imgR2) after thresholding.
         *
         * Motivation: some workflows produce a non-biological, axis-aligned line artifact that can run from the
         * very top to the very bottom of an image/tile, which should not be considered part of the mask.
         */
        public boolean removeFullHeightVerticalLine = false;
        /** Minimum fraction of rows that must be nonzero in a column for that column to be removed (1.0 = strict full height). */
        public double removeFullHeightVerticalLineMinCoverage = 1.0;

        /**
         * Optional artifact suppression: remove any perfectly-horizontal line that spans (nearly) the full image width
         * in the reconstructed mask (imgR2) after thresholding.
         *
         * Motivation: some pipelines yield axis-aligned banding artifacts that can appear as full-width horizontal lines.
         */
        public boolean removeFullWidthHorizontalLine = false;
        /** Minimum fraction of columns that must be nonzero in a row for that row to be removed (1.0 = strict full width). */
        public double removeFullWidthHorizontalLineMinCoverage = 1.0;

        /**
         * Optional artifact suppression: if SOL has a strong spike at an exact angle (e.g. 90deg),
         * attenuate that bin (and optionally a small neighborhood) BEFORE peak finding and mask creation.
         *
         * Motivation: some imaging workflows introduce non-biological, perfectly axial line artifacts
         * that manifest as an unnaturally sharp peak at exactly 0/90 degrees.
         */
        public boolean suppressAngleSpike = false;
        /** center angle (deg in [0,179]) to suppress; typical: 90 */
        public int suppressAngleDeg = 90;
        /** Optional additional angles (deg in [0,179]) to suppress. If set, all angles in this list are processed. */
        public int[] suppressAnglesDeg = null;
        /** half-width in degrees around suppressAngleDeg to attenuate; 0 means only the exact bin */
        public int suppressHalfWidthDeg = 0;
        /** Only suppress if SOL[angle] is greater than this ratio times median(SOL). */
        public double suppressIfOverMedianRatio = 6.0;
    }

    public static final class Result {
        /** N = 2*w (image side length) */
        public int n;
        public int w;

        /** orientation signal, length 180, angles 0..179 */
        public double[] sol;
        public double meanSol;
        public double stdSol;

        /** Weighted-average peak fiber angle (0..179), matching fiba.m */
        public int pAng;
        /** 1 if there might be more than one peak */
        public int warnPk;

        /** statistically significant peak width (deg) */
        public int spWid;
        /** peak strength (0..1), i.e., BandH in MATLAB */
        public double bandStrength;
        /** 30% peak bandwidth in degrees */
        public int pWidth;

        /** peak boundary angles in [0,179] for reconstruction; may wrap */
        public int ang1;
        public int ang2;

        /** indices (angles) used for highlighting reconstructed band on plot */
        public int[] aind1;
        public int[] aind2;

        /** display panels (all scaled to 0..1) */
        public double[][] origNorm;
        public double[][] imgS;
        public double[][] imgFDisp;
        /** the frequency-domain reconstruction mask used for inverse FFT (0..1-ish) */
        public double[][] reconMask;
        public double[][] imgR2;
        public double[][][] overlay;
    }

    /**
     * Run the full algorithm on a square grayscale image.
     *
     * @param j input image intensities (N x N)
     */
    public static Result process(double[][] j, Params params) {
        if (j == null || j.length == 0 || j[0].length == 0) {
            throw new IllegalArgumentException("Input image is empty");
        }
        final int n = j.length;
        final int m = j[0].length;
        if (params.requireSquare && n != m) {
            throw new IllegalArgumentException("Input image must be square (got " + n + "x" + m + ")");
        }

        final Result out = new Result();
        out.n = n;
        out.w = n / 2;

        final int w = out.w;
        final int rmin = Math.max(0, params.rmin);
        final int rmax = (params.rmax > 0) ? params.rmax : (w - 1);
        if (rmax <= rmin) {
            throw new IllegalArgumentException("rmax must be > rmin");
        }

        // ======================= Step I: FFT2 Analysis =====================
        // MATLAB: J1 = imadjust(J)
        final double[][] j1 = imadjustStretch01(j, 0.01);

        final double avgJ = mean(j1);

        // MATLAB: generate edge blending (Tukey window). Equivalent to tukeywin(n, alpha)
        // and applying a separable 2D window S = s*s'.
        final double[] s = tukeyWin1D(n, params.alpha);
        final double[][] imgS = new double[n][n];
        double maxImgS = 0;
        for (int row = 0; row < n; row++) {
            for (int col = 0; col < n; col++) {
                final double S = s[row] * s[col];
                final double v = S * (j1[row][col] - avgJ) + avgJ;
                imgS[row][col] = v;
                if (v > maxImgS) maxImgS = v;
            }
        }
        if (maxImgS <= 0) maxImgS = 1;
        for (int row = 0; row < n; row++) {
            for (int col = 0; col < n; col++) {
                imgS[row][col] /= maxImgS;
            }
        }
        out.imgS = imgS;

        // MATLAB: K = fft2(ImgS)
        final DoubleFFT_2D fft2 = new DoubleFFT_2D(n, n);
        final double[][] K = new double[n][2 * n];
        for (int row = 0; row < n; row++) {
            for (int col = 0; col < n; col++) {
                // JTransforms realForwardFull expects the real input in a[row][0..n-1].
                // The full complex output is written back interleaved into a[row][0..2n-1].
                K[row][col] = imgS[row][col];
            }
        }
        fft2.realForwardFull(K);

        final double[][] amp = new double[n][n];
        final double[][] phase = new double[n][n];
        for (int row = 0; row < n; row++) {
            for (int col = 0; col < n; col++) {
                final double re = K[row][2 * col];
                final double im = K[row][2 * col + 1];
                amp[row][col] = Math.hypot(re, im);
                phase[row][col] = Math.atan2(im, re);
            }
        }
        // MATLAB: ImgF = fftshift(abs(K))
        final double[][] imgF = fftShift(amp);

        // Display spectrum (visualization only): robust log-power scaling.
        // Using a simple max() normalization often makes axis-aligned components look "clipped"
        // when a few pixels (especially DC) dominate the dynamic range.
        final double[][] imgFDisp = spectrumDisplay01(imgF);
        out.imgFDisp = imgFDisp;

        // MATLAB: original normalized as double(J)/max
        out.origNorm = normalizeByMax(j);

        // ====================== Step II: Data Analysis =====================
        // Extract orientation data as in imfft() in MATLAB.
        final double[] B = extractOrientationSignal(imgF, w, rmin, rmax);
        final double sumB = sum(B);
        final double[] sol = new double[180];
        for (int t = 0; t < 90; t++) {
            sol[t] = B[t + 90] / sumB;
        }
        for (int t = 90; t < 180; t++) {
            sol[t] = B[t - 90] / sumB;
        }
        final double[] solFiltered = maybeSuppressAngleSpike(sol, params);
        out.sol = solFiltered;
        out.meanSol = mean(solFiltered);
        out.stdSol = std(solFiltered, out.meanSol);

        // Identify statistically significant peak following MATLAB logic.
        final PeakInfo peak = findPeak(solFiltered, out.meanSol, out.stdSol);
        out.pAng = peak.pAng;
        out.warnPk = peak.warnPk;
        out.spWid = peak.spWid;
        out.bandStrength = peak.bandStrength;
        out.pWidth = peak.pWidth;
        out.ang1 = peak.ang1;
        out.ang2 = peak.ang2;

        // ================= Step III: Inverse FFT2 Reconstruction ===========
        final Reconstruction recon = reconstruct(j1, imgS, imgF, phase, w, rmin, rmax,
            params.beta, params.gamma, out.ang1, out.ang2,
            params.removeFullHeightVerticalLine, params.removeFullHeightVerticalLineMinCoverage,
            params.removeFullWidthHorizontalLine, params.removeFullWidthHorizontalLineMinCoverage);
        out.reconMask = recon.mask;
        out.imgR2 = recon.imgR2;
        out.overlay = recon.overlay;
        out.aind1 = recon.aind1;
        out.aind2 = recon.aind2;

        return out;
    }

    private static double[] maybeSuppressAngleSpike(double[] sol, Params params) {
        if (sol == null || sol.length != 180) return sol;
        if (params == null || !params.suppressAngleSpike) return sol;

        final int halfWidth = Math.max(0, params.suppressHalfWidthDeg);
        final int[] centers;
        if (params.suppressAnglesDeg != null && params.suppressAnglesDeg.length > 0) {
            centers = new int[params.suppressAnglesDeg.length];
            for (int i = 0; i < centers.length; i++) centers[i] = mod180(params.suppressAnglesDeg[i]);
        } else {
            centers = new int[] { mod180(params.suppressAngleDeg) };
        }

        final double med = median(sol);
        if (!(med > 0)) {
            // If SOL is all zeros (or NaN), there's nothing sensible to do.
            return sol;
        }

        double[] out = sol;
        boolean changed = false;
        for (int c = 0; c < centers.length; c++) {
            final int center = centers[c];

            final double ratio = out[center] / med;
            if (!(ratio >= params.suppressIfOverMedianRatio)) {
                continue;
            }

            // Additional check: make sure this looks like a narrow spike (higher than immediate neighbors).
            final double n1 = out[mod180(center - 1)];
            final double n2 = out[mod180(center + 1)];
            final double neighborAvg = 0.5 * (n1 + n2);
            if (!(out[center] > neighborAvg)) {
                continue;
            }

            if (!changed) {
                out = out.clone();
                changed = true;
            }
            for (int d = -halfWidth; d <= halfWidth; d++) {
                out[mod180(center + d)] = neighborAvg;
            }
        }

        if (!changed) return sol;

        // Re-normalize to keep SOL as a probability-like distribution.
        final double s = sum(out);
        if (s > 0) {
            for (int i = 0; i < out.length; i++) out[i] /= s;
        }
        return out;
    }

    private static int mod180(int a) {
        int m = a % 180;
        if (m < 0) m += 180;
        return m;
    }

    private static double median(double[] a) {
        final double[] copy = new double[a.length];
        int n = 0;
        for (int i = 0; i < a.length; i++) {
            final double v = a[i];
            if (Double.isNaN(v) || Double.isInfinite(v)) continue;
            copy[n++] = v;
        }
        if (n == 0) return Double.NaN;
        java.util.Arrays.sort(copy, 0, n);
        final int mid = n / 2;
        if ((n & 1) == 1) return copy[mid];
        return 0.5 * (copy[mid - 1] + copy[mid]);
    }

    // ---------------------------------------------------------------------
    // Orientation extraction
    // ---------------------------------------------------------------------

    private static double[] extractOrientationSignal(double[][] imgFShifted, int w, int rmin, int rmax) {
        // MATLAB code builds A(i,j) for i=1..w-1, j=1..180 with bilinear sampling.
        // In the MATLAB reference, theta is measured from the vertical axis (row direction):
        // theta = 0 follows +row, theta = 90 follows +col.
        // This matches later reconstruction code using atan(x/y) where x is column-offset and y is row-offset.
        final double[] B = new double[180];

        final int maxR = Math.min(w - 2, rmax - 1);
        final int minR = Math.max(0, rmin);

        for (int thetaDeg = 0; thetaDeg < 180; thetaDeg++) {
            double sum = 0;
            final double theta = Math.toRadians(thetaDeg);
            for (int r = minR; r <= maxR; r++) {

                // row and col coordinates in the fftshifted image.
                final double row = r * Math.cos(theta) + w; // 0-based center
                final double col = r * Math.sin(theta) + w;

                final double v = bilinearSampleXY(imgFShifted, row, col);
                sum += v;
            }
            B[thetaDeg] = sum;
        }

        return B;
    }

    private static double bilinearSampleXY(double[][] img, double xRow, double yCol) {
        final int nRows = img.length;
        final int nCols = img[0].length;

        int x0 = (int) Math.floor(xRow);
        int x1 = (int) Math.ceil(xRow);
        int y0 = (int) Math.floor(yCol);
        int y1 = (int) Math.ceil(yCol);

        // Clamp to bounds to avoid out-of-range (MATLAB indexing implicitly assumes in-range)
        if (x0 < 0) x0 = 0;
        if (y0 < 0) y0 = 0;
        if (x1 >= nRows) x1 = nRows - 1;
        if (y1 >= nCols) y1 = nCols - 1;

        final double dx = xRow - x0;
        final double dy = yCol - y0;

        final double v00 = img[x0][y0];
        final double v10 = img[x1][y0];
        final double v01 = img[x0][y1];
        final double v11 = img[x1][y1];

        return (1 - dx) * (1 - dy) * v00
                + dx * (1 - dy) * v10
                + (1 - dx) * dy * v01
                + dx * dy * v11;
    }

    private static final class PeakInfo {
        int pAng;
        int warnPk;
        int spWid;
        double bandStrength;
        int pWidth;
        int ang1;
        int ang2;
    }

    private static PeakInfo findPeak(double[] sol, double mean, double std) {
        // MATLAB:
        // SOL3 = [SOL SOL SOL]
        // peakbd = find(SOL3 > (mSOL + stdSOL))
        // [~, peakid] = max(SOL3)
        // ind = find(peakbd == peakid) + length(peakbd)/3
        // find contiguous bounds around that position
        // index = [peakbd(peakl) peakbd(peakr)]
        // weighted-median angle within the peak band

        final double[] sol3 = new double[180 * 3];
        for (int i = 0; i < sol3.length; i++) sol3[i] = sol[i % 180];

        final double thresh = mean + std;
        final List<Integer> peakbd = new ArrayList<Integer>();
        for (int i = 0; i < sol3.length; i++) {
            if (sol3[i] > thresh) peakbd.add(i);
        }

        // If nothing exceeds threshold, fall back to global max as a 1-degree peak.
        int peakId = 0;
        double best = sol3[0];
        for (int i = 1; i < sol3.length; i++) {
            if (sol3[i] > best) {
                best = sol3[i];
                peakId = i;
            }
        }
        // Force peakId into the middle SOL3 copy (indices 180..359) to avoid wrap ambiguities.
        peakId = 180 + (peakId % 180);

        int indInPeakbd = -1;
        for (int i = 0; i < peakbd.size(); i++) {
            if (peakbd.get(i) == peakId) {
                indInPeakbd = i;
                break;
            }
        }

        // If the peak itself isn't above threshold, treat the peak as a single point.
        int leftIdx = peakId;
        int rightIdx = peakId;
        int warnPk = 0;

        if (!peakbd.isEmpty() && indInPeakbd >= 0) {
            // peakId is already in the middle copy; use its exact position in peakbd.
            final int centerPos = indInPeakbd;

            // Scan right for break
            int r = centerPos;
            while (r + 1 < peakbd.size() && peakbd.get(r + 1) - peakbd.get(r) == 1) r++;
            // Scan left for break
            int l = centerPos;
            while (l - 1 >= 0 && peakbd.get(l) - peakbd.get(l - 1) == 1) l--;

            leftIdx = peakbd.get(l);
            rightIdx = peakbd.get(r);

            // MATLAB warning heuristic
            final int peakWidth = (rightIdx - leftIdx + 1);
            if (peakWidth < peakbd.size() / 3) {
                warnPk = 1;
            }
        }

        final PeakInfo out = new PeakInfo();
        out.warnPk = warnPk;

        // Weighted-median within the band
        double bandH = 0;
        for (int i = leftIdx; i <= rightIdx; i++) bandH += sol3[i];
        if (bandH <= 0) bandH = 1;

        double bandS = 0;
        int pAng = 0;
        for (int i = leftIdx; i <= rightIdx; i++) {
            bandS += sol3[i] / bandH;
            if (bandS >= 0.5) {
                // 0-based equivalent of MATLAB's index - 180 - 1 (1-based)
                pAng = i - 180;
                break;
            }
        }

        // Peak boundary angles for reconstruction (MATLAB: ang = index - 180)
        int ang1 = leftIdx - 180;
        int ang2 = rightIdx - 180;
        final int spWid = ang2 - ang1 + 1;

        // 30% bandwidth around pAng
        // MATLAB uses SOL3(pAng-i+180 : pAng+i+180)
        int pWidth = 0;
        for (int i = 1; i <= 90; i++) {
            final int center = pAng + 180;
            final int lo = Math.max(0, center - i);
            final int hi = Math.min(sol3.length - 1, center + i);
            double band = 0;
            for (int k = lo; k <= hi; k++) band += sol3[k];
            if (band >= 0.3) {
                pWidth = 2 * (i - 1);
                break;
            }
        }

        // Fold pAng into [0,179]
        while (pAng < 0) pAng += 180;
        while (pAng >= 180) pAng -= 180;

        // Fold reconstruction bounds into [0,179] while preserving wrap (ang2 < ang1 indicates wrap).
        ang1 = mod180(ang1);
        ang2 = mod180(ang2);

        out.pAng = pAng;
        out.spWid = spWid;
        out.bandStrength = bandH;
        out.pWidth = pWidth;
        out.ang1 = ang1;
        out.ang2 = ang2;
        return out;
    }

    // ---------------------------------------------------------------------
    // Reconstruction (inverse FFT)
    // ---------------------------------------------------------------------

    private static final class Reconstruction {
        double[][] mask;
        double[][] imgR2;
        double[][][] overlay;
        int[] aind1;
        int[] aind2;
    }

        private static Reconstruction reconstruct(
            double[][] jAdjusted01,
            double[][] imgS,
            double[][] imgFShifted,
            double[][] phase,
            int w,
            int rmin,
            int rmax,
            double beta,
            double gamma,
            int ang1,
            int ang2,
            boolean removeFullHeightVerticalLine,
            double removeFullHeightVerticalLineMinCoverage,
            boolean removeFullWidthHorizontalLine,
            double removeFullWidthHorizontalLineMinCoverage
    ) {
        final int n = imgS.length;
        final int[] ainds = computeAinds(ang1, ang2);
        final int[] aind1;
        final int[] aind2;
        aind1 = extractAind1(ang1, ang2);
        aind2 = extractAind2(ang1, ang2);

        final double[][] mask = buildReconstructionMask(n, w, rmin, rmax, beta, gamma, ang1, ang2);

        // Apply mask in shifted domain, then inverse shift back.
        final double[][] maskedShifted = new double[n][n];
        for (int row = 0; row < n; row++) {
            for (int col = 0; col < n; col++) {
                maskedShifted[row][col] = imgFShifted[row][col] * mask[row][col];
            }
        }
        // Convert shifted spectrum back to the unshifted layout used by the FFT output arrays.
        // MATLAB equivalent: masked = ifftshift(maskedShifted)
        final double[][] masked = ifftShift(maskedShifted);

        // Combine masked amplitude with original phase and inverse FFT.
        final DoubleFFT_2D fft2 = new DoubleFFT_2D(n, n);
        final double[][] complex = new double[n][2 * n];
        for (int row = 0; row < n; row++) {
            for (int col = 0; col < n; col++) {
                final double a = masked[row][col];
                final double ph = phase[row][col];
                complex[row][2 * col] = a * Math.cos(ph);
                complex[row][2 * col + 1] = a * Math.sin(ph);
            }
        }
        fft2.complexInverse(complex, true);

        final double[][] imgR1 = new double[n][n];
        double min = Double.POSITIVE_INFINITY;
        for (int row = 0; row < n; row++) {
            for (int col = 0; col < n; col++) {
                final double v = complex[row][2 * col];
                imgR1[row][col] = v;
                if (v < min) min = v;
            }
        }
        // MATLAB: ImgR2 = imadjust((ImgR1-min(min(ImgR1))))
        for (int row = 0; row < n; row++) {
            for (int col = 0; col < n; col++) {
                imgR1[row][col] -= min;
            }
        }
        double[][] imgR2 = imadjustMinMax01(imgR1);

        // MATLAB: ImgR2 = ImgR2.*ImgS
        for (int row = 0; row < n; row++) {
            for (int col = 0; col < n; col++) {
                imgR2[row][col] *= imgS[row][col];
            }
        }
        // MATLAB: ImgR2(ImgR2<mean(mean(ImgR2))*1.1) = 0;
        final double meanR2 = mean(imgR2);
        final double thr = meanR2 * 1.1;
        for (int row = 0; row < n; row++) {
            for (int col = 0; col < n; col++) {
                if (imgR2[row][col] < thr) imgR2[row][col] = 0;
            }
        }

        if (removeFullHeightVerticalLine) {
            suppressFullHeightVerticalLinesInPlace(imgR2, removeFullHeightVerticalLineMinCoverage);
        }

        if (removeFullWidthHorizontalLine) {
            suppressFullWidthHorizontalLinesInPlace(imgR2, removeFullWidthHorizontalLineMinCoverage);
        }

        // Overlay: ImgS2(:,:,1) = imadjust(J)/255; ImgS2(:,:,2/3) = base - ImgR2
        final double[][][] overlay = new double[n][n][3];
        for (int row = 0; row < n; row++) {
            for (int col = 0; col < n; col++) {
                final double base = jAdjusted01[row][col];
                overlay[row][col][0] = clamp01(base);
                overlay[row][col][1] = clamp01(base - imgR2[row][col]);
                overlay[row][col][2] = clamp01(base - imgR2[row][col]);
            }
        }

        final Reconstruction out = new Reconstruction();
        out.mask = mask;
        out.imgR2 = imgR2;
        out.overlay = overlay;
        out.aind1 = aind1;
        out.aind2 = aind2;
        return out;
    }

    private static void suppressFullHeightVerticalLinesInPlace(double[][] img, double minCoverageFrac) {
        if (img == null || img.length == 0 || img[0].length == 0) return;
        final int nRows = img.length;
        final int nCols = img[0].length;

        final double frac = Math.max(0.0, Math.min(1.0, minCoverageFrac));
        final int needed = (int) Math.ceil(frac * nRows);
        if (needed <= 0) return;

        for (int col = 0; col < nCols; col++) {
            int nonZero = 0;
            for (int row = 0; row < nRows; row++) {
                if (img[row][col] > 0) nonZero++;
            }
            if (nonZero >= needed) {
                for (int row = 0; row < nRows; row++) {
                    img[row][col] = 0;
                }
            }
        }
    }

    private static void suppressFullWidthHorizontalLinesInPlace(double[][] img, double minCoverageFrac) {
        if (img == null || img.length == 0 || img[0].length == 0) return;
        final int nRows = img.length;
        final int nCols = img[0].length;

        final double frac = Math.max(0.0, Math.min(1.0, minCoverageFrac));
        final int needed = (int) Math.ceil(frac * nCols);
        if (needed <= 0) return;

        for (int row = 0; row < nRows; row++) {
            int nonZero = 0;
            for (int col = 0; col < nCols; col++) {
                if (img[row][col] > 0) nonZero++;
            }
            if (nonZero >= needed) {
                for (int col = 0; col < nCols; col++) {
                    img[row][col] = 0;
                }
            }
        }
    }

    // aind helpers (for plot highlighting)
    private static int[] computeAinds(int ang1, int ang2) {
        // kept for parity with MATLAB structure; actual arrays computed below
        return new int[0];
    }

    private static int[] extractAind1(int ang1, int ang2) {
        if (ang2 >= ang1) {
            int len = (ang2 - ang1 + 1);
            int[] a = new int[len];
            for (int i = 0; i < len; i++) a[i] = (ang1 + i) + 1; // 1-based like MATLAB
            return a;
        }
        // wrap case
        int len = (ang2 + 1);
        int[] a = new int[len];
        for (int i = 0; i < len; i++) a[i] = i + 1;
        return a;
    }

    private static int[] extractAind2(int ang1, int ang2) {
        if (ang2 >= ang1) {
            return new int[] { 0 };
        }
        int len = (179 - ang1 + 1);
        int[] a = new int[len];
        for (int i = 0; i < len; i++) a[i] = (ang1 + i) + 1;
        return a;
    }

    private static double[][] buildReconstructionMask(
            int n,
            int w,
            int rmin,
            int rmax,
            double beta,
            double gamma,
            int ang1,
            int ang2
    ) {
        // Ported directly from fiba.m / imifft(), keeping MATLAB variable names and rot90 semantics.
        int ind = 0;
        int a1 = ang1;
        int a2 = ang2;
        if (a2 < a1) {
            a2 = a2 + 180;
            ind = 1;
        }

        final int theta1 = a1 - 90;
        final int theta2 = a2 - 90;

        final double[][] maskA = new double[n][n];
        final double[][] maskB = new double[n][n];

        final double rd = (rmax - rmin + 1);

        // MATLAB:
        // for i=(w+1):(2*w)
        //   for j=1:(2*w)
        //     R = sqrt((j-w-0.5)^2 + (i-w-0.5)^2);
        //     x = (R-rmin)/Rd;
        //     theta = atan((j-w-0.5)/(i-w-0.5))/pi*180;
        //     ... build MaskA, MaskB ...
        for (int i = w; i < n; i++) {
            final double iOff = (i - w + 0.5);
            for (int j = 0; j < n; j++) {
                final double jOff = (j - w + 0.5);

                final double r = Math.sqrt(jOff * jOff + iOff * iOff);
                final double x = (r - rmin) / rd;

                final double theta = Math.toDegrees(Math.atan(jOff / iOff));

                // MaskA window in r-direction (Tukey-like)
                if (x >= 0 && x <= (beta / 2.0)) {
                    maskA[i][j] = 0.5 * (1.0 + Math.cos((2.0 * Math.PI / beta) * (x - beta / 2.0)));
                } else if (x > (beta / 2.0) && x <= (1.0 - beta / 2.0)) {
                    maskA[i][j] = 1.0;
                } else if (x > (1.0 - beta / 2.0) && x <= 1.0) {
                    maskA[i][j] = 0.5 * (1.0 + Math.cos((2.0 * Math.PI / beta) * (x - 1.0 + beta / 2.0)));
                }

                // MaskB window in theta-direction (Tukey-like)
                final double th = (theta - theta1) / (double) (theta2 - theta1);
                if (th >= 0 && th <= (gamma / 2.0)) {
                    maskB[i][j] = 0.5 * (1.0 + Math.cos((2.0 * Math.PI / gamma) * (th - gamma / 2.0)));
                } else if (th > (gamma / 2.0) && th <= (1.0 - gamma / 2.0)) {
                    maskB[i][j] = 1.0;
                } else if (th > (1.0 - gamma / 2.0) && th <= 1.0) {
                    maskB[i][j] = 0.5 * (1.0 + Math.cos((2.0 * Math.PI / gamma) * (th - 1.0 + gamma / 2.0)));
                }
            }
        }

        // MATLAB:
        // MaskA = MaskA + rot90(MaskA,2);
        // MaskB = MaskB + rot90(MaskB,2);
        final double[][] maskA2 = add(maskA, rot90(maskA, 2));
        final double[][] maskB2 = add(maskB, rot90(maskB, 2));

        // MATLAB:
        // if ind == 0
        //   Mask = MaskA.*MaskB;
        // else
        //   Mask = MaskA.*rot90(MaskB);
        // end
        if (ind == 0) {
            return multiply(maskA2, maskB2);
        }
        return multiply(maskA2, rot90(maskB2, 1));
    }

    private static double[][] rot90(double[][] in, int k) {
        // MATLAB rot90(A,k): rotates CCW by 90 degrees k times.
        int kk = k % 4;
        if (kk < 0) kk += 4;
        if (kk == 0) {
            final int n = in.length;
            final int m = in[0].length;
            final double[][] out = new double[n][m];
            for (int i = 0; i < n; i++) {
                System.arraycopy(in[i], 0, out[i], 0, m);
            }
            return out;
        }

        final int n = in.length;
        final int m = in[0].length;
        if (n != m) {
            throw new IllegalArgumentException("rot90 expects a square matrix, got " + n + "x" + m);
        }

        final double[][] out = new double[n][n];
        if (kk == 1) {
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    out[n - 1 - j][i] = in[i][j];
                }
            }
        } else if (kk == 2) {
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    out[n - 1 - i][n - 1 - j] = in[i][j];
                }
            }
        } else { // kk == 3
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    out[j][n - 1 - i] = in[i][j];
                }
            }
        }
        return out;
    }

    // ---------------------------------------------------------------------
    // Small numeric utilities (kept private to keep surface area small)
    // ---------------------------------------------------------------------

    private static double[][] normalizeByMax(double[][] in) {
        final int n = in.length;
        final int m = in[0].length;
        double max = 0;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                if (in[i][j] > max) max = in[i][j];
            }
        }
        if (max <= 0) max = 1;
        final double[][] out = new double[n][m];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                out[i][j] = in[i][j] / max;
            }
        }
        return out;
    }

    private static double mean(double[][] a) {
        double s = 0;
        long c = 0;
        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < a[0].length; j++) {
                s += a[i][j];
                c++;
            }
        }
        return (c == 0) ? 0 : s / c;
    }

    private static double mean(double[] a) {
        double s = 0;
        for (int i = 0; i < a.length; i++) s += a[i];
        return (a.length == 0) ? 0 : s / a.length;
    }

    private static double std(double[] a, double mean) {
        double s2 = 0;
        for (int i = 0; i < a.length; i++) {
            final double d = a[i] - mean;
            s2 += d * d;
        }
        return (a.length <= 1) ? 0 : Math.sqrt(s2 / (a.length - 1));
    }

    private static double sum(double[] a) {
        double s = 0;
        for (int i = 0; i < a.length; i++) s += a[i];
        return s;
    }

    private static double sum(double[][] a) {
        double s = 0;
        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < a[0].length; j++) s += a[i][j];
        }
        return s;
    }

    private static double clamp01(double v) {
        if (v < 0) return 0;
        if (v > 1) return 1;
        return v;
    }

    /**
     * MATLAB tukeywin(n, alpha) equivalent (0<=alpha<=1), as a 1D window.
     *
     * Notes:
     * - alpha=0 => rectangular window (all ones)
     * - alpha=1 => Hann window
     */
    private static double[] tukeyWin1D(int n, double alpha) {
        if (n <= 0) return new double[0];
        if (n == 1) return new double[] { 1.0 };

        final double a = Math.max(0.0, Math.min(1.0, alpha));
        final double[] w = new double[n];
        if (a <= 0.0) {
            for (int i = 0; i < n; i++) w[i] = 1.0;
            return w;
        }

        final double N = (double) (n - 1);
        final double a2 = a / 2.0;
        for (int i = 0; i < n; i++) {
            final double t = i / N; // 0..1
            final double v;
            if (t < a2) {
                v = 0.5 * (1.0 + Math.cos(Math.PI * ((2.0 * t / a) - 1.0)));
            } else if (t <= (1.0 - a2)) {
                v = 1.0;
            } else {
                v = 0.5 * (1.0 + Math.cos(Math.PI * ((2.0 * t / a) - (2.0 / a) + 1.0)));
            }
            w[i] = v;
        }
        return w;
    }

    /**
     * Approximation of MATLAB imadjust(I) for grayscale images: stretch contrast
     * to [0,1] with a small saturation fraction.
     */
    private static double[][] imadjustStretch01(double[][] in, double saturateFraction) {
        final int n = in.length;
        final int m = in[0].length;
        // Copy and find min/max first.
        double min = Double.POSITIVE_INFINITY;
        double max = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                final double v = in[i][j];
                if (v < min) min = v;
                if (v > max) max = v;
            }
        }
        if (!(max > min)) {
            // constant image
            final double[][] out = new double[n][m];
            return out;
        }

        // A simple (and fast) robust stretch: clamp to percentile-like bounds
        // computed from a coarse histogram.
        final int bins = 1024;
        final long[] hist = new long[bins];
        final double scale = (bins - 1) / (max - min);
        long total = 0;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                int b = (int) ((in[i][j] - min) * scale);
                if (b < 0) b = 0;
                if (b >= bins) b = bins - 1;
                hist[b]++;
                total++;
            }
        }

        final long cut = (long) Math.floor(total * saturateFraction);
        long c = 0;
        int loBin = 0;
        for (int b = 0; b < bins; b++) {
            c += hist[b];
            if (c > cut) {
                loBin = b;
                break;
            }
        }
        c = 0;
        int hiBin = bins - 1;
        for (int b = bins - 1; b >= 0; b--) {
            c += hist[b];
            if (c > cut) {
                hiBin = b;
                break;
            }
        }
        if (hiBin <= loBin) {
            loBin = 0;
            hiBin = bins - 1;
        }

        final double lo = min + loBin / scale;
        final double hi = min + hiBin / scale;
        final double denom = (hi - lo);

        final double[][] out = new double[n][m];
        if (!(denom > 0)) {
            return out;
        }
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                double v = (in[i][j] - lo) / denom;
                if (v < 0) v = 0;
                if (v > 1) v = 1;
                out[i][j] = v;
            }
        }
        return out;
    }

    private static double[][] imadjustMinMax01(double[][] in) {
        final int n = in.length;
        final int m = in[0].length;
        double min = Double.POSITIVE_INFINITY;
        double max = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                final double v = in[i][j];
                if (v < min) min = v;
                if (v > max) max = v;
            }
        }
        if (!(max > min)) {
            return new double[n][m];
        }
        final double denom = max - min;
        final double[][] out = new double[n][m];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                out[i][j] = (in[i][j] - min) / denom;
            }
        }
        return out;
    }

    private static double[][] fftShift(double[][] in) {
        final int n = in.length;
        final int m = in[0].length;
        final double[][] out = new double[n][m];
        final int hn = n / 2;
        final int hm = m / 2;

        for (int i = 0; i < n; i++) {
            final int ii = (i + hn) % n;
            for (int j = 0; j < m; j++) {
                final int jj = (j + hm) % m;
                out[ii][jj] = in[i][j];
            }
        }
        return out;
    }

    /**
     * Inverse of {@link #fftShift(double[][])}. MATLAB equivalent: ifftshift().
     *
     * For even sizes, fftshift and ifftshift are the same; for odd sizes they differ by 1 sample.
     */
    private static double[][] ifftShift(double[][] in) {
        final int n = in.length;
        final int m = in[0].length;
        final double[][] out = new double[n][m];
        final int hn = (n + 1) / 2;
        final int hm = (m + 1) / 2;

        for (int i = 0; i < n; i++) {
            final int ii = (i + hn) % n;
            for (int j = 0; j < m; j++) {
                final int jj = (j + hm) % m;
                out[ii][jj] = in[i][j];
            }
        }
        return out;
    }

    private static double[][] rotate180(double[][] in) {
        final int n = in.length;
        final int m = in[0].length;
        final double[][] out = new double[n][m];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                out[n - 1 - i][m - 1 - j] = in[i][j];
            }
        }
        return out;
    }

    private static double[][] rotate90CCW(double[][] in) {
        final int n = in.length;
        final int m = in[0].length;
        final double[][] out = new double[m][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                out[m - 1 - j][i] = in[i][j];
            }
        }
        return out;
    }

    /**
     * Build a stable 0..1 visualization of an fftshifted magnitude spectrum.
     *
     * Strategy:
     * - use log(1 + power) to compress dynamic range
     * - normalize by a high percentile (not max) to avoid a single spike dominating
     * - suppress the DC center pixel for display
     */
    private static double[][] spectrumDisplay01(double[][] imgFShiftedMag) {
        final int n = imgFShiftedMag.length;
        final int m = imgFShiftedMag[0].length;
        final double[][] v = new double[n][m];

        double min = Double.POSITIVE_INFINITY;
        double max = Double.NEGATIVE_INFINITY;

        // Log-power transform
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                final double a = imgFShiftedMag[i][j];
                final double p = a * a;
                final double lv = Math.log1p(p);
                v[i][j] = lv;
                if (lv < min) min = lv;
                if (lv > max) max = lv;
            }
        }

        // Suppress DC for display so it doesn't dominate the contrast.
        final int ci = n / 2;
        final int cj = m / 2;
        if (ci >= 0 && ci < n && cj >= 0 && cj < m) {
            v[ci][cj] = min;
        }

        if (!(max > min)) {
            return new double[n][m];
        }

        // Percentile-based contrast stretch (histogram) to avoid "clipping" effects.
        final int bins = 2048;
        final long[] hist = new long[bins];
        final double scale = (bins - 1) / (max - min);
        long total = 0;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                int b = (int) ((v[i][j] - min) * scale);
                if (b < 0) b = 0;
                if (b >= bins) b = bins - 1;
                hist[b]++;
                total++;
            }
        }

        // Use low/high percentiles as black/white points.
        final double loFrac = 0.01;   // 1%
        final double hiFrac = 0.999; // 99.9%
        final long loCount = (long) Math.floor(total * loFrac);
        final long hiCount = (long) Math.floor(total * hiFrac);

        long c = 0;
        int loBin = 0;
        for (int b = 0; b < bins; b++) {
            c += hist[b];
            if (c >= loCount) {
                loBin = b;
                break;
            }
        }

        c = 0;
        int hiBin = bins - 1;
        for (int b = 0; b < bins; b++) {
            c += hist[b];
            if (c >= hiCount) {
                hiBin = b;
                break;
            }
        }

        if (hiBin <= loBin) {
            loBin = 0;
            hiBin = bins - 1;
        }

        final double lo = min + loBin / scale;
        final double hi = min + hiBin / scale;
        final double denom = (hi - lo);
        final double[][] out = new double[n][m];
        if (!(denom > 0)) {
            return out;
        }
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                double t = (v[i][j] - lo) / denom;
                if (t < 0) t = 0;
                if (t > 1) t = 1;
                out[i][j] = t;
            }
        }
        return out;
    }

    private static double[][] add(double[][] a, double[][] b) {
        final int n = a.length;
        final int m = a[0].length;
        final double[][] out = new double[n][m];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                out[i][j] = a[i][j] + b[i][j];
            }
        }
        return out;
    }

    private static double[][] multiply(double[][] a, double[][] b) {
        final int n = a.length;
        final int m = a[0].length;
        final double[][] out = new double[n][m];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                out[i][j] = a[i][j] * b[i][j];
            }
        }
        return out;
    }
}
