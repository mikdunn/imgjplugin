package fftanalysis.imagej;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

public class FibaMatlabProcessorTest {

    @Test
    void process_returnsExpectedShapesAndRanges() {
        int n = 128;
        double[][] img = new double[n][n];
        for (int y = 0; y < n; y++) {
            for (int x = 0; x < n; x++) {
                img[y][x] = (x + y) % 256;
            }
        }

        FibaMatlabProcessor.Params p = new FibaMatlabProcessor.Params();
        p.rmin = 4;
        p.rmax = (n / 2) - 1;
        p.alpha = 0.4;
        p.beta = 0.3;
        p.gamma = 0.3;

        FibaMatlabProcessor.Result r = FibaMatlabProcessor.process(img, p);

        assertEquals(n, r.n);
        assertEquals(n / 2, r.w);
        assertNotNull(r.sol);
        assertEquals(180, r.sol.length);
        assertNotNull(r.imgS);
        assertEquals(n, r.imgS.length);
        assertEquals(n, r.imgS[0].length);
        assertNotNull(r.imgFDisp);
        assertEquals(n, r.imgFDisp.length);
        assertEquals(n, r.imgFDisp[0].length);
        assertNotNull(r.imgR2);
        assertEquals(n, r.imgR2.length);
        assertEquals(n, r.imgR2[0].length);
        assertNotNull(r.overlay);
        assertEquals(n, r.overlay.length);
        assertEquals(n, r.overlay[0].length);
        assertEquals(3, r.overlay[0][0].length);

        // Ranges
        for (int i = 0; i < r.sol.length; i++) {
            assertTrue(r.sol[i] >= 0.0, "SOL must be non-negative");
        }
        assertTrue(r.pAng >= 0 && r.pAng < 180);
    }

    @Test
    void process_detectsOrientedSinusoidWithinTolerance() {
        // Create a synthetic oriented grating. The FFT magnitude has peaks
        // orthogonal to the grating direction; the MATLAB algorithm's pAng
        // should be consistent for this pattern.
        int n = 256;
        double[][] img = new double[n][n];

        // Desired fiber direction ~30 degrees.
        // Note: the MATLAB-modeled implementation measures angles from the vertical axis
        // (row direction): 0deg = vertical, 90deg = horizontal.
        double thetaDeg = 30.0;
        double theta = Math.toRadians(thetaDeg);
        // Construct a grating whose wave-vector is at thetaDeg from vertical.
        // In image coords (x=col, y=row), this means:
        //   x component = sin(theta), y component = cos(theta).
        double fx = Math.sin(theta);
        double fy = Math.cos(theta);
        double freq = 8.0 / n; // cycles per pixel

        for (int y = 0; y < n; y++) {
            for (int x = 0; x < n; x++) {
                double t = (x * fx + y * fy) * (2.0 * Math.PI * freq);
                double v = 128.0 + 80.0 * Math.sin(t);
                img[y][x] = v;
            }
        }

        FibaMatlabProcessor.Params p = new FibaMatlabProcessor.Params();
        p.rmin = 4;
        p.rmax = (n / 2) - 1;
        p.alpha = 0.4;
        p.beta = 0.3;
        p.gamma = 0.3;

        FibaMatlabProcessor.Result r = FibaMatlabProcessor.process(img, p);

        // The FFT magnitude has peaks orthogonal to the grating (fiber) direction.
        // The MATLAB-modeled pipeline reports the fiber direction, so we expect
        // ~thetaDeg + 90 (mod 180). We allow a generous tolerance and accept the
        // 180-degree symmetry.
        int pang = r.pAng;
        double expectedFiberDeg = thetaDeg + 90.0;
        int diff = angleDiffDeg(pang, (int) Math.round(expectedFiberDeg));
        int diffSym = angleDiffDeg(pang, (int) Math.round(expectedFiberDeg + 180.0));
        int best = Math.min(diff, diffSym);
        assertTrue(best <= 20, "Expected pAng near " + expectedFiberDeg + "deg, got " + pang);
    }

    private static int angleDiffDeg(int a, int b) {
        int d = Math.abs(a - b) % 360;
        if (d > 180) d = 360 - d;
        return d;
    }
}
