package fftanalysis.imagej;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

public class TileUtilsTest {

    @Test
    void findContentBox_detectsSimpleCenteredRectangle() {
        int w = 20;
        int h = 30;
        byte[] px = new byte[w * h];

        // Draw a bright rectangle from x=5..14, y=7..25
        for (int y = 7; y <= 25; y++) {
            for (int x = 5; x <= 14; x++) {
                px[y * w + x] = (byte) 200;
            }
        }

        TileUtils.ContentBox box = TileUtils.findContentBox(px, w, h, 0, 0.5, 0.5);
        assertEquals(5, box.left);
        assertEquals(14, box.right);
        assertEquals(7, box.top);
    }

    @Test
    void computeTileTops_reachesBottomAndIsMonotone() {
        int imageHeight = 100;
        int tileSize = 20;
        int startY = 10;

        int[] tops = TileUtils.computeTileTops(startY, imageHeight, tileSize);

        assertTrue(tops.length >= 1);
        assertEquals(startY, tops[0]);
        assertEquals(imageHeight - tileSize, tops[tops.length - 1]);
        for (int i = 1; i < tops.length; i++) {
            assertTrue(tops[i] >= tops[i - 1]);
            assertTrue(tops[i] <= imageHeight - tileSize);
        }
    }

    @Test
    void forceEvenSize_matchesMatlabBehavior() {
        assertEquals(100, TileUtils.forceEvenSize(101));
        assertEquals(100, TileUtils.forceEvenSize(100));
        assertEquals(0, TileUtils.forceEvenSize(0));
    }

    @Test
    void findBandEdgesByGradient_detectsJumpDropEdges() {
        int w = 20;
        int h = 10;
        byte[] px = new byte[w * h];

        // Background ~10 everywhere.
        for (int i = 0; i < px.length; i++) px[i] = (byte) 10;

        // Specimen band from x=5..14 with strong jump/drop.
        for (int y = 0; y < h; y++) {
            for (int x = 5; x <= 14; x++) {
                px[y * w + x] = (byte) 200;
            }
        }
        // Add a bright internal spike to ensure we don't lock onto an internal gradient.
        for (int y = 0; y < h; y++) {
            px[y * w + 9] = (byte) 255;
        }

        int[] e = TileUtils.findBandEdgesByGradient(px, w, h, 0, h, 2);
        assertTrue(Math.abs(e[0] - 5) <= 1, "left edge should be near 5 but was " + e[0]);
        assertTrue(Math.abs(e[1] - 14) <= 1, "right edge should be near 14 but was " + e[1]);
        assertTrue(Math.abs(e[2] - ((5 + 14) / 2)) <= 1, "center should be near 9 but was " + e[2]);
    }
}
