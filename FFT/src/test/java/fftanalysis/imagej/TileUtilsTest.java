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
}
