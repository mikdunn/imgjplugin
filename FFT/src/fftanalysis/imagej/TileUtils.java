package fftanalysis.imagej;

import java.util.ArrayList;
import java.util.List;

/**
 * Helpers for tiling a mostly-black image into overlapping square crops.
 *
 * The intent is to mimic the spirit of the MATLAB cropimage.m workflow, but
 * without manual ROI selection:
 *  - detect the "content" region by non-black pixels
 *  - choose a square tile width from the content's left/right edges
 *  - create vertically-overlapping tiles that reach the bottom of the image
 */
public final class TileUtils {

    private TileUtils() {
    }

    public static final class ContentBox {
        /** inclusive */
        public final int left;
        /** inclusive */
        public final int right;
        /** inclusive */
        public final int top;

        public ContentBox(int left, int right, int top) {
            this.left = left;
            this.right = right;
            this.top = top;
        }

        public int width() {
            return right - left + 1;
        }
    }

    /**
     * Find left/right content bounds by looking for columns with "enough" non-black pixels,
     * and a top bound by looking for rows with "enough" non-black pixels.
     *
     * @param pixels     8-bit pixels, row-major, length = width*height
     * @param threshold  pixel values > threshold are treated as non-black
     * @param colFrac    fraction of max column count required to be considered content (0..1]
     * @param rowFrac    fraction of max row count required to be considered content (0..1]
     */
    public static ContentBox findContentBox(byte[] pixels, int width, int height,
                                           int threshold,
                                           double colFrac,
                                           double rowFrac) {
        if (pixels == null) throw new IllegalArgumentException("pixels is null");
        if (width <= 0 || height <= 0) throw new IllegalArgumentException("invalid dimensions");
        if (pixels.length < width * height) throw new IllegalArgumentException("pixels length mismatch");

        if (!(colFrac > 0)) colFrac = 0.1;
        if (colFrac > 1) colFrac = 1;
        if (!(rowFrac > 0)) rowFrac = 0.1;
        if (rowFrac > 1) rowFrac = 1;

        final int[] colCounts = new int[width];
        final int[] rowCounts = new int[height];

        for (int y = 0; y < height; y++) {
            final int rowOff = y * width;
            for (int x = 0; x < width; x++) {
                final int v = pixels[rowOff + x] & 0xFF;
                if (v > threshold) {
                    colCounts[x]++;
                    rowCounts[y]++;
                }
            }
        }

        int maxCol = 0;
        for (int x = 0; x < width; x++) if (colCounts[x] > maxCol) maxCol = colCounts[x];
        int maxRow = 0;
        for (int y = 0; y < height; y++) if (rowCounts[y] > maxRow) maxRow = rowCounts[y];

        if (maxCol == 0 || maxRow == 0) {
            // all black (or below threshold)
            return new ContentBox(0, width - 1, 0);
        }

        final int colMin = Math.max(1, (int) Math.ceil(maxCol * colFrac));
        final int rowMin = Math.max(1, (int) Math.ceil(maxRow * rowFrac));

        int left = 0;
        while (left < width && colCounts[left] < colMin) left++;
        int right = width - 1;
        while (right >= 0 && colCounts[right] < colMin) right--;

        if (left >= width || right < 0 || right < left) {
            left = 0;
            right = width - 1;
        }

        int top = 0;
        while (top < height && rowCounts[top] < rowMin) top++;
        if (top >= height) top = 0;

        return new ContentBox(left, right, top);
    }

    /**
     * Make a size even (as in cropimage.m) by subtracting 1 if odd.
     */
    public static int forceEvenSize(int size) {
        if (size <= 0) return size;
        return (size % 2 == 0) ? size : (size - 1);
    }

    public static int clamp(int v, int lo, int hi) {
        if (v < lo) return lo;
        if (v > hi) return hi;
        return v;
    }

    /**
     * Compute tile top positions from a content-derived startY down to the bottom edge.
     *
     * The number of tiles is chosen minimally (like ceil(height/tileSize)) while still
     * reaching the bottom, using evenly-spaced overlap.
     *
     * @param startY      suggested top of first tile (e.g., top of content)
     * @param imageHeight image height in pixels
     * @param tileSize    square tile size (pixels)
     */
    public static int[] computeTileTops(int startY, int imageHeight, int tileSize) {
        if (imageHeight <= 0) throw new IllegalArgumentException("imageHeight must be > 0");
        if (tileSize <= 0) throw new IllegalArgumentException("tileSize must be > 0");
        if (tileSize > imageHeight) throw new IllegalArgumentException("tileSize exceeds image height");

        final int endTop = imageHeight - tileSize;

        // Clamp start so all tiles can fit.
        startY = clamp(startY, 0, endTop);

        final int span = imageHeight - startY; // pixels from startY to bottom
        int count = (int) Math.ceil(span / (double) tileSize);
        if (count < 1) count = 1;

        // If only one tile, align it to bottom so it "reaches the bottom".
        if (count == 1) {
            return new int[]{endTop};
        }

        final int[] tops = new int[count];
        for (int i = 0; i < count; i++) {
            double t = (double) i / (double) (count - 1);
            int y = startY + (int) Math.round(t * (endTop - startY));
            tops[i] = y;
        }

        // Ensure exact endpoints (avoid rounding drift)
        tops[0] = startY;
        tops[count - 1] = endTop;

        // Ensure non-decreasing and within bounds.
        int prev = tops[0];
        for (int i = 1; i < tops.length; i++) {
            if (tops[i] < prev) tops[i] = prev;
            if (tops[i] > endTop) tops[i] = endTop;
            prev = tops[i];
        }

        // Remove accidental duplicates at the end (rare, but can happen with tiny spans).
        // Keep at least 1.
        List<Integer> uniq = new ArrayList<Integer>();
        int last = Integer.MIN_VALUE;
        for (int i = 0; i < tops.length; i++) {
            if (tops[i] != last) {
                uniq.add(tops[i]);
                last = tops[i];
            }
        }
        int[] out = new int[uniq.size()];
        for (int i = 0; i < uniq.size(); i++) out[i] = uniq.get(i);
        return out;
    }
}
