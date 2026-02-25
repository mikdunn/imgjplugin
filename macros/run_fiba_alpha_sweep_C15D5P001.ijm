// Alpha sweep on a single real image to visualize edge-artifact suppression by the Tukey window.
// Uses ONLY this file (per instruction):
//   C:\Users\dunnmk\Downloads\C15D5P001.jpg
//
// Saves outputs into the SAME folder (Downloads) without creating a new folder.
// Each alpha value gets a unique title so *_rec.jpg and *_dat.jpg do not overwrite.

path = "C:/Users/dunnmk/Downloads/C15D5P001.jpg";
outDir = "C:/Users/dunnmk/Downloads";

open(path);

// Alpha values to sweep.
// IMPORTANT: Use string tags that avoid floating-point formatting weirdness
// (e.g., 0.30000000000000004), which can break window selection + file naming.
alphaTags = newArray("0p0", "0p1", "0p2", "0p3");
alphas    = newArray(0.0,   0.1,   0.2,   0.3);

// Disable artifact-suppression features so we can see the pure Tukey (alpha) effect on edge peaks.
// NOTE: This macro runs the TILE workflow (cropboxes) so FFT/polar/mask/iFFT are done PER TILE.
baseOpts = "outputDir="+outDir+" overwrite=false tilesY=7 trackX=true ";
baseOpts = baseOpts + "suppressAngleSpike=false removeFullHeightVerticalLine=false removeFullWidthHorizontalLine=false ";

for (i=0; i<lengthOf(alphas); i++) {
  a = alphas[i];
  tag = alphaTags[i];

  // Duplicate so we can set a unique, stable title per alpha.
  // Duplicate becomes the active image, so no selectWindow() needed.
  run("Duplicate...", "title=C15D5P001_alpha_"+tag);

  // Run the cropbox/tile pipeline
  run("FIBA Tile Montage (MATLAB)", baseOpts + "alpha="+a);

  // Close the duplicate
  close();
}

// Close the original
selectWindow("C15D5P001.jpg");
close();
