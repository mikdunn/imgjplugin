// Run the TILE montage workflow once (single alpha) and save per-tile FFT panels
// so we can inspect the power spectrum image (tile_*_fft.jpg) for axis clipping.
//
// Input:
//   C:/Users/dunnmk/Downloads/C15D5P001.jpg
// Output:
//   C:/Users/dunnmk/Downloads/C15D5P001_tilemontage_tukeyfix_a0p3/
//
// Notes:
// - This uses the cropbox/tile pipeline, so FFT/mask/iFFT are done ONLY on cropped tiles.
// - No hard-coded artifact suppression is enabled; artifacts are handled ONLY via alpha.

setBatchMode(true);

path = "C:/Users/dunnmk/Downloads/C15D5P001.jpg";
outDir = "C:/Users/dunnmk/Downloads/C15D5P001_tilemontage_tukeyfix_a0p3";
File.makeDirectory(outDir);

open(path);

// Duplicate to avoid touching the original window title.
run("Duplicate...", "title=C15D5P001_tilemontage_a0p3");

opts = "outputDir="+outDir+" overwrite=true tilesY=7 trackX=true ";
opts = opts + "saveOverlay=true saveMontage=true savePlot=true savePerTilePanels=true saveSolPlots=true saveCsv=true ";
opts = opts + "alpha=0.3 beta=0.3 gamma=0.3 rmin=4";

run("FIBA Tile Montage (MATLAB)", opts);

// Close windows
close();
selectWindow("C15D5P001.jpg");
close();

setBatchMode(false);
