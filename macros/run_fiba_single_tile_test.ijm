// Run the single-tile (already-cropped) FIBA pipeline on one image.
// Input (selected by user):
//   C:\Users\dunnmk\Downloads\C15D5P001_alpha_0p0_tile1_crop.jpg
// Outputs written into Downloads as:
//   tile_test_rec.jpg
//   tile_test_dat.jpg
// And a debug log in Downloads:
//   fiba_plugin_debug.txt

path = "C:/Users/dunnmk/Downloads/C15D5P001_alpha_0p0_tile1_crop.jpg";
outDir = "C:/Users/dunnmk/Downloads";

open(path);
rename("tile_test");

opts = "outputDir="+outDir+" save=true showComposite=false showPlot=false debug=true ";
// For MATLAB-faithful comparison (no extra suppression beyond fiba.m):
opts = opts + "suppressAngleSpike=false removeFullHeightVerticalLine=false removeFullWidthHorizontalLine=false ";
// Keep the canonical params unless you want to tweak:
opts = opts + "alpha=0.4 beta=0.3 gamma=0.3 rmin=4";

run("FIBA Orientation (MATLAB)", opts);

close();
