// Single-tile FIBA run on a specified crop tile.
// Input (requested):
//   C:\Users\dunnmk\Downloads\C15D5P001_alpha_0p3_tile4_crop.jpg
// Outputs written into Downloads as:
//   tile4_test_rec.jpg
//   tile4_test_dat.jpg
//   tile4_test_sol.csv

path = "C:/Users/dunnmk/Downloads/C15D5P001_alpha_0p3_tile4_crop.jpg";
outDir = "C:/Users/dunnmk/Downloads";

open(path);
rename("tile4_test");

opts = "outputDir="+outDir+" save=true saveSolCsv=true showComposite=false showPlot=false debug=true ";
// MATLAB-faithful comparison (no extra suppression beyond fiba.m):
opts = opts + "suppressAngleSpike=false removeFullHeightVerticalLine=false removeFullWidthHorizontalLine=false ";
opts = opts + "alpha=0.3 beta=0.3 gamma=0.3 rmin=4";

run("FIBA Orientation (MATLAB)", opts);

close();
