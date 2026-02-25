// Generate polar graphs across alpha values and save a montage.
// Input: C:/Users/dunnmk/Downloads/C15D5P001_tiles/4.jpg
// Output folder: C:/Users/dunnmk/Downloads/tile4_polar_alpha_0_to_4/
// Produces: tile4_a_*_polar.jpg (plus *_rec.jpg and *_dat.jpg per run)

setBatchMode(true);

path = "C:/Users/dunnmk/Downloads/C15D5P001_tiles/4.jpg";
outDir = "C:/Users/dunnmk/Downloads/tile4_polar_alpha_0_to_4";
File.makeDirectory(outDir);

// Sweep alpha from 0.0 to 0.4 in steps of 0.1 (5 plots)
for (a = 0.0; a <= 0.4001; a += 0.1) {
    open(path);

    aStr = d2s(a, 1); // fixed 1 decimal
    aName = replace(aStr, ".", "p");
    rename("tile4_a_" + aName);

    opts = "outputDir=" + outDir + " save=true saveSolCsv=false showComposite=false showPlot=false debug=false ";
    // Make alpha differences visible in polar outputs (fixed scaling + optional Tukey on plotted SOL)
    opts = opts + "plotApplyTukey=true polarScaleMax=0.03 ";
    // Do NOT enable any hard-coded artifact suppression here.
    // Artifacts must be handled ONLY via alpha values.
    opts = opts + "alpha=" + aStr + " beta=0.3 gamma=0.3 rmin=4";

    run("FIBA Orientation (MATLAB)", opts);

    close();
}

setBatchMode(false);
