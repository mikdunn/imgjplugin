// Generate ONLY the crop boxes overlay for the specified image.
// This does NOT run any FFT on the full image.

path = "C:/Users/dunnmk/Downloads/C15D5P001.jpg";
outDir = "C:/Users/dunnmk/Downloads";

open(path);

// Run tile montage with outputs mostly off; keep overlay on.
// trackX=true makes crop boxes follow the specimen across X.
run("FIBA Tile Montage (MATLAB)",
    "outputDir="+outDir+" overwrite=false " +
    "saveOverlay=true saveMontage=false savePlot=false savePerTile=false saveSol=false saveCsv=false " +
    "labelOverlay=true overlayLineWidth=2 tilesY=7 distributeYToEdges=true trackX=true");

close();
