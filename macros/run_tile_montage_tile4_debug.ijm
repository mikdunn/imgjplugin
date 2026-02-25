// Headless-friendly macro to run tile montage on a single test image.
// Outputs will be written to Downloads unless outputDir is overridden.

setBatchMode(true);

// Update this path if your test input is elsewhere.
input = "C:/Users/dunnmk/Downloads/C15D5P001_tiles/4.jpg";

open(input);

// Force output dir explicitly for troubleshooting consistency.
outDir = "C:/Users/dunnmk/Downloads";

// Note: option parsing is whitespace-delimited key=value.
run("FIBA Tile Montage (MATLAB)",
    "outputDir=" + outDir +
    " overwrite=true" +
    " savePerTile=true" +
    " saveSol=true" +
    " saveCsv=true" +
    " polarAxes=true" +
    " saveFftTif=true");

close();
run("Quit");
