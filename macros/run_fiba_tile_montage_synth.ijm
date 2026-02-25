// Synthetic regression macro for the tile-only montage workflow.
// Creates a mostly-black image with a bright content block, then runs the plugin headlessly.

outDir = "C:/Users/dunnmk/fiba_out/tile_montage_test";
File.makeDirectory(outDir);

// Create a black image and draw some content.
newImage("tile_montage_synth", "8-bit black", 256, 256, 1);
setForegroundColor(255, 255, 255);
makeRectangle(60, 20, 120, 220);
run("Fill");

// Add some diagonal structure within the content box.
setLineWidth(2);
for (i=0; i<12; i++) {
  x1 = 60;
  y1 = 25 + i*18;
  x2 = 180;
  y2 = y1 + 35;
  makeLine(x1, y1, x2, y2);
  run("Draw");
}

// Run the tile montage plugin.
run("FIBA Tile Montage (MATLAB)", "outputDir="+outDir+" tilesY=6 threshold=0 colFrac=0.1 rowFrac=0.1 overwrite=true saveOverlay=true saveMontage=true savePlot=true savePerTile=true saveCsv=true wrap90=true debug=true");

// Cleanup (optional)
close();
