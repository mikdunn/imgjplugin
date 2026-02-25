// Synthetic regression macro for the tile-only montage workflow (variant output dir).

outDir = "C:/Users/dunnmk/fiba_out/tile_montage_test2";
File.makeDirectory(outDir);

newImage("tile_montage_synth", "8-bit black", 256, 256, 1);
setForegroundColor(255, 255, 255);
makeRectangle(60, 20, 120, 220);
run("Fill");

setLineWidth(2);
for (i=0; i<12; i++) {
  x1 = 60;
  y1 = 25 + i*18;
  x2 = 180;
  y2 = y1 + 35;
  makeLine(x1, y1, x2, y2);
  run("Draw");
}

run("FIBA Tile Montage (MATLAB)", "outputDir="+outDir+" tilesY=6 threshold=0 colFrac=0.1 rowFrac=0.1 overwrite=true saveOverlay=true saveMontage=true savePlot=true savePerTile=true saveSol=true saveCsv=true wrap90=true debug=true");

close();
